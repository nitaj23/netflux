"""
identity.py - What makes two submitted GEMs "the same model".

Three identifiers, three jobs (§3.1.1):

``upload_hash``   sha256 of the raw bytes. Provenance only.
``file_hash``     sha256 of a canonical serialization. Exact-duplicate key
                  *within* a format -- immune to key order, float spelling and
                  annotation noise, but NOT stable across SBML <-> JSON.
``structure_fp``  semantic identity. Stable across formats, because it is built
                  only from what survives the round trip: the stoichiometric
                  matrix, metabolite compartments, the gene set, and the
                  canonicalized GPR of each reaction.

The tiering is deliberate (§3.2.1). Identity = stoichiometry, compartments,
genes, reaction->GPR. Configuration = bounds and objective (``config_fp``),
because a re-parameterized copy is a *version*, not a duplicate. Annotation =
names, subsystems, formulas, charges -- never compared, and exactly the fields
measured as unstable across a format change.

GPR canonicalization is syntactic, not semantic: ``b1 and (b2 or b3)`` and
``(b1 and b2) or (b1 and b3)`` hash differently. Canonical DNF would close that
gap and can blow up exponentially on real genome-scale rules. The failure
direction is a false *negative* -- a re-spelled rule reads as a different model
and gets reviewed, rather than being auto-rejected as a duplicate.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
import unicodedata
from pathlib import Path

import cobra
from cobra.core.gene import GPR
from cobra.io.dict import model_to_dict

from src.schema import (Gene, Metabolite, ModelGene, ModelMetabolite,
                                ModelReaction, Reaction, ReactionMetabolite)

#: Bump on any change to the canonical form; stored per submission so old and
#: new hashes are never silently compared.
CANONICALIZER_VERSION = "canon-1"
FINGERPRINT_VERSION = "fp-2"

#: Dropped from the canonical form: serializer-dependent or pure provenance.
_VOLATILE = {"notes", "annotation", "version"}

_UNIT = "\x1f"          # field separator; no BiGG-style id contains it
_GZIP_MAGIC = b"\x1f\x8b"


# ── primitives ───────────────────────────────────────────────────────────────
def _norm(value) -> str:
    """NFC-normalize an id so a decomposed twin hashes identically."""
    return unicodedata.normalize("NFC", str(value))


def _sortkey(value: str) -> bytes:
    """Sort by UTF-8 bytes, never by SQL collation: SQLite's BINARY and
    PostgreSQL's locale collation disagree on the mixed-case ids this corpus
    actually contains (EXCH_cellb_e, R_Biomass, BETA_D_FRUCTOSE)."""
    return value.encode("utf-8")


def _num(value) -> str:
    """12 significant digits, with -0.0 folded to 0.0.

    Coefficients round-trip through two serializers as doubles; 12 digits
    absorbs representation noise while keeping 1.0 and 1.0000001 distinct.
    """
    return repr(float(f"{float(value):.12g}") + 0.0)


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


# ── GPR canonicalization ─────────────────────────────────────────────────────
def _tree(node):
    """cobra GPR ast -> nested tuple, flattened, deduped and sorted."""
    if isinstance(node, ast.Name):
        return ("leaf", node.id)
    if isinstance(node, ast.BoolOp):
        op = "and" if isinstance(node.op, ast.And) else "or"
        kids: list = []
        for value in node.values:
            child = _tree(value)
            if child[0] == op:                       # a and (b and c) -> a, b, c
                kids.extend(child[1])
            else:
                kids.append(child)
        unique = {_render(k): k for k in kids}       # b1 and b1 -> b1
        if len(unique) == 1:
            # Collapse the single survivor. Without this "b1 and b1" renders
            # "(b1)" while "b1" renders "b1", and one rule gets two hashes.
            return next(iter(unique.values()))
        return (op, tuple(unique[k] for k in sorted(unique, key=_sortkey)))
    raise ValueError(f"unsupported GPR node {type(node).__name__}")


def _render(node) -> str:
    if node[0] == "leaf":
        return node[1]
    glue = "&" if node[0] == "and" else "|"
    return "(" + glue.join(_render(k) for k in node[1]) + ")"


def canon_gpr(rule: str | None) -> str:
    """Canonical form of a gene-reaction rule; "" for an empty rule.

    An unparseable rule hashes to a value derived from its text -- never to "",
    so a broken rule can never read as "no GPR".
    """
    text = (rule or "").strip()
    if not text:
        return ""
    try:
        body = GPR.from_string(text).body
    except Exception:
        body = None
    if body is None:
        # cobra does not raise on a malformed rule -- it warns and hands back an
        # empty GPR. Returning "" here would make a broken rule indistinguishable
        # from a reaction that genuinely has none, so a non-empty rule that did
        # not parse gets its own value.
        return "UNPARSEABLE:" + hashlib.sha256(text.encode()).hexdigest()[:16]
    return _render(_tree(body))


# ── upload handling ──────────────────────────────────────────────────────────
def sniff_format(raw: bytes) -> str:
    """Detect the real format from content, never the filename (FR-1b)."""
    if raw[:2] == _GZIP_MAGIC:
        return "gzip"
    head = raw[:4096].lstrip()
    if head[:1] in (b"{", b"["):
        return "json"
    if head[:1] == b"<" or b"<sbml" in raw[:4096]:
        return "sbml"
    raise ValueError("not recognizable as cobra-JSON or SBML")


def decompress(raw: bytes, max_bytes: int) -> bytes:
    """Bounded gunzip -- a cap here is the zip-bomb guard (C-4)."""
    with gzip.GzipFile(fileobj=__import__("io").BytesIO(raw)) as handle:
        out = handle.read(max_bytes + 1)
    if len(out) > max_bytes:
        raise ValueError(f"decompressed upload exceeds {max_bytes} bytes")
    return out


def load_model(path: Path, fmt: str) -> cobra.Model:
    """Parse an on-disk upload. SBML goes through libSBML; JSON through cobra."""
    if fmt == "sbml":
        model = cobra.io.read_sbml_model(str(path))
        _restore_subsystems(model)
        return model
    return cobra.io.load_json_model(str(path))


def _restore_subsystems(model: cobra.Model) -> None:
    """Put subsystems back from SBML groups (FR-1d).

    cobra usually does this on read; when it has not, the group name is the
    subsystem. Measured: an SBML written without groups loses every subsystem,
    and ``reactions.subsystem`` is an indexed catalogue column the atlas
    colours by -- so V15 reports what is left after this repair.
    """
    for group in getattr(model, "groups", []) or []:
        label = group.name or group.id
        for member in group.members:
            if isinstance(member, cobra.Reaction) and not member.subsystem:
                member.subsystem = label


# ── canonical serialization (file_hash) ──────────────────────────────────────
def _strip(obj):
    if isinstance(obj, dict):
        return {k: _strip(v) for k, v in obj.items() if k not in _VOLATILE}
    if isinstance(obj, list):
        return [_strip(v) for v in obj]
    return obj


def canonical_bytes(model: cobra.Model) -> bytes:
    """A defined normalization of the model, not raw ``to_json`` output.

    ``to_json`` serializes in parser order with full float repr, so two loads
    of the same file can differ. This form sorts every collection, normalizes
    every float, canonicalizes every GPR and drops volatile keys.
    """
    doc = _strip(model_to_dict(model))
    for key in ("reactions", "metabolites", "genes"):
        doc[key] = sorted(doc.get(key, []), key=lambda e: _sortkey(_norm(e["id"])))
    for rxn in doc["reactions"]:
        rxn["id"] = _norm(rxn["id"])
        mets = rxn.get("metabolites") or {}
        rxn["metabolites"] = {
            _norm(m): _num(c)
            for m, c in sorted(mets.items(), key=lambda kv: _sortkey(_norm(kv[0])))
        }
        for field in ("lower_bound", "upper_bound", "objective_coefficient"):
            if rxn.get(field) is not None:
                rxn[field] = _num(rxn[field])
        rxn["gene_reaction_rule"] = canon_gpr(rxn.get("gene_reaction_rule"))
    for met in doc["metabolites"]:
        met["id"] = _norm(met["id"])
        if met.get("charge") is not None:
            met["charge"] = int(met["charge"])
    doc["canonicalizer"] = CANONICALIZER_VERSION
    return json.dumps(doc, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")


def file_hash(model: cobra.Model) -> str:
    return sha256_bytes(canonical_bytes(model))


# ── fingerprints ─────────────────────────────────────────────────────────────
def _digest(reactions: list[str], metabolites: list[str], genes: list[str]) -> str:
    payload = "\n".join([FINGERPRINT_VERSION,
                         "R:", *reactions, "M:", *metabolites, "G:", *genes])
    return sha256_bytes(payload.encode("utf-8"))


def structure_fingerprint(model: cobra.Model) -> str:
    """Semantic identity: S-matrix + compartments + genes + reaction->GPR.

    Genes appear twice on purpose -- once as a set (a gene added or dropped)
    and once inside each rule (the same gene set re-mapped to different
    reactions, which changes every knockout prediction).
    """
    reactions = []
    for rxn in sorted(model.reactions, key=lambda r: _sortkey(_norm(r.id))):
        column = ",".join(
            f"{_norm(met.id)}:{_num(coeff)}"
            for met, coeff in sorted(rxn.metabolites.items(),
                                     key=lambda kv: _sortkey(_norm(kv[0].id))))
        reactions.append(_UNIT.join((_norm(rxn.id), column,
                                     canon_gpr(rxn.gene_reaction_rule))))
    metabolites = [
        f"{_norm(met.id)}{_UNIT}{_norm(met.compartment or '')}"
        for met in sorted(model.metabolites, key=lambda m: _sortkey(_norm(m.id)))]
    genes = sorted((_norm(g.id) for g in model.genes), key=_sortkey)
    return _digest(reactions, metabolites, genes)


def config_fingerprint(model: cobra.Model) -> str:
    """Environment: bounds, objective coefficients, objective direction."""
    bounds = [f"{_norm(r.id)}{_UNIT}{_num(r.lower_bound)}{_UNIT}{_num(r.upper_bound)}"
              for r in sorted(model.reactions, key=lambda r: _sortkey(_norm(r.id)))]
    objective = [f"{_norm(r.id)}{_UNIT}{_num(r.objective_coefficient)}"
                 for r in sorted(model.reactions, key=lambda r: _sortkey(_norm(r.id)))
                 if r.objective_coefficient]
    payload = "\n".join([FINGERPRINT_VERSION, "B:", *bounds, "O:", *objective,
                         "D:" + str(model.objective.direction)])
    return sha256_bytes(payload.encode("utf-8"))


def db_structure_fingerprint(db, model_id: int) -> str:
    """Rebuild ``structure_fingerprint`` from the catalogue -- no file parsing.

    Uses the *same* helpers as the in-memory path, because if the two paths
    differ by one byte the duplicate check silently never fires.

    One caveat (C-7): ``reactions.gene_rule`` is global and first-writer-wins,
    so for a reaction shared between models it reflects whichever model
    registered it first. A miss here is therefore re-checked against the
    candidate's own file before declaring "not a duplicate" (see service).
    """
    rxn_rows = (db.query(Reaction.id, Reaction.bigg_id, Reaction.gene_rule)
                  .join(ModelReaction, ModelReaction.reaction_id == Reaction.id)
                  .filter(ModelReaction.model_id == model_id).all())
    met_rows = (db.query(Metabolite.id, Metabolite.bigg_id, Metabolite.compartment)
                  .join(ModelMetabolite, ModelMetabolite.metabolite_id == Metabolite.id)
                  .filter(ModelMetabolite.model_id == model_id).all())
    gene_rows = (db.query(Gene.bigg_id)
                   .join(ModelGene, ModelGene.gene_id == Gene.id)
                   .filter(ModelGene.model_id == model_id).all())
    stoich_rows = (db.query(ReactionMetabolite.reaction_id,
                            ReactionMetabolite.metabolite_id,
                            ReactionMetabolite.stoichiometry)
                     .filter(ReactionMetabolite.model_id == model_id).all())

    met_name = {mid: bigg for mid, bigg, _ in met_rows}
    columns: dict[int, list[tuple[str, float]]] = {}
    for rxn_id, met_id, coeff in stoich_rows:
        columns.setdefault(rxn_id, []).append((met_name.get(met_id, ""), coeff))

    reactions = []
    for rxn_id, bigg, rule in rxn_rows:
        column = ",".join(
            f"{_norm(m)}:{_num(c)}"
            for m, c in sorted(columns.get(rxn_id, []), key=lambda kv: _sortkey(_norm(kv[0]))))
        reactions.append(_UNIT.join((_norm(bigg), column, canon_gpr(rule))))
    reactions.sort(key=lambda s: _sortkey(s.split(_UNIT)[0]))

    metabolites = sorted(
        (f"{_norm(bigg)}{_UNIT}{_norm(comp or '')}" for _, bigg, comp in met_rows),
        key=_sortkey)
    genes = sorted((_norm(g[0]) for g in gene_rows), key=_sortkey)
    return _digest(reactions, metabolites, genes)
