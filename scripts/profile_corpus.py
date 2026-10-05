"""Corpus profile: per-organism sizes at each reduction stage of the extraction
pipeline. Written once to corpus_profile.json and consumed by report_corpus_sweep.

Stages
------
  GEM         reactions in the source gem_models/*.json (metadata.json)
  cached      the pre-built FASTCC-reduced model the sweep loads
              (extraction/fastcc/models/, built at eps=1e-6 with growth preserved)
  consistent  what extract_bioproduction_pathway's OWN in-pipeline FASTCC keeps
              when re-run on the cached model at FASTCC_EPS = 1e-4

The cached/consistent gap is the reason a target drawn as "producible" from the
cached model can still raise TargetNotFound: the stricter in-pipeline eps drops
reactions, remove_orphans=True then deletes the now-orphaned target metabolite.

Run: .venv/Scripts/python.exe -u -m scripts.profile_corpus
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.extraction.fastcc import available_models, fastcc_from_cobra
from src.extraction.pathway import FASTCC_EPS
from scripts.validate_gudie_all_organisms import load_reduced

OUT = PROJECT_ROOT / "corpus_profile.json"


def main():
    meta = json.loads((PROJECT_ROOT / "gem_models" / "metadata.json").read_text(encoding="utf-8"))
    profile = {}
    for key in sorted(available_models()):
        rec = {"organism": meta.get(key, {}).get("organism", ""),
               "gem_reactions": meta.get(key, {}).get("num_reactions"),
               "gem_metabolites": meta.get(key, {}).get("num_metabolites"),
               "gem_genes": meta.get(key, {}).get("num_genes")}
        try:
            m = load_reduced(key)
            rec["cached_reactions"] = len(m.reactions)
            rec["cached_metabolites"] = len(m.metabolites)
            _, cons = fastcc_from_cobra(m, FASTCC_EPS)
            rec["consistent_reactions"] = len(set(cons))
            rec["consistent_frac"] = round(len(set(cons)) / max(1, len(m.reactions)), 4)
        except Exception as e:
            rec["error"] = f"{type(e).__name__}: {str(e)[:80]}"
        profile[key] = rec
        print(f"{key:10} {rec.get('cached_reactions','-'):>6} cached -> "
              f"{rec.get('consistent_reactions','-'):>6} consistent "
              f"({rec.get('consistent_frac','-')})", flush=True)
    OUT.write_text(json.dumps(profile, indent=1, sort_keys=True), encoding="utf-8")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
