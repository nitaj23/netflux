import { useQuery, useQueryClient } from '@tanstack/react-query'
import { useMemo, useRef, useState } from 'react'
import { ApiError } from '../api/client'
import { compoundOrganisms, fetchDefaults, searchCompounds, uploadGem } from '../api/endpoints'
import type { CompoundGroup, CompoundOut, ModelDefaults, Overrides, SourceSpec } from '../api/types'
import { useDebounced } from '../lib/hooks'
import type { RunSnapshot } from '../state'
import { Combobox, type Option } from './Combobox'

const UNIT = 'mmol/gDW/h'

interface Gem {
  token: string
  name: string
  warnings: string[]
  summary: Record<string, unknown>
}

/** The server returns substring matches alphabetically; surface exact and prefix matches first. */
function rankMatches(list: CompoundOut[], text: string, preferProducible: boolean): CompoundOut[] {
  const t = text.toLowerCase()
  const rank = (c: CompoundOut) => {
    const names = [c.base_id, c.name ?? ''].map((x) => x.toLowerCase())
    return names.some((n) => n === t) ? 0 : names.some((n) => n.startsWith(t)) ? 1 : 2
  }
  const producible = (c: CompoundOut) => (preferProducible && !c.suggested_target ? 1 : 0)
  return [...list].sort((a, b) => producible(a) - producible(b) || rank(a) - rank(b))
}

/** Fold entries that share a KEGG id into one group; the richest entry (most compartments) names it. */
function groupByKegg(list: CompoundOut[]): CompoundGroup[] {
  const out: CompoundGroup[] = []
  const byKegg = new Map<string, CompoundGroup>()
  for (const c of list) {
    const g = c.kegg_id ? byKegg.get(c.kegg_id) : undefined
    if (!g) {
      const fresh: CompoundGroup = { ...c, base_ids: [c.base_id] }
      out.push(fresh)
      if (c.kegg_id) byKegg.set(c.kegg_id, fresh)
    } else {
      g.base_ids.push(c.base_id)
      if (c.compartments.length > g.compartments.length) {
        g.base_id = c.base_id
        g.name = c.name ?? g.name
        g.formula = c.formula ?? g.formula
        g.compartments = c.compartments
      }
    }
  }
  return out
}

const describe = (e: unknown) => (e instanceof ApiError ? e : new ApiError(0, { detail: String(e) }))

export function ParamPanel({ onGenerate }: { onGenerate: (s: RunSnapshot) => void }) {
  const qc = useQueryClient()
  const fileInput = useRef<HTMLInputElement>(null)
  const seq = useRef(0)

  const [gem, setGem] = useState<Gem | null>(null)
  const [uploading, setUploading] = useState(false)
  const [uploadError, setUploadError] = useState<ApiError | null>(null)

  const [compoundText, setCompoundText] = useState('')
  const [compound, setCompound] = useState<CompoundGroup | null>(null)
  const [modelId, setModelId] = useState<number | null>(null)
  const [defaults, setDefaults] = useState<ModelDefaults | null>(null)
  const [target, setTarget] = useState<string | null>(null)
  const [defaultsLoading, setDefaultsLoading] = useState(false)
  const [defaultsError, setDefaultsError] = useState<ApiError | null>(null)

  const [substrate, setSubstrate] = useState('')
  const [lb, setLb] = useState('')
  const [ub, setUb] = useState('')
  const [objective, setObjective] = useState('')
  const [formError, setFormError] = useState<string | null>(null)

  const gemToken = gem?.token ?? null
  const debounced = useDebounced(compoundText.trim(), 250)

  const compounds = useQuery({
    queryKey: ['compounds', gemToken, debounced],
    queryFn: ({ signal }) => searchCompounds(debounced, gemToken, signal),
    enabled: debounced.length > 0 && !compound,
    staleTime: 60_000,
  })

  const organisms = useQuery({
    queryKey: ['organisms', compound?.base_ids],
    queryFn: () => compoundOrganisms(compound!.base_ids),
    enabled: !!compound && !gemToken,
    staleTime: 60_000,
  })

  const groupedMatches = useMemo(() => {
    const ranked = rankMatches(compounds.data ?? [], debounced, !!gemToken)
    return gemToken ? ranked.map((c) => ({ ...c, base_ids: [c.base_id] })) : groupByKegg(ranked)
  }, [compounds.data, debounced, gemToken])

  const compoundOptions: Option[] = compound
    ? []
    : groupedMatches.map((c) => ({
        value: c.base_id,
        label: c.name ?? c.base_id,
        sub: [[...new Set(c.base_ids)].join(' / '), c.formula, gemToken && !c.suggested_target ? 'no exchange' : null]
          .filter(Boolean)
          .join(' · '),
      }))

  const substrateOptions: Option[] = useMemo(() => {
    if (!defaults) return []
    const uptake = Object.entries(defaults.exchange_bounds)
      .filter(([, [lo]]) => lo < 0)
      .map(([id]) => id)
    const ids = [...new Set([...defaults.carbon_source_options, ...uptake.sort()])]
    const t = substrate.trim().toLowerCase()
    return ids
      .filter((id) => id.toLowerCase().includes(t))
      .slice(0, 50)
      .map((id) => ({ value: id, label: id }))
  }, [defaults, substrate])

  const objectiveOptions: Option[] = (defaults?.objective_rxns ?? []).map((id) => ({ value: id, label: id }))

  function fillBounds(d: ModelDefaults, exchangeId: string | null | undefined) {
    const b = exchangeId ? d.exchange_bounds[exchangeId] : undefined
    setLb(b ? String(b[0]) : '')
    setUb(b ? String(b[1]) : '')
  }

  function resetSelection() {
    seq.current++
    setCompoundText('')
    setCompound(null)
    setModelId(null)
    setDefaults(null)
    setTarget(null)
    setDefaultsLoading(false)
    setDefaultsError(null)
    setSubstrate('')
    setLb('')
    setUb('')
    setObjective('')
    setFormError(null)
  }

  async function loadDefaults(baseId: string, mid: number | null, token: string | null, fallbackTarget?: string | null) {
    const n = ++seq.current
    setDefaultsLoading(true)
    setDefaultsError(null)
    try {
      const d = await qc.fetchQuery({
        queryKey: ['defaults', token ?? mid, baseId],
        queryFn: () => fetchDefaults(baseId, mid, token),
        staleTime: 60_000,
      })
      if (n !== seq.current) return
      setDefaults(d)
      setTarget(d.target ?? fallbackTarget ?? null)
      setObjective(d.objective_rxn ?? '')
      setSubstrate('')
      fillBounds(d, d.carbon_source)
    } catch (e) {
      if (n === seq.current) setDefaultsError(describe(e))
    } finally {
      if (n === seq.current) setDefaultsLoading(false)
    }
  }

  function chooseCompound(c: CompoundGroup) {
    setCompound(c)
    setCompoundText(c.name ?? c.base_id)
    setModelId(null)
    setDefaults(null)
    setTarget(null)
    setSubstrate('')
    setLb('')
    setUb('')
    setObjective('')
    setFormError(null)
    if (gemToken) void loadDefaults(c.base_id, null, gemToken, c.suggested_target)
  }

  function chooseOrganism(id: number) {
    setModelId(id)
    setFormError(null)
    if (!compound) return
    const row = organisms.data?.find((o) => o.model.id === id)
    // Each model stores the compound under its own id scheme; use the row's id.
    void loadDefaults(row?.base_id ?? compound.base_id, id, null, row?.suggested_target)
  }

  async function onFile(file: File | undefined) {
    if (!file) return
    setUploading(true)
    setUploadError(null)
    try {
      const r = await uploadGem(file)
      resetSelection()
      setGem({ token: r.gem_token, name: file.name, warnings: r.warnings, summary: r.summary })
    } catch (e) {
      setUploadError(describe(e))
    } finally {
      setUploading(false)
      if (fileInput.current) fileInput.current.value = ''
    }
  }

  function removeGem() {
    setGem(null)
    setUploadError(null)
    resetSelection()
  }

  function generate() {
    const fail = (m: string) => setFormError(m)
    if (!compound) return fail('Choose a compound first.')
    if (!gem && modelId == null) return fail('Choose an organism.')
    if (!defaults) return fail('Model defaults have not loaded yet.')
    if (!target) return fail('This model has no exchange reaction for that compound, so it cannot be produced here.')
    const obj = objective.trim()
    if (!obj) return fail('Set an objective reaction (normally the biomass reaction).')

    const overrides: Overrides = { objective_rxn: obj }
    if (lb.trim() !== '' || ub.trim() !== '') {
      const lo = Number(lb)
      const hi = Number(ub)
      const exId = substrate.trim() || defaults.carbon_source
      if (lb.trim() === '' || ub.trim() === '' || !Number.isFinite(lo) || !Number.isFinite(hi))
        return fail('Flux LB and UB must both be numbers.')
      if (lo > hi) return fail('Flux LB cannot be greater than UB.')
      if (!exId) return fail('Enter a substrate so the flux bounds know which exchange they apply to.')
      if (!(exId in defaults.exchange_bounds)) return fail(`"${exId}" is not an exchange reaction in this model.`)
      overrides.custom_exchange_bounds = { [exId]: [lo, hi] }
    }

    const source: SourceSpec = gem
      ? { strategy: 'custom', gem_token: gem.token, overrides }
      : { strategy: 'repo', model_id: modelId!, overrides }
    const organism = gem ? gem.name : (organisms.data?.find((o) => o.model.id === modelId)?.model.organism ?? '')
    setFormError(null)
    onGenerate({
      source,
      target,
      biomass: obj,
      substrate: substrate.trim() || undefined,
      carbonSource: substrate.trim() || defaults.carbon_source || undefined,
      label: `${organism} → ${compound.name ?? compound.base_id}`,
    })
  }

  const boundsExchange = substrate.trim() || defaults?.carbon_source
  const nReactions = gem?.summary.n_reactions ?? gem?.summary.reactions

  return (
    <form
      className="params"
      onSubmit={(e) => {
        e.preventDefault()
        generate()
      }}
    >
      <div className="row">
        <label htmlFor="compound">Compound</label>
        <Combobox
          id="compound"
          value={compoundText}
          loading={compounds.isFetching}
          options={compoundOptions}
          onChange={(t) => {
            setCompoundText(t)
            if (compound) {
              setCompound(null)
              setModelId(null)
              setDefaults(null)
              setTarget(null)
            }
          }}
          onSelect={(o) => {
            const c = groupedMatches.find((x) => x.base_id === o.value)
            if (c) chooseCompound(c)
          }}
        />
      </div>

      <div className="row">
        <label htmlFor="organism">Organism</label>
        <select
          id="organism"
          className="field"
          disabled={!!gem || !compound || organisms.isLoading}
          value={gem ? 'gem' : (modelId ?? '')}
          onChange={(e) => chooseOrganism(Number(e.target.value))}
        >
          {gem ? (
            <option value="gem">Uploaded: {gem.name}</option>
          ) : (
            <>
              <option value="" disabled>
                {!compound ? 'Choose a compound first' : organisms.isLoading ? 'Loading…' : 'Select organism'}
              </option>
              {organisms.data?.map((o) => (
                <option key={o.model.id} value={o.model.id}>
                  {o.model.organism}
                  {!o.has_exchange ? ' (no exchange)' : !o.carries_flow ? ' (no flow)' : ''}
                </option>
              ))}
            </>
          )}
        </select>
      </div>
      {organisms.isError && <p className="note err">{describe(organisms.error).message}</p>}
      {compounds.isError && <p className="note err">{describe(compounds.error).message}</p>}

      <div className="row">
        <label htmlFor="substrate">Substrate</label>
        <Combobox
          id="substrate"
          value={substrate}
          disabled={!defaults}
          placeholder={defaults?.carbon_source ?? 'auto-detected'}
          options={substrateOptions}
          onChange={(t) => {
            setSubstrate(t)
            if (!defaults) return
            const id = t.trim()
            if (id === '') fillBounds(defaults, defaults.carbon_source)
            else if (id in defaults.exchange_bounds) fillBounds(defaults, id)
          }}
          onSelect={(o) => {
            setSubstrate(o.value)
            if (defaults) fillBounds(defaults, o.value)
          }}
        />
      </div>

      <div className="row flux">
        <span className="label">Flux</span>
        <div className="flux-grid">
          <label htmlFor="lb">LB</label>
          <input id="lb" className="field" inputMode="decimal" value={lb} disabled={!defaults} onChange={(e) => setLb(e.target.value)} />
          <input className="field unit" value={UNIT} readOnly tabIndex={-1} aria-label="LB unit" />
          <label htmlFor="ub">UB</label>
          <input id="ub" className="field" inputMode="decimal" value={ub} disabled={!defaults} onChange={(e) => setUb(e.target.value)} />
          <input className="field unit" value={UNIT} readOnly tabIndex={-1} aria-label="UB unit" />
        </div>
      </div>
      {boundsExchange && defaults && <p className="note">Bounds apply to {boundsExchange}</p>}

      <div className="row">
        <label htmlFor="objective">Objective</label>
        <Combobox
          id="objective"
          value={objective}
          disabled={!defaults}
          options={objectiveOptions}
          onChange={setObjective}
          onSelect={(o) => setObjective(o.value)}
        />
      </div>

      {defaultsLoading && <p className="note">Loading model defaults…</p>}
      {defaultsError && <p className="note err">{defaultsError.message}</p>}
      {defaults && !target && !defaultsLoading && (
        <p className="note err">This model has no exchange reaction for that compound.</p>
      )}

      {gem && (
        <div className="gem-chip">
          <div>
            <b>{gem.name}</b>
            {typeof nReactions === 'number' && <span> · {nReactions} reactions</span>}
            <button type="button" className="link" onClick={removeGem}>
              Remove
            </button>
          </div>
          {gem.warnings.map((w) => (
            <p key={w} className="note warn">
              {w === 'no_carbon_source'
                ? 'No carbon source detected. Set a substrate.'
                : w === 'no_objective'
                  ? 'No objective detected. Set one below.'
                  : w}
            </p>
          ))}
        </div>
      )}
      {uploadError && (
        <p className="note err">
          {uploadError.message}
          {uploadError.hint ? ` ${uploadError.hint}` : ''}
        </p>
      )}
      {formError && (
        <p className="note err" role="alert">
          {formError}
        </p>
      )}

      <div className="actions">
        <input
          ref={fileInput}
          type="file"
          hidden
          accept=".json,.xml,.sbml,application/json,text/xml"
          onChange={(e) => void onFile(e.target.files?.[0])}
        />
        <button type="button" className="pill" disabled={uploading} onClick={() => fileInput.current?.click()} title="Upload a COBRA .json or SBML .xml model">
          {uploading ? '…' : 'Upload'}
        </button>
        <button type="submit" className="pill primary" disabled={defaultsLoading}>
          Generate
        </button>
      </div>
    </form>
  )
}
