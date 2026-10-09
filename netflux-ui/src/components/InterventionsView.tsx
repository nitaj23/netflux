import { useQuery } from '@tanstack/react-query'
import { useEffect } from 'react'
import { ApiError } from '../api/client'
import { fetchJob } from '../api/endpoints'
import type { ResultTable } from '../api/types'
import { useNow } from '../lib/hooks'
import type { SubgraphTab } from '../state'

const TERMINAL = ['completed', 'failed']

function cell(v: unknown): string {
  if (v === null || v === undefined || v === '') return '—'
  if (Array.isArray(v)) return v.length ? v.join(' + ') : '—'
  if (typeof v === 'number') return Number.isInteger(v) ? String(v) : String(Number(v.toPrecision(4)))
  if (typeof v === 'boolean') return v ? 'yes' : 'no'
  if (typeof v === 'object') return JSON.stringify(v)
  return String(v)
}

/** Render a "g1 + g2 (partial)" gene cell with each gene id linked; suffix text kept as-is. */
function GeneCell({ value, href }: { value: string; href: (gene: string) => string }) {
  if (value === '—') return <>—</>
  return (
    <>
      {value.split(' + ').map((seg, i) => {
        const sp = seg.indexOf(' ') // split a leading gene id from any trailing " (partial)"
        const gene = sp === -1 ? seg : seg.slice(0, sp)
        const rest = sp === -1 ? '' : seg.slice(sp)
        return (
          <span key={i}>
            {i > 0 && ' + '}
            <a className="gene-link" href={href(gene)} target="_blank" rel="noreferrer">
              {gene}
            </a>
            {rest}
          </span>
        )
      })}
    </>
  )
}

function Table({ t, empty, geneHref }: { t: ResultTable; empty: string; geneHref: (gene: string) => string }) {
  if (!t.rows.length) return <p className="muted">{empty}</p>
  return (
    <div className="table-scroll">
      <table>
        <thead>
          <tr>
            {t.columns.map((c) => (
              <th key={c.key}>{c.label}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {t.rows.map((r, i) => (
            <tr key={i}>
              {t.columns.map((c) => (
                <td key={c.key}>{c.key === 'genes' ? <GeneCell value={cell(r[c.key])} href={geneHref} /> : cell(r[c.key])}</td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}

function ErrorBox({ title, e }: { title: string; e: { error?: string; type?: string; detail?: string; message?: string; hint?: string | null } }) {
  return (
    <div className="state err" role="alert">
      <h3>{title}</h3>
      <p>{e.detail ?? e.message}</p>
      {(e.error ?? e.type) && <p className="muted">{e.error ?? e.type}</p>}
      {e.hint && <p className="hint">{e.hint}</p>}
    </div>
  )
}

interface Props {
  tab: SubgraphTab | undefined
  onStart: (tabId: number) => void
  onRetry: (tabId: number) => void
}

export function InterventionsView({ tab, onStart, onRetry }: Props) {
  const tabId = tab?.id
  const needsStart = !!tab && tab.status === 'done' && !tab.jobId && !tab.jobSubmitError

  // Opening the tab is what submits the job (once per Subgraph).
  useEffect(() => {
    if (tabId !== undefined && needsStart) onStart(tabId)
  }, [tabId, needsStart, onStart])

  const job = useQuery({
    queryKey: ['job', tab?.jobId],
    queryFn: () => fetchJob(tab!.jobId!),
    enabled: !!tab?.jobId,
    refetchInterval: (q) => (q.state.data && TERMINAL.includes(q.state.data.status) ? false : 3000),
    staleTime: Infinity,
    gcTime: Infinity,
  })
  const running = !!tab?.jobId && !(job.data && TERMINAL.includes(job.data.status))
  const now = useNow(running)

  if (!tab) return <div className="state">Generate a Subgraph first, then open Interventions to rank strain designs for it.</div>
  if (tab.status === 'loading') return <div className="state">Waiting for {tab.title} to finish…</div>
  if (tab.status === 'error') return <div className="state">{tab.title} failed, so there is nothing to rank.</div>

  // label is "Organism → Compound"; italicise the organism (scientific-name convention).
  const [org, ...rest] = tab.snapshot.label.split(' → ')
  const header = (
    <p className="muted for">
      For <b>{tab.title}</b> · <i>{org}</i>
      {rest.length > 0 && <> → {rest.join(' → ')}</>}
    </p>
  )

  if (tab.jobSubmitError)
    return (
      <div className="pane-scroll">
        {header}
        <ErrorBox title="Could not start the design job" e={tab.jobSubmitError} />
        <button className="pill" onClick={() => onRetry(tab.id)}>
          Retry
        </button>
      </div>
    )

  if (!tab.jobId || (job.isPending && !job.isError))
    return (
      <div className="pane-scroll">
        {header}
        <div className="state">
          <span className="spinner" /> Submitting job…
        </div>
      </div>
    )

  if (job.isError)
    return (
      <div className="pane-scroll">
        {header}
        <ErrorBox title="Lost track of the design job" e={job.error instanceof ApiError ? job.error : { message: String(job.error) }} />
        <button className="pill" onClick={() => onRetry(tab.id)}>
          Run again
        </button>
      </div>
    )

  const j = job.data!
  if (j.status === 'queued' || j.status === 'running') {
    const s = Math.max(0, Math.floor((now - (tab.jobStartedAt ?? now)) / 1000))
    return (
      <div className="pane-scroll">
        {header}
        <div className="state">
          <span className="spinner" />
          <h3>{j.status === 'queued' ? 'Queued' : 'Ranking designs'}</h3>
          <p className="muted">
            {Math.floor(s / 60)}:{String(s % 60).padStart(2, '0')} elapsed. This usually takes a few minutes; you can
            switch tabs and come back.
          </p>
        </div>
      </div>
    )
  }

  if (j.status === 'failed')
    return (
      <div className="pane-scroll">
        {header}
        <ErrorBox title="Design run failed" e={j.error ?? { message: 'The job failed without details.' }} />
        <button className="pill" onClick={() => onRetry(tab.id)}>
          Run again
        </button>
      </div>
    )

  const r = j.result!
  const rejected = r.rejected && r.rejected.n > 0 ? Object.entries(r.rejected.reasons).map(([k, v]) => `${v} ${k.replace(/_/g, ' ')}`).join(', ') : null

  // Repo (BiGG) models: link each gene to its exact BiGG gene page; otherwise an NCBI gene search
  // scoped by organism (uploaded GEMs can use any id scheme).
  const { biggId, organism } = tab.snapshot
  const geneHref = (gene: string) =>
    biggId
      ? `https://bigg.ucsd.edu/models/${encodeURIComponent(biggId)}/genes/${encodeURIComponent(gene)}`
      : `https://www.ncbi.nlm.nih.gov/gene/?term=${encodeURIComponent(`${gene} ${organism ?? ''}`.trim())}`
  return (
    <div className="pane-scroll">
      {header}
      {r.baseline && (
        <div className="chips">
          {Object.entries(r.baseline).map(([k, v]) => (
            <span key={k} className="chip">
              <small>{k.replace(/_/g, ' ')}</small> {cell(v)}
            </span>
          ))}
        </div>
      )}
      <h3>Ranked designs</h3>
      <Table t={r.table} empty={`No verified designs were found${rejected ? ` (rejected: ${rejected})` : ''}. See the leads below.`} geneHref={geneHref} />
      <h3>Hypotheses (leads)</h3>
      <Table t={r.hypotheses_table} empty="No hypotheses were produced." geneHref={geneHref} />
      {r.notes && r.notes.length > 0 && (
        <ul className="muted">
          {r.notes.map((n) => (
            <li key={n}>{n}</li>
          ))}
        </ul>
      )}
    </div>
  )
}
