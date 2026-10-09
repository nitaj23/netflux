import { useQuery } from '@tanstack/react-query'
import { useState } from 'react'
import { ApiError } from '../api/client'
import { fetchJob, modelDownloadUrl } from '../api/endpoints'
import { base64ToBytes, downloadUrl, safeName, saveBlob } from '../lib/download'
import type { ActiveKey, SubgraphTab } from '../state'
import { InterventionsView } from './InterventionsView'
import { SubgraphView } from './SubgraphView'

interface Props {
  tabs: SubgraphTab[]
  active: ActiveKey
  designTab: SubgraphTab | undefined
  onActivate: (k: ActiveKey) => void
  onClose: (id: number) => void
  onStartDesign: (id: number) => void
  onRetryDesign: (id: number) => void
}

export function ResultsPane({ tabs, active, designTab, onActivate, onClose, onStartDesign, onRetryDesign }: Props) {
  const [menu, setMenu] = useState(false)
  const current = typeof active === 'number' ? tabs.find((t) => t.id === active) : undefined

  // Same key as InterventionsView; enabled:false so this only reads the shared cache.
  const job = useQuery({
    queryKey: ['job', designTab?.jobId],
    queryFn: () => fetchJob(designTab!.jobId!),
    enabled: false,
  })
  const designResult = job.data?.status === 'completed' ? job.data.result : null

  type Item = { label: string; run: () => void }
  const items: Item[] = []
  if (active === 'interventions') {
    if (designTab && designResult) {
      const base = safeName(`netflux_${designTab.snapshot.target}`)
      items.push({ label: 'Designs (.csv)', run: () => saveBlob(`${base}_designs.csv`, designResult.csv, 'text/csv') })
      if (designResult.hypotheses_csv)
        items.push({
          label: 'Hypotheses (.csv)',
          run: () => saveBlob(`${base}_hypotheses.csv`, designResult.hypotheses_csv!, 'text/csv'),
        })
    }
  } else if (current?.data) {
    const d = current.data
    const base = safeName(`netflux_${current.snapshot.target}_subgraph${current.id}`)
    items.push({
      label: 'Subgraph (.json)',
      run: () => saveBlob(`${base}.json`, JSON.stringify(d.subgraph, null, 2), 'application/json'),
    })
    if (d.graphml_base64)
      items.push({
        label: 'GraphML (.graphml)',
        run: () => saveBlob(`${base}.graphml`, base64ToBytes(d.graphml_base64!), 'application/graphml+xml'),
      })
    const modelUrl = modelDownloadUrl(current.snapshot.source)
    if (modelUrl) items.push({ label: 'Base model (.json)', run: () => downloadUrl(modelUrl) })
  }

  const showEmptySubgraph = tabs.length === 0 && active !== 'interventions'

  return (
    <section className="results">
      <div className="tabbar" role="tablist">
        {showEmptySubgraph && (
          <button role="tab" aria-selected className="tab on">
            Subgraph
          </button>
        )}
        {tabs.map((t) => (
          <div key={t.id} className={active === t.id ? 'tab on' : 'tab'}>
            <button role="tab" aria-selected={active === t.id} onClick={() => onActivate(t.id)}>
              {t.title}
              {t.status === 'loading' && <span className="spinner sm" />}
            </button>
            <button className="tab-x" aria-label={`Close ${t.title}`} onClick={() => onClose(t.id)}>
              ×
            </button>
          </div>
        ))}
        <button
          role="tab"
          aria-selected={active === 'interventions'}
          className={active === 'interventions' ? 'tab on' : 'tab'}
          onClick={() => onActivate('interventions')}
        >
          Interventions
        </button>
      </div>

      <div className="canvas">
        {active === 'interventions' ? (
          <InterventionsView tab={designTab} onStart={onStartDesign} onRetry={onRetryDesign} />
        ) : !current ? (
          <div className="state muted">Choose a compound and organism, then press Generate.</div>
        ) : current.status === 'loading' ? (
          <div className="state">
            <span className="spinner" /> Extracting pathway…
          </div>
        ) : current.status === 'error' ? (
          <div className="state err" role="alert">
            <h3>Could not extract a subgraph</h3>
            <p>{current.error?.message}</p>
            {current.error instanceof ApiError && current.error.hint && <p className="hint">{current.error.hint}</p>}
          </div>
        ) : (
          <SubgraphView data={current.data!} target={current.snapshot.target} carbonSource={current.snapshot.carbonSource} biggId={current.snapshot.biggId} organism={current.snapshot.organism} />
        )}
      </div>

      <div className="footer">
        <div className="dl">
          <button className="pill dark" disabled={items.length === 0} aria-haspopup="menu" aria-expanded={menu} onClick={() => setMenu((m) => !m)}>
            Download
          </button>
          {menu && items.length > 0 && (
            <ul className="menu" role="menu" onMouseLeave={() => setMenu(false)}>
              {items.map((i) => (
                <li key={i.label} role="none">
                  <button
                    role="menuitem"
                    onClick={() => {
                      i.run()
                      setMenu(false)
                    }}
                  >
                    {i.label}
                  </button>
                </li>
              ))}
            </ul>
          )}
        </div>
      </div>
    </section>
  )
}
