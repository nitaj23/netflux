import { useCallback, useEffect, useRef, useState } from 'react'
import { ApiError } from './api/client'
import { extractPathway, submitDesign } from './api/endpoints'
import { Header } from './components/Header'
import { ParamPanel } from './components/ParamPanel'
import { ResultsPane } from './components/ResultsPane'
import type { ActiveKey, RunSnapshot, SubgraphTab } from './state'

const asApiError = (e: unknown) => (e instanceof ApiError ? e : new ApiError(0, { detail: String(e) }))

export default function App() {
  const [collapsed, setCollapsed] = useState(false)
  const [tabs, setTabs] = useState<SubgraphTab[]>([])
  const [active, setActive] = useState<ActiveKey>(0)
  const [lastSubgraph, setLastSubgraph] = useState<number | null>(null)
  const counter = useRef(0)
  const started = useRef(new Set<number>())
  const tabsRef = useRef(tabs)
  useEffect(() => {
    tabsRef.current = tabs
  }, [tabs])

  const patch = useCallback((id: number, p: Partial<SubgraphTab>) => setTabs((ts) => ts.map((t) => (t.id === id ? { ...t, ...p } : t))), [])

  const generate = useCallback(
    (snapshot: RunSnapshot) => {
      const id = ++counter.current
      setTabs((ts) => [...ts, { id, title: `Subgraph ${id}`, snapshot, status: 'loading' }])
      setActive(id)
      setLastSubgraph(id)
      extractPathway(snapshot.source, snapshot.target).then(
        (data) => patch(id, { status: 'done', data }),
        (e) => patch(id, { status: 'error', error: asApiError(e) }),
      )
    },
    [patch],
  )

  const startDesign = useCallback(
    (id: number) => {
      const tab = tabsRef.current.find((t) => t.id === id)
      if (!tab || started.current.has(id)) return
      started.current.add(id)
      const s = tab.snapshot
      submitDesign({ source: s.source, target: s.target, biomass: s.biomass, ...(s.substrate ? { substrate: s.substrate } : {}) }).then(
        (r) => patch(id, { jobId: r.job_id, jobStartedAt: Date.now() }),
        (e) => patch(id, { jobSubmitError: asApiError(e) }),
      )
    },
    [patch],
  )

  const retryDesign = useCallback(
    (id: number) => {
      started.current.delete(id)
      patch(id, { jobId: undefined, jobStartedAt: undefined, jobSubmitError: undefined })
    },
    [patch],
  )

  const activate = (k: ActiveKey) => {
    setActive(k)
    if (typeof k === 'number') setLastSubgraph(k)
  }

  const close = (id: number) => {
    const rest = tabs.filter((t) => t.id !== id)
    setTabs(rest)
    const fallback = rest[rest.length - 1]?.id
    if (active === id) setActive(fallback ?? 0)
    if (lastSubgraph === id) setLastSubgraph(fallback ?? null)
  }

  const designTab = tabs.find((t) => t.id === lastSubgraph)

  return (
    <div className="app">
      <Header />
      <main className={collapsed ? 'layout collapsed' : 'layout'}>
        <div className="side">
          <div className="side-head">
            <h2>Parameters</h2>
            <button className="collapse" onClick={() => setCollapsed((c) => !c)} aria-label={collapsed ? 'Show parameters' : 'Hide parameters'}>
              {collapsed ? '>>>' : '<<<'}
            </button>
          </div>
          {/* Kept mounted while collapsed so the form keeps its state. */}
          <div className="panel" hidden={collapsed}>
            <ParamPanel onGenerate={generate} />
          </div>
        </div>
        <ResultsPane
          tabs={tabs}
          active={active}
          designTab={designTab}
          onActivate={activate}
          onClose={close}
          onStartDesign={startDesign}
          onRetryDesign={retryDesign}
        />
      </main>
    </div>
  )
}
