import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import type { PathwayResponse } from '../api/types'
import { RW, computeLayout, type Arrange, type PlacedNode, type Rxn } from '../lib/pathwayLayout'

const TAGS: Record<string, { label: string; color: string; note: string }> = {
  protected: { label: 'Protected', color: '#f2798c', note: 'backbone, coupled to target' },
  nontarget_essential: { label: 'Nontarget essential', color: '#e0a44f', note: 'backbone, not coupled' },
  alt_route: { label: 'Alt route', color: '#6bb2f0', note: 'flexible, coupled' },
  ko_candidate: { label: 'KO candidate', color: '#43c493', note: 'flexible, not coupled' },
}
const PLAIN = '#8fb8c9'
const SRC = '#6aa4ff'
const TGT = '#f2798c'
const LINK = '#5b6b78'
const COMPARTMENT: Record<string, string> = { c: 'cytosol', e: 'extracellular', p: 'periplasm', m: 'mitochondrion', x: 'peroxisome', n: 'nucleus', v: 'vacuole', g: 'golgi', r: 'ER', h: 'chloroplast', l: 'lysosome' }

const fmt = (v: number | null | undefined) => {
  if (v === null || v === undefined) return '—'
  const a = Math.abs(v)
  if (a === 0) return '0'
  return a < 0.001 || a >= 10000 ? v.toExponential(2) : String(Number(v.toPrecision(4)))
}
const clip = (s: string, n: number) => (s.length > n ? `${s.slice(0, n - 1)}…` : s)

interface View {
  k: number
  tx: number
  ty: number
}

interface Props {
  data: PathwayResponse
  target: string
  carbonSource?: string
}

export function SubgraphView({ data, target, carbonSource }: Props) {
  const [arrange, setArrange] = useState<Arrange>('subsystem')
  const [corridor, setCorridor] = useState(true)
  const [cofactors, setCofactors] = useState(false)
  const [metLabels, setMetLabels] = useState(true)
  const [picked, setPicked] = useState<string | null>(null)
  const [view, setView] = useState<View>({ k: 1, tx: 0, ty: 0 })
  const box = useRef<HTMLDivElement>(null)
  const drag = useRef<{ x: number; y: number; tx: number; ty: number; moved: boolean } | null>(null)

  const layout = useMemo(
    () => computeLayout(data, target, carbonSource, { arrange, corridor, cofactors }),
    [data, target, carbonSource, arrange, corridor, cofactors],
  )

  const fit = useCallback(() => {
    const el = box.current
    if (!el) return
    const { x0, x1, y0, y1 } = layout.bbox
    const w = el.clientWidth
    const h = el.clientHeight
    const k = Math.min(1.3, Math.max(0.05, Math.min(w / (x1 - x0), h / (y1 - y0))))
    setView({ k, tx: (w - (x1 - x0) * k) / 2 - x0 * k, ty: (h - (y1 - y0) * k) / 2 - y0 * k })
  }, [layout])

  useEffect(() => {
    fit()
  }, [fit])

  // Wheel needs a non-passive listener to stop the page scrolling.
  useEffect(() => {
    const el = box.current
    if (!el) return
    const onWheel = (e: WheelEvent) => {
      e.preventDefault()
      const r = el.getBoundingClientRect()
      const px = e.clientX - r.left
      const py = e.clientY - r.top
      const f = Math.exp(-e.deltaY * 0.0015)
      setView((v) => {
        const k = Math.min(4, Math.max(0.05, v.k * f))
        return { k, tx: px - ((px - v.tx) / v.k) * k, ty: py - ((py - v.ty) / v.k) * k }
      })
    }
    el.addEventListener('wheel', onWheel, { passive: false })
    return () => el.removeEventListener('wheel', onWheel)
  }, [])

  const zoom = (f: number) => {
    const el = box.current
    if (!el) return
    const px = el.clientWidth / 2
    const py = el.clientHeight / 2
    setView((v) => {
      const k = Math.min(4, Math.max(0.05, v.k * f))
      return { k, tx: px - ((px - v.tx) / v.k) * k, ty: py - ((py - v.ty) / v.k) * k }
    })
  }

  const selected = picked ? layout.nodes.get(picked) : undefined
  const colorOf = (r: Rxn) => TAGS[r.design]?.color ?? PLAIN

  return (
    <div className={selected ? 'graph-wrap side' : 'graph-wrap'}>
      <div
        ref={box}
        className="graph"
        onPointerDown={(e) => {
          drag.current = { x: e.clientX, y: e.clientY, tx: view.tx, ty: view.ty, moved: false }
          e.currentTarget.setPointerCapture(e.pointerId)
        }}
        onPointerMove={(e) => {
          const d = drag.current
          if (!d) return
          const dx = e.clientX - d.x
          const dy = e.clientY - d.y
          if (Math.abs(dx) + Math.abs(dy) > 4) d.moved = true
          if (d.moved) setView((v) => ({ ...v, tx: d.tx + dx, ty: d.ty + dy }))
        }}
        onPointerUp={() => {
          // A press that never moved is a click on the background.
          if (drag.current && !drag.current.moved) setPicked(null)
          drag.current = null
        }}
      >
        <svg width="100%" height="100%" role="img" aria-label="Pathway map">
          <defs>
            <marker id="ah" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="5" markerHeight="5" orient="auto">
              <path d="M0,0 L8,4 L0,8 z" fill={LINK} />
            </marker>
            <marker id="ahs" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="5" markerHeight="5" orient="auto">
              <path d="M0,0 L8,4 L0,8 z" fill="#c9d6de" />
            </marker>
          </defs>
          <g transform={`translate(${view.tx} ${view.ty}) scale(${view.k})`}>
            {arrange === 'subsystem' &&
              layout.lanes.map((l, i) => (
                <g key={l.name}>
                  <rect x={layout.bbox.x0} width={layout.bbox.x1 - layout.bbox.x0} y={l.top} height={l.bottom - l.top} fill={i % 2 ? 'rgba(255,255,255,0.045)' : 'rgba(255,255,255,0.015)'} />
                  <LaneLabel x={layout.bbox.x0 + 10} y={(l.top + l.bottom) / 2} name={l.name} />
                </g>
              ))}

            {layout.edges.map((e) => {
              const a = layout.nodes.get(e.s)!
              const b = layout.nodes.get(e.t)!
              const x1 = a.x + (a.kind === 'r' ? RW + 5 : 6)
              const x2 = b.x - (b.kind === 'r' ? RW : 7)
              const mx = (x1 + x2) / 2
              const w = 0.7 + 2.8 * Math.sqrt(e.mass / layout.maxMass)
              return (
                <path
                  key={`${e.s}>${e.t}`}
                  d={`M${x1},${a.y} C${mx},${a.y} ${mx},${b.y} ${x2},${b.y}`}
                  fill="none"
                  stroke={e.spine ? '#c9d6de' : LINK}
                  strokeWidth={e.spine ? Math.max(2, w) : w}
                  strokeOpacity={e.spine ? 0.9 : 0.7}
                  markerEnd={e.spine ? 'url(#ahs)' : 'url(#ah)'}
                />
              )
            })}

            {[...layout.nodes.values()]
              .filter((n) => n.kind === 'm')
              .map((n) => (
                <g key={n.key} className="mnode" onClick={(ev) => (ev.stopPropagation(), setPicked(n.key))} onPointerDown={(ev) => ev.stopPropagation()}>
                  <title>{`${n.id}${n.met?.compartment ? ` · ${COMPARTMENT[n.met.compartment] ?? n.met.compartment}` : ''}${n.met?.currency ? ' · cofactor' : ''}`}</title>
                  <circle cx={n.x} cy={n.y} r={n.met?.currency ? 3.5 : 5.5} fill={n.met?.currency ? '#2a3640' : '#a9bac4'} stroke={picked === n.key ? '#fff' : '#33454f'} strokeWidth={picked === n.key ? 2 : 1} />
                  {metLabels && (
                    <text className="mlabel" x={n.x} y={n.y + 16}>
                      {clip(n.id, 14)}
                    </text>
                  )}
                </g>
              ))}

            {[...layout.nodes.values()]
              .filter((n) => n.kind === 'r')
              .map((n) => (
                <ReactionGlyph key={n.key} n={n} color={colorOf(n.rxn!)} selected={picked === n.key} onPick={() => setPicked(n.key)} />
              ))}
          </g>
        </svg>
      </div>

      <div className="graph-tools">
        <div className="seg" role="group" aria-label="Arrangement">
          <button aria-pressed={arrange === 'subsystem'} onClick={() => setArrange('subsystem')}>
            By subsystem
          </button>
          <button aria-pressed={arrange === 'flow'} onClick={() => setArrange('flow')}>
            Flow only
          </button>
        </div>
        <button className="chip-btn" aria-pressed={corridor} onClick={() => setCorridor((c) => !c)} title="Show the dominant substrate-to-target route plus reactions carrying at least 10% of the target flux">
          Main route only
        </button>
        <button className="chip-btn" aria-pressed={cofactors} onClick={() => setCofactors((c) => !c)} title="Draw links that exist only through ATP, NADH, H⁺ and similar">
          Cofactors
        </button>
        <button className="chip-btn" aria-pressed={metLabels} onClick={() => setMetLabels((c) => !c)}>
          Metabolite labels
        </button>
      </div>

      <div className="zoom">
        <button onClick={() => zoom(1.25)} aria-label="Zoom in">
          +
        </button>
        <button onClick={() => zoom(0.8)} aria-label="Zoom out">
          –
        </button>
        <button onClick={fit} aria-label="Fit to view" title="Fit to view" className="fit">
          ⤢
        </button>
      </div>

      <div className="graph-legend">
        <span>
          Showing {layout.shown} of {layout.total} reactions · flow runs left to right · line width = mass flow
        </span>
        <span className="keys">
          <i className="ring" style={{ borderColor: SRC }} /> substrate <i className="ring" style={{ borderColor: TGT }} /> target
          {layout.designs.map((d) => (
            <span key={d} title={TAGS[d]?.note}>
              <i style={{ background: TAGS[d]?.color ?? PLAIN }} /> {TAGS[d]?.label ?? d}
            </span>
          ))}
          {layout.compartments.length > 0 && <span>Metabolites: {layout.compartments.map((c) => `${c} ${COMPARTMENT[c] ?? ''}`.trim()).join(', ')}</span>}
        </span>
      </div>

      {selected && <Details node={selected} onClose={() => setPicked(null)} />}
    </div>
  )
}

/** Subsystem name in the left gutter, wrapped onto at most two lines. */
function LaneLabel({ x, y, name }: { x: number; y: number; name: string }) {
  const words = name.split(/\s+/)
  const lines: string[] = ['']
  for (const w of words) {
    const cur = lines[lines.length - 1]
    if (cur && (cur + ' ' + w).length > 22 && lines.length < 2) lines.push(w)
    else lines[lines.length - 1] = cur ? `${cur} ${w}` : w
  }
  return (
    <text className="lane-label" x={x} y={y - (lines.length - 1) * 7 + 4}>
      <title>{name}</title>
      {lines.map((ln, i) => (
        <tspan key={i} x={x} dy={i ? 14 : 0}>
          {clip(ln, 24)}
        </tspan>
      ))}
    </text>
  )
}

function ReactionGlyph({ n, color, selected, onPick }: { n: PlacedNode; color: string; selected: boolean; onPick: () => void }) {
  const r = n.rxn!
  return (
    <g className={n.cofactorOnly ? 'rgroup cof' : 'rgroup'} onClick={(e) => (e.stopPropagation(), onPick())} onPointerDown={(e) => e.stopPropagation()}>
      <title>{`${r.name} (${r.id})`}</title>
      <rect x={n.x - RW - 8} y={n.y - 20} width={2 * RW + 16} height={30} fill="transparent" />
      {n.onSpine && <line x1={n.x - RW} x2={n.x + RW} y1={n.y} y2={n.y} stroke="#c9d6de" strokeWidth={9} strokeLinecap="round" opacity={0.16} />}
      {(r.isSource || r.isTarget) && <rect x={n.x - RW - 5} y={n.y - 9} width={2 * RW + 10} height={18} rx={4} fill="none" stroke={r.isTarget ? TGT : SRC} strokeWidth={1.4} />}
      {selected && <rect x={n.x - RW - 7} y={n.y - 11} width={2 * RW + 14} height={22} rx={5} fill="none" stroke="#fff" strokeWidth={1.6} />}
      <line x1={n.x - RW} x2={n.x + RW} y1={n.y} y2={n.y} stroke={color} strokeWidth={3.6} strokeLinecap="round" />
      <path d={`M${n.x + RW - 1},${n.y - 4} L${n.x + RW + 6},${n.y} L${n.x + RW - 1},${n.y + 4} z`} fill={color} />
      <text className="rlabel" x={n.x} y={n.y - 9}>
        {clip(r.name, 24)}
      </text>
    </g>
  )
}

function Details({ node, onClose }: { node: PlacedNode; onClose: () => void }) {
  return (
    <aside className="node-info" aria-live="polite">
      <button className="icon-btn" onClick={onClose} aria-label="Close details">
        ×
      </button>
      {node.rxn ? <ReactionDetails r={node.rxn} /> : <MetaboliteDetails id={node.id} comp={node.met?.compartment ?? ''} currency={!!node.met?.currency} />}
    </aside>
  )
}

function ReactionDetails({ r }: { r: Rxn }) {
  const tag = TAGS[r.design]
  const flags = [r.isSource && 'substrate uptake', r.isTarget && 'target'].filter(Boolean).join(' · ')
  return (
    <>
      <h3>{r.name}</h3>
      <p className="muted id">{r.id}</p>
      {flags && <p className="flag">{flags}</p>}
      {tag && (
        <p className="tag" style={{ color: tag.color }}>
          {tag.label} <span className="muted">· {tag.note}</span>
        </p>
      )}
      <dl>
        <dt>Flux</dt>
        <dd>{fmt(r.flux)}</dd>
        {(r.vMin !== null || r.vMax !== null) && (
          <>
            <dt>FVA range</dt>
            <dd>
              {fmt(r.vMin)} .. {fmt(r.vMax)}
              {r.cls && <span className="muted"> ({r.cls})</span>}
            </dd>
          </>
        )}
        {r.coupling !== null && (
          <>
            <dt>Coupling</dt>
            <dd>{fmt(r.coupling)}</dd>
          </>
        )}
        <dt>Bounds</dt>
        <dd>
          {fmt(r.lb)} .. {fmt(r.ub)}
        </dd>
        <dt>Subsystem</dt>
        <dd>{r.subsystem}</dd>
        <dt>Genes</dt>
        <dd>{r.genes.length ? r.genes.join(', ') : '—'}</dd>
        {r.geneRule && (
          <>
            <dt>Gene rule</dt>
            <dd className="mono">{r.geneRule}</dd>
          </>
        )}
        {r.ins.length > 0 && (
          <>
            <dt>Takes in</dt>
            <dd className="mono">{r.ins.join(', ')}</dd>
          </>
        )}
        {r.outs.length > 0 && (
          <>
            <dt>Passes on</dt>
            <dd className="mono">{r.outs.join(', ')}</dd>
          </>
        )}
      </dl>
    </>
  )
}

function MetaboliteDetails({ id, comp, currency }: { id: string; comp: string; currency: boolean }) {
  return (
    <>
      <h3>{id}</h3>
      <dl>
        <dt>Compartment</dt>
        <dd>{comp ? `${comp} (${COMPARTMENT[comp] ?? 'unknown'})` : '—'}</dd>
        <dt>Role</dt>
        <dd>{currency ? 'cofactor' : 'carbon carrier'}</dd>
      </dl>
    </>
  )
}

