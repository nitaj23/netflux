import { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState } from 'react'
import type { PathwayResponse } from '../api/types'
import { RW, computeLayout, type PlacedNode, type Rxn } from '../lib/pathwayLayout'
import { geneUrl } from '../lib/links'

// Almost-transparent subsystem films, cycled by index. Easy on the eyes against the dark teal canvas.
const FILM_FILL = ['rgba(247,111,142,0.06)', 'rgba(244,241,187,0.055)', 'rgba(143,212,166,0.06)', 'rgba(230,235,224,0.045)']
const FILM_STROKE = ['rgba(247,111,142,0.24)', 'rgba(244,241,187,0.22)', 'rgba(143,212,166,0.24)', 'rgba(230,235,224,0.20)']

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
  biggId?: string
  organism?: string
}

export function SubgraphView({ data, target, carbonSource, biggId, organism }: Props) {
  const [cofactors, setCofactors] = useState(false)
  const [metLabels, setMetLabels] = useState(true)
  const [hover, setHover] = useState<string | null>(null)
  const [view, setView] = useState<View>({ k: 1, tx: 0, ty: 0 })
  const box = useRef<HTMLDivElement>(null)
  const drag = useRef<{ x: number; y: number; tx: number; ty: number; moved: boolean } | null>(null)
  // Hover open/close with a short grace period so the cursor can travel from node to tooltip
  // (the tooltip is interactive, so its gene links are clickable).
  const hoverTimer = useRef<number>(0)
  const openHover = useCallback((key: string) => {
    clearTimeout(hoverTimer.current)
    setHover(key)
  }, [])
  const closeHover = useCallback(() => {
    clearTimeout(hoverTimer.current)
    hoverTimer.current = window.setTimeout(() => setHover(null), 140)
  }, [])

  const layout = useMemo(
    () => computeLayout(data, target, carbonSource, { cofactors }),
    [data, target, carbonSource, cofactors],
  )

  const fit = useCallback(() => {
    const el = box.current
    if (!el) return
    const { x0, x1, y0, y1 } = layout.bbox
    const w = el.clientWidth
    const h = el.clientHeight
    const dx = x1 - x0
    const dy = y1 - y0
    // Guard against an unlaid-out element or an empty/degenerate bbox (Infinity from Math.min of
    // no nodes) — either would make tx/ty NaN and the whole graph <g> vanish.
    if (w <= 1 || h <= 1 || !(dx > 0) || !(dy > 0) || !Number.isFinite(dx) || !Number.isFinite(dy)) return
    const k = Math.min(2.4, Math.max(0.05, Math.min(w / dx, h / dy)))
    setView({ k, tx: (w - dx * k) / 2 - x0 * k, ty: (h - dy * k) / 2 - y0 * k })
  }, [layout])

  // Refit when the layout changes (filters) and whenever the graph box resizes — opening the
  // details panel, collapsing Parameters, or resizing the window. ResizeObserver fires after
  // layout settles, so fit() always reads the real post-change size (no blank-view race).
  useEffect(() => {
    fit()
  }, [fit])

  useEffect(() => {
    const el = box.current
    if (!el) return
    let raf = 0
    const ro = new ResizeObserver(() => {
      cancelAnimationFrame(raf)
      raf = requestAnimationFrame(() => fit())
    })
    ro.observe(el)
    return () => {
      cancelAnimationFrame(raf)
      ro.disconnect()
    }
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

  const hoverNode = hover ? layout.nodes.get(hover) : undefined
  const colorOf = (r: Rxn) => TAGS[r.design]?.color ?? PLAIN

  return (
    <div className="graph-wrap">
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
            {layout.films.map((f, i) => (
              <g key={f.name} pointerEvents="none">
                <rect x={f.x0} y={f.y0} width={f.x1 - f.x0} height={f.y1 - f.y0} rx={16} fill={FILM_FILL[i % FILM_FILL.length]} stroke={FILM_STROKE[i % FILM_STROKE.length]} strokeWidth={1} />
                <text className="film-label" x={f.x0 + 14} y={f.y0 + 18}>
                  <title>{f.name}</title>
                  {clip(f.name, 30)}
                </text>
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
                <g key={n.key} className="mnode" onMouseEnter={() => openHover(n.key)} onMouseLeave={closeHover}>
                  <circle cx={n.x} cy={n.y} r={n.met?.currency ? 6 : 10} fill={n.met?.currency ? '#294746' : '#bfd0c9'} stroke={hover === n.key ? '#fff' : '#3a5c5a'} strokeWidth={hover === n.key ? 2.5 : 1.3} />
                  {metLabels && (
                    <text className="mlabel" x={n.x} y={n.y + 30}>
                      {clip(n.id, 18)}
                    </text>
                  )}
                </g>
              ))}

            {[...layout.nodes.values()]
              .filter((n) => n.kind === 'r')
              .map((n) => (
                <ReactionGlyph
                  key={n.key}
                  n={n}
                  color={colorOf(n.rxn!)}
                  selected={hover === n.key}
                  onEnter={() => openHover(n.key)}
                  onLeave={closeHover}
                />
              ))}
          </g>
        </svg>
      </div>

      {layout.nodes.size === 0 && (
        <div className="state muted" style={{ position: 'absolute', inset: 0, pointerEvents: 'none' }}>
          No reactions to display for this subgraph.
        </div>
      )}

      <div className="graph-tools">
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
          <span className="key-hi">
            <i className="ring" style={{ borderColor: SRC }} /> <b>Carbon source / substrate</b>
          </span>
          <span className="key-hi">
            <i className="ring" style={{ borderColor: TGT }} /> <b>Target</b>
          </span>
          {layout.designs.map((d) => (
            <span key={d} title={TAGS[d]?.note}>
              <i style={{ background: TAGS[d]?.color ?? PLAIN }} /> {TAGS[d]?.label ?? d}
            </span>
          ))}
          {layout.compartments.length > 0 && <span>Metabolites: {layout.compartments.map((c) => `${c} ${COMPARTMENT[c] ?? ''}`.trim()).join(', ')}</span>}
        </span>
      </div>

      {hoverNode && (
        <Tooltip
          node={hoverNode}
          view={view}
          box={box.current}
          geneHref={(g) => geneUrl(g, biggId, organism)}
          onEnter={() => openHover(hoverNode.key)}
          onLeave={closeHover}
        />
      )}
    </div>
  )
}

function ReactionGlyph({ n, color, selected, onEnter, onLeave }: { n: PlacedNode; color: string; selected: boolean; onEnter: () => void; onLeave: () => void }) {
  const r = n.rxn!
  return (
    <g className={n.cofactorOnly ? 'rgroup cof' : 'rgroup'} onMouseEnter={onEnter} onMouseLeave={onLeave}>
      <title>{`${r.name} (${r.id})`}</title>
      <rect x={n.x - RW - 10} y={n.y - 30} width={2 * RW + 20} height={46} fill="transparent" />
      {n.onSpine && <line x1={n.x - RW} x2={n.x + RW} y1={n.y} y2={n.y} stroke="#c9d6de" strokeWidth={15} strokeLinecap="round" opacity={0.16} />}
      {(r.isSource || r.isTarget) && <rect x={n.x - RW - 7} y={n.y - 13} width={2 * RW + 14} height={26} rx={7} fill="none" stroke={r.isTarget ? TGT : SRC} strokeWidth={2.2} />}
      {selected && <rect x={n.x - RW - 10} y={n.y - 16} width={2 * RW + 20} height={32} rx={8} fill="none" stroke="#fff" strokeWidth={1.8} />}
      <line x1={n.x - RW} x2={n.x + RW} y1={n.y} y2={n.y} stroke={color} strokeWidth={6.5} strokeLinecap="round" />
      <path d={`M${n.x + RW - 1},${n.y - 7} L${n.x + RW + 11},${n.y} L${n.x + RW - 1},${n.y + 7} z`} fill={color} />
      <text className="rlabel" x={n.x} y={n.y - 17}>
        {clip(r.name, 28)}
      </text>
    </g>
  )
}

/** Floating annotation near the hovered node, clamped so the whole box stays inside the canvas. */
function Tooltip({
  node,
  view,
  box,
  geneHref,
  onEnter,
  onLeave,
}: {
  node: PlacedNode
  view: View
  box: HTMLDivElement | null
  geneHref: (gene: string) => string
  onEnter: () => void
  onLeave: () => void
}) {
  const ref = useRef<HTMLDivElement>(null)
  const [pos, setPos] = useState<{ left: number; top: number } | null>(null)

  // Measure the real box size, then place it beside the node but fully within the canvas.
  useLayoutEffect(() => {
    const el = ref.current
    if (!el || !box) return
    const M = 8
    const bw = box.clientWidth
    const bh = box.clientHeight
    const tw = el.offsetWidth
    const th = el.offsetHeight
    const sx = node.x * view.k + view.tx
    const sy = node.y * view.k + view.ty
    let left = sx + 16
    if (left + tw > bw - M) left = sx - 16 - tw // flip to the node's left
    left = Math.min(Math.max(M, left), Math.max(M, bw - tw - M))
    const top = Math.min(Math.max(M, sy + 14), Math.max(M, bh - th - M))
    setPos({ left, top })
  }, [node, view, box])

  return (
    <div
      ref={ref}
      className="node-tip"
      style={pos ? { left: pos.left, top: pos.top } : { left: 0, top: 0, visibility: 'hidden' }}
      aria-live="polite"
      onMouseEnter={onEnter}
      onMouseLeave={onLeave}
    >
      {node.rxn ? <ReactionDetails r={node.rxn} geneHref={geneHref} /> : <MetaboliteDetails id={node.id} comp={node.met?.compartment ?? ''} currency={!!node.met?.currency} />}
    </div>
  )
}

function ReactionDetails({ r, geneHref }: { r: Rxn; geneHref: (gene: string) => string }) {
  const tag = TAGS[r.design]
  const flags = [r.isSource && 'Substrate uptake', r.isTarget && 'Target product'].filter(Boolean).join(' · ')
  return (
    <>
      <h3>{r.name}</h3>
      <p className="muted id mono">{r.id}</p>
      {flags && <p className="flag">{flags}.</p>}
      {tag && (
        <p className="tag" style={{ color: tag.color }}>
          {tag.label} <span className="muted">— {tag.note}.</span>
        </p>
      )}
      <dl>
        <dt>Flux</dt>
        <dd>{fmt(r.flux)}</dd>
        {(r.vMin !== null || r.vMax !== null) && (
          <>
            <dt>FVA range</dt>
            <dd>
              {fmt(r.vMin)} to {fmt(r.vMax)}
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
          {fmt(r.lb)} to {fmt(r.ub)}
        </dd>
        <dt>Subsystem</dt>
        <dd>{r.subsystem || 'Unknown'}</dd>
        <dt>Genes</dt>
        <dd>
          {r.genes.length
            ? r.genes.map((g, i) => (
                <span key={g}>
                  {i > 0 && ', '}
                  <a className="gene-link" href={geneHref(g)} target="_blank" rel="noreferrer">
                    {g}
                  </a>
                </span>
              ))
            : '—'}
        </dd>
        {r.geneRule && (
          <>
            <dt>Gene rule</dt>
            <dd className="mono">{r.geneRule}</dd>
          </>
        )}
        {r.ins.length > 0 && (
          <>
            <dt>Consumes</dt>
            <dd className="mono">{r.ins.join(', ')}</dd>
          </>
        )}
        {r.outs.length > 0 && (
          <>
            <dt>Produces</dt>
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
      <h3 className="mono">{id}</h3>
      <dl>
        <dt>Compartment</dt>
        <dd>{comp ? `${COMPARTMENT[comp] ?? 'Unknown'} (${comp})` : '—'}</dd>
        <dt>Role</dt>
        <dd>{currency ? 'Cofactor / currency metabolite' : 'Carbon carrier'}</dd>
      </dl>
    </>
  )
}

