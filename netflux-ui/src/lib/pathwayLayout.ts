// Escher-style pathway map layout. Reactions become labelled arrows between metabolite
// circles; flow runs left to right. Ported from scripts/pathway_map_template.html, plus
// subsystem lanes. Pure functions: no DOM, no React.
import type { PathwayResponse } from '../api/types'

export type Arrange = 'subsystem' | 'flow'
export interface LayoutOpts {
  arrange: Arrange
  /** Keep only reactions on a source -> target route. */
  corridor: boolean
  /** Draw links that exist only through cofactors (ATP, NADH, H+ ...). */
  cofactors: boolean
}

export interface Rxn {
  id: string
  name: string
  flux: number
  subsystem: string
  geneRule: string
  genes: string[]
  lb: number
  ub: number
  design: string
  cls: string
  coupling: number | null
  vMin: number | null
  vMax: number | null
  isTarget: boolean
  isSource: boolean
  /** Metabolites consumed / produced, from the incident mass-flow edges. */
  ins: string[]
  outs: string[]
}

export interface Met {
  id: string
  compartment: string
  currency: boolean
}

export interface PlacedNode {
  key: string
  kind: 'r' | 'm'
  id: string
  x: number
  y: number
  rxn?: Rxn
  met?: Met
  /** Reaction that shares only cofactors with the rest of the map. */
  cofactorOnly?: boolean
  onSpine?: boolean
}

export interface PlacedEdge {
  s: string
  t: string
  mass: number
  spine: boolean
}

export interface Lane {
  name: string
  top: number
  bottom: number
}

export interface Layout {
  nodes: Map<string, PlacedNode>
  edges: PlacedEdge[]
  lanes: Lane[]
  bbox: { x0: number; y0: number; x1: number; y1: number }
  maxMass: number
  shown: number
  total: number
  compartments: string[]
  designs: string[]
}

export const COLW = 150
export const ROWH = 40
export const LABEL_W = 180
export const RW = 22
const PAD = 6
/** "Main route" keeps reactions carrying at least this share of the target's flux. */
export const MAIN_FLUX_SHARE = 0.1

const CURRENCY = new Set([
  'h', 'h2o', 'atp', 'adp', 'amp', 'pi', 'ppi', 'nad', 'nadh', 'nadp', 'nadph', 'co2', 'o2', 'coa', 'nh4', 'fad',
  'fadh2', 'q8', 'q8h2', 'h2o2', 'gtp', 'gdp', 'udp', 'utp', 'ctp', 'cdp', 'so4', 'na1', 'k', 'mg2', 'fe2', 'fe3',
  'cl', 'thf', 'mlthf', '10fthf', 'so3', 'no3', 'no2', 'n2', 'h2', 'hco3', 'acp',
])

const SUFFIX = /^(.*)_([a-z]{1,3})$/

/** BiGG ids end in a compartment (`succ_c`); `M_`-style ids (`M_SUC`) have none. */
function parseMet(id: string): Met {
  const s = id.replace(/^M_/, '')
  const m = SUFFIX.exec(s)
  const base = (m ? m[1] : s).toLowerCase()
  return { id, compartment: m ? m[2] : '', currency: CURRENCY.has(base) }
}

export function computeLayout(data: PathwayResponse, target: string, carbonSource: string | undefined, opts: LayoutOpts): Layout {
  const per = data.annotation?.per_reaction ?? {}
  const rx = new Map<string, Rxn>()
  for (const n of data.subgraph.nodes) {
    const a = per[n.reaction_id]
    const op = a?.per_op?.[a.per_op.length - 1]
    rx.set(n.reaction_id, {
      id: n.reaction_id,
      name: n.reaction_name || n.reaction_id,
      flux: n.flux ?? 0,
      subsystem: n.subsystem || 'Unknown',
      geneRule: n.gene_rule ?? '',
      genes: n.genes ?? [],
      lb: n.lower_bound ?? 0,
      ub: n.upper_bound ?? 0,
      design: op?.design ?? '',
      cls: op?.cls ?? '',
      coupling: a?.coupling ?? null,
      vMin: op?.v_min ?? null,
      vMax: op?.v_max ?? null,
      isTarget: n.reaction_id === target,
      isSource: false,
      ins: [],
      outs: [],
    })
  }
  const edgesIn = data.subgraph.edges.filter((e) => rx.has(e.source) && rx.has(e.target))
  const metMeta = new Map<string, Met>()
  const meta = (id: string) => {
    let m = metMeta.get(id)
    if (!m) metMeta.set(id, (m = parseMet(id)))
    return m
  }

  const fwd = new Map<string, string[]>()
  const rev = new Map<string, string[]>()
  for (const id of rx.keys()) {
    fwd.set(id, [])
    rev.set(id, [])
  }
  const freq = new Map<string, number>()
  for (const e of edgesIn) {
    fwd.get(e.source)!.push(e.target)
    rev.get(e.target)!.push(e.source)
    for (const m of e.metabolites ?? []) {
      freq.set(m, (freq.get(m) ?? 0) + 1)
      const a = rx.get(e.source)!
      const b = rx.get(e.target)!
      if (!a.outs.includes(m)) a.outs.push(m)
      if (!b.ins.includes(m)) b.ins.push(m)
    }
  }

  // Sources: the chosen carbon source if it is on the map, else every reaction nothing feeds.
  let sources = carbonSource && rx.has(carbonSource) ? [carbonSource] : [...rx.keys()].filter((id) => rev.get(id)!.length === 0 && id !== target)
  if (!sources.length) sources = [[...rx.keys()][0]].filter(Boolean)
  for (const s of sources) rx.get(s)!.isSource = true

  const bfs = (starts: string[], g: Map<string, string[]>) => {
    const seen = new Set(starts.filter((s) => g.has(s)))
    const q = [...seen]
    while (q.length) {
      const n = q.shift()!
      for (const m of g.get(n) ?? [])
        if (!seen.has(m)) {
          seen.add(m)
          q.push(m)
        }
    }
    return seen
  }
  const toTgt = bfs([target], rev)
  const routes = (from: string[]) => {
    const fromSrc = bfs(from, fwd)
    return new Set([...fromSrc].filter((id) => toTgt.has(id)))
  }
  let corridor = routes(sources)
  if (!corridor.has(target)) {
    // The chosen substrate does not reach the target here: fall back to every natural entry point.
    const entries = [...rx.keys()].filter((id) => rev.get(id)!.length === 0 && id !== target)
    sources = entries.length ? entries : sources
    for (const r of rx.values()) r.isSource = sources.includes(r.id)
    corridor = routes(sources)
  }
  if (!corridor.has(target)) corridor = new Set(rx.keys()) // truly disconnected: show everything rather than nothing

  // Dominant-flux source -> target route (weights log(max/mf)); always drawn.
  const spine = dominantRoute(edgesIn, sources, target, [...rx.keys()])
  const spineEdge = new Set<string>()
  for (let i = 0; i < spine.length - 1; i++) spineEdge.add(`${spine[i]}>${spine[i + 1]}`)
  const spineNode = new Set(spine)

  // Ranks use every edge so positions stay put when filters toggle.
  const rank = new Map<string, number>()
  let q = [...sources]
  for (const s of sources) rank.set(s, 0)
  while (q.length) {
    const n = q.shift()!
    for (const m of fwd.get(n)!)
      if (!rank.has(m)) {
        rank.set(m, rank.get(n)! + 1)
        q.push(m)
      }
  }
  const maxF = Math.max(0, ...rank.values())
  const back = new Map<string, number>([[target, 0]])
  q = [target]
  while (q.length) {
    const n = q.shift()!
    for (const m of rev.get(n)!)
      if (!back.has(m)) {
        back.set(m, back.get(n)! + 1)
        q.push(m)
      }
  }
  for (const id of rx.keys()) if (!rank.has(id)) rank.set(id, back.has(id) ? Math.max(0, maxF - back.get(id)!) : maxF)
  rank.set(target, Math.max(...rank.values()) + 1)

  // Main route: the dominant route plus reactions carrying a real share of the product flux,
  // kept only while they stay connected to that route (no orphan islands).
  const tFlux = Math.abs(rx.get(target)!.flux) || Math.max(...[...rx.values()].map((r) => Math.abs(r.flux)))
  const significant = new Set([...corridor].filter((id) => spineNode.has(id) || Math.abs(rx.get(id)!.flux) >= MAIN_FLUX_SHARE * tFlux))
  let main = significant
  if (spine.length) {
    const adj = new Map<string, string[]>([...significant].map((id) => [id, []]))
    for (const e of edgesIn)
      if (significant.has(e.source) && significant.has(e.target)) {
        adj.get(e.source)!.push(e.target)
        adj.get(e.target)!.push(e.source)
      }
    main = bfs(spine, adj)
  }
  const kept = opts.corridor ? main : new Set(rx.keys())

  // One key metabolite per edge: the rarest non-cofactor (hubs like CoA would collapse the map).
  const links: { a: string; b: string; met: string; mass: number; spine: boolean }[] = []
  for (const e of edgesIn) {
    if (!kept.has(e.source) || !kept.has(e.target)) continue
    const onSpine = spineEdge.has(`${e.source}>${e.target}`)
    const all = e.metabolites ?? []
    let c = all.filter((m) => !meta(m).currency)
    if (!c.length) c = opts.cofactors || onSpine ? all.slice() : []
    if (!c.length) continue
    c.sort((x, y) => (freq.get(x)! - freq.get(y)!) || x.localeCompare(y))
    links.push({ a: e.source, b: e.target, met: c[0], mass: e.mass_flow ?? 0, spine: onSpine })
  }
  const linked = new Set<string>()
  for (const l of links) {
    linked.add(l.a)
    linked.add(l.b)
  }

  // Bipartite nodes: reactions on integer ranks, metabolites on half-integer ranks between them.
  type Tmp = { key: string; kind: 'r' | 'm'; id: string; rank: number; lane: string }
  const tmp = new Map<string, Tmp>()
  const laneOf = (r: Rxn) => (opts.arrange === 'subsystem' ? r.subsystem : '')
  for (const id of kept) {
    // A reaction linked to the map only through cofactors would float alone; show it only with "Cofactors" on.
    if (!linked.has(id) && !spineNode.has(id)) continue
    tmp.set(`R:${id}`, { key: `R:${id}`, kind: 'r', id, rank: rank.get(id)!, lane: laneOf(rx.get(id)!) })
  }
  const metRanks = new Map<string, number[]>()
  const metLane = new Map<string, string>()
  for (const l of links) {
    ;(metRanks.get(l.met) ?? metRanks.set(l.met, []).get(l.met)!).push((rank.get(l.a)! + rank.get(l.b)!) / 2)
    if (!metLane.has(l.met)) metLane.set(l.met, laneOf(rx.get(l.a)!)) // lives with its producer
  }
  for (const [m, arr] of metRanks) {
    const r = arr.reduce((x, y) => x + y, 0) / arr.length
    tmp.set(`M:${m}`, { key: `M:${m}`, kind: 'm', id: m, rank: Math.floor(r) + 0.5, lane: metLane.get(m)! })
  }

  const emap = new Map<string, PlacedEdge>()
  const addEdge = (s: string, t: string, mass: number, spineFlag: boolean) => {
    const k = `${s}>${t}`
    const cur = emap.get(k)
    if (!cur) emap.set(k, { s, t, mass, spine: spineFlag })
    else {
      cur.mass = Math.max(cur.mass, mass)
      cur.spine = cur.spine || spineFlag
    }
  }
  for (const l of links) {
    addEdge(`R:${l.a}`, `M:${l.met}`, l.mass, l.spine)
    addEdge(`M:${l.met}`, `R:${l.b}`, l.mass, l.spine)
  }
  const edges = [...emap.values()]

  // Lanes ordered by mean reaction rank so the flow still reads roughly top-left to bottom-right.
  const laneRanks = new Map<string, number[]>()
  for (const n of tmp.values()) if (n.kind === 'r') (laneRanks.get(n.lane) ?? laneRanks.set(n.lane, []).get(n.lane)!).push(n.rank)
  for (const n of tmp.values()) if (!laneRanks.has(n.lane)) laneRanks.set(n.lane, [n.rank])
  const mean = (a: number[]) => a.reduce((x, y) => x + y, 0) / a.length
  const laneNames = [...laneRanks.keys()].sort((a, b) => mean(laneRanks.get(a)!) - mean(laneRanks.get(b)!) || a.localeCompare(b))

  // Cells = (lane, rank). Nodes stack vertically inside a cell; barycenter sweeps order them.
  const cellKey = (n: Tmp) => `${n.lane}|${n.rank}`
  const cells = new Map<string, Tmp[]>()
  for (const n of [...tmp.values()].sort((a, b) => a.id.localeCompare(b.id))) (cells.get(cellKey(n)) ?? cells.set(cellKey(n), []).get(cellKey(n))!).push(n)
  const laneRows = new Map<string, number>()
  for (const [k, arr] of cells) {
    const lane = k.slice(0, k.lastIndexOf('|'))
    laneRows.set(lane, Math.max(laneRows.get(lane) ?? 1, arr.length))
  }
  const laneBox = new Map<string, Lane>()
  let yCursor = 0
  const hasLabel = opts.arrange === 'subsystem'
  for (const name of laneNames) {
    const h = 2 * PAD + laneRows.get(name)! * ROWH
    laneBox.set(name, { name, top: yCursor, bottom: yCursor + h })
    yCursor += h
  }
  const nbr = new Map<string, string[]>()
  for (const k of tmp.keys()) nbr.set(k, [])
  for (const e of edges) {
    nbr.get(e.s)!.push(e.t)
    nbr.get(e.t)!.push(e.s)
  }
  const y = new Map<string, number>()
  const place = () => {
    for (const arr of cells.values()) {
      const lane = laneBox.get(arr[0].lane)!
      const off = (laneRows.get(arr[0].lane)! - arr.length) / 2
      arr.forEach((n, i) => y.set(n.key, lane.top + PAD + (i + off) * ROWH + ROWH / 2))
    }
  }
  place()
  for (let it = 0; it < 6; it++) {
    for (const arr of cells.values()) {
      const bary = new Map(
        arr.map((n) => {
          const ns = nbr.get(n.key)!.filter((k) => tmp.get(k)!.rank !== n.rank)
          return [n.key, ns.length ? ns.reduce((s, k) => s + y.get(k)!, 0) / ns.length : y.get(n.key)!] as const
        }),
      )
      arr.sort((a, b) => bary.get(a.key)! - bary.get(b.key)! || a.id.localeCompare(b.id))
    }
    place()
  }

  const r0 = Math.min(...[...tmp.values()].map((n) => n.rank))
  const nodes = new Map<string, PlacedNode>()
  for (const n of tmp.values()) {
    const base: PlacedNode = { key: n.key, kind: n.kind, id: n.id, x: (n.rank - r0) * COLW, y: y.get(n.key)! }
    if (n.kind === 'r') {
      base.rxn = rx.get(n.id)
      base.cofactorOnly = !linked.has(n.id)
      base.onSpine = spineNode.has(n.id)
    } else base.met = meta(n.id)
    nodes.set(n.key, base)
  }
  const xs = [...nodes.values()].map((n) => n.x)
  const bbox = { x0: Math.min(...xs) - 110 - (hasLabel ? LABEL_W : 0), x1: Math.max(...xs) + 110, y0: 0, y1: Math.max(yCursor, 1) }
  const lanes = laneNames.map((n) => laneBox.get(n)!)
  const comps = new Set<string>()
  for (const n of nodes.values()) if (n.met?.compartment) comps.add(n.met.compartment)
  const designs = new Set<string>()
  for (const n of nodes.values()) if (n.rxn?.design) designs.add(n.rxn.design)

  return {
    nodes,
    edges,
    lanes,
    bbox,
    maxMass: Math.max(1e-12, ...edges.map((e) => e.mass)),
    shown: [...nodes.values()].filter((n) => n.kind === 'r').length,
    total: rx.size,
    compartments: [...comps].sort(),
    designs: [...designs],
  }
}

function dominantRoute(edges: { source: string; target: string; mass_flow?: number }[], sources: string[], target: string, ids: string[]): string[] {
  const adj = new Map<string, [string, number][]>(ids.map((id) => [id, []]))
  const mx = Math.max(1e-12, ...edges.map((e) => e.mass_flow ?? 0))
  for (const e of edges) adj.get(e.source)!.push([e.target, Math.log(mx / Math.max(e.mass_flow ?? 0, 1e-12))])
  const dist = new Map<string, number>(sources.filter((s) => adj.has(s)).map((s) => [s, 0]))
  const prev = new Map<string, string>()
  const done = new Set<string>()
  for (;;) {
    let u: string | null = null
    let best = Infinity
    for (const [k, v] of dist)
      if (!done.has(k) && v < best) {
        best = v
        u = k
      }
    if (u === null || u === target) break
    done.add(u)
    for (const [v, w] of adj.get(u)!) {
      const nd = best + w
      if (!dist.has(v) || nd < dist.get(v)!) {
        dist.set(v, nd)
        prev.set(v, u)
      }
    }
  }
  if (!dist.has(target)) return []
  const path = [target]
  for (let g = 0; prev.has(path[0]) && g < 5000; g++) path.unshift(prev.get(path[0])!)
  return path
}
