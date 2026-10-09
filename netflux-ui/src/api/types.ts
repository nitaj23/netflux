// Wire types for the NetFlux API (/api/v1). Mirrors src/api/schemas.py; the
// compute endpoints return free-form dicts, so only the fields the UI reads are typed.

export interface ApiErrorBody {
  error?: string
  detail?: string
  hint?: string | null
}

export interface ModelOut {
  id: number
  organism: string
  source_db?: string | null
  bigg_id?: string | null
}

export interface CompoundOut {
  base_id: string
  name?: string | null
  formula?: string | null
  kegg_id?: string | null
  compartments: string[]
  suggested_target?: string | null
}

export interface CompoundOrganismOut {
  model: ModelOut
  has_exchange: boolean
  carries_flow: boolean
  suggested_target?: string | null
}

/** An organism row tagged with the compound id under which that model stores the compound. */
export type OrganismRow = CompoundOrganismOut & { base_id: string }

/** Repo compounds that share a KEGG id but use different id schemes (e.g. `succ` / `M_SUC`). */
export type CompoundGroup = CompoundOut & { base_ids: string[] }

export interface ModelDefaults {
  objective_rxn?: string | null
  objective_rxns: string[]
  direction: string
  carbon_source?: string | null
  carbon_source_options: string[]
  exchange_bounds: Record<string, [number, number]>
  target?: string | null
}

export interface GemRegistration {
  gem_token: string
  summary: Record<string, unknown>
  warnings: string[]
}

export type SourceSpec =
  | { strategy: 'repo'; model_id: number; overrides?: Overrides }
  | { strategy: 'custom'; gem_token: string; overrides?: Overrides }

export interface Overrides {
  custom_exchange_bounds?: Record<string, [number, number]>
  objective_rxn?: string
}

export interface SubgraphNode {
  reaction_id: string
  reaction_name?: string
  flux?: number
  subsystem?: string
  lower_bound?: number
  upper_bound?: number
  gene_rule?: string
  genes?: string[]
}

export interface SubgraphEdge {
  source: string
  target: string
  mass_flow?: number
  metabolites?: string[]
}

export interface Subgraph {
  directed: boolean
  nodes: SubgraphNode[]
  edges: SubgraphEdge[]
}

export interface ReactionAnnotation {
  coupling?: number | null
  per_op?: { v_min?: number | null; v_max?: number | null; cls?: string; design?: string }[]
}

export interface PathwayResponse {
  target_id: string
  reactions: string[]
  n_pathway: number
  n_consistent: number
  n_active: number
  subgraph: Subgraph
  annotation?: { per_reaction?: Record<string, ReactionAnnotation>; v_star_target?: number } | null
  graphml_base64?: string | null
}

export type JobStatus = 'queued' | 'running' | 'completed' | 'failed'

export interface JobSubmit {
  job_id: string
  status: JobStatus
}

export interface ResultTable {
  kind: string
  columns: { key: string; label: string }[]
  rows: Record<string, unknown>[]
}

export interface DesignResult {
  designs: unknown[]
  hypotheses: unknown[]
  table: ResultTable
  hypotheses_table: ResultTable
  csv: string
  hypotheses_csv?: string
  baseline?: Record<string, unknown>
  rejected?: { n: number; reasons: Record<string, number> }
  notes?: string[]
}

export interface JobState {
  job_id: string
  status: JobStatus
  result?: DesignResult | null
  error?: ApiErrorBody | null
}
