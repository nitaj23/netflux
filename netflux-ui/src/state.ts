import type { ApiError } from './api/client'
import type { PathwayResponse, SourceSpec } from './api/types'

/** Everything a Generate click froze: the request for the Subgraph and for the design job. */
export interface RunSnapshot {
  source: SourceSpec
  target: string
  biomass: string
  substrate?: string
  /** Feedstock uptake exchange the map is drawn from (Substrate, else the model default). */
  carbonSource?: string
  label: string
  /** BiGG model id for repo models (e.g. 'iML1515'); enables exact BiGG gene links. */
  biggId?: string
  /** Organism name, for the gene-search fallback link. */
  organism?: string
}

export interface SubgraphTab {
  id: number
  title: string
  snapshot: RunSnapshot
  status: 'loading' | 'done' | 'error'
  data?: PathwayResponse
  error?: ApiError
  jobId?: string
  jobStartedAt?: number
  jobSubmitError?: ApiError
}

export type ActiveKey = number | 'interventions'
