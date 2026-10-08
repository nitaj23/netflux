import { apiUrl, get, post, upload } from './client'
import type {
  CompoundOrganismOut,
  OrganismRow,
  CompoundOut,
  GemRegistration,
  JobState,
  JobSubmit,
  ModelDefaults,
  PathwayResponse,
  SourceSpec,
} from './types'

const q = (params: Record<string, string | number>) =>
  new URLSearchParams(Object.entries(params).map(([k, v]) => [k, String(v)])).toString()

export const searchCompounds = (text: string, gemToken: string | null, signal?: AbortSignal) =>
  get<CompoundOut[]>(
    gemToken
      ? `/gems/${encodeURIComponent(gemToken)}/compounds?${q({ q: text, limit: 50 })}`
      : `/compounds?${q({ q: text, limit: 50 })}`,
    signal,
  )

/** Organisms for every id in the group, merged by model (a model that has an exchange wins), flow-carrying first. */
export async function compoundOrganisms(baseIds: string[]): Promise<OrganismRow[]> {
  const lists = await Promise.all(
    baseIds.map(async (id) =>
      (await get<CompoundOrganismOut[]>(`/compounds/${encodeURIComponent(id)}/organisms`)).map((r) => ({ ...r, base_id: id })),
    ),
  )
  const byModel = new Map<number, OrganismRow>()
  for (const r of lists.flat()) {
    const prev = byModel.get(r.model.id)
    const score = (x: OrganismRow) => Number(x.has_exchange) * 2 + Number(x.carries_flow)
    if (!prev || score(r) > score(prev)) byModel.set(r.model.id, r)
  }
  return [...byModel.values()].sort(
    (a, b) =>
      Number(b.carries_flow) - Number(a.carries_flow) ||
      Number(b.has_exchange) - Number(a.has_exchange) ||
      a.model.organism.localeCompare(b.model.organism),
  )
}

export const fetchDefaults = (compound: string, modelId: number | null, gemToken: string | null) =>
  get<ModelDefaults>(
    gemToken
      ? `/gems/${encodeURIComponent(gemToken)}/defaults?${q({ compound })}`
      : `/models/${modelId}/defaults?${q({ compound })}`,
  )

export const uploadGem = (file: File) => {
  const form = new FormData()
  form.append('file', file)
  return upload<GemRegistration>('/gems/upload', form)
}

/** Target pinned at this fraction of its optimum for the FVA annotation (same as the corpus sweep). */
const FVA_FRACTION = 0.5

export const extractPathway = (source: SourceSpec, target: string) =>
  post<PathwayResponse>('/pathways', { source, target, include_graphml: true, fva_operating_points: [FVA_FRACTION] })

export interface DesignBody {
  source: SourceSpec
  target: string
  biomass: string
  substrate?: string
}

export const submitDesign = (body: DesignBody) => post<JobSubmit>('/designs', body)

export const fetchJob = (jobId: string) => get<JobState>(`/jobs/${encodeURIComponent(jobId)}`)

/** Direct-download URL for the base GEM behind a run (repo file or uploaded GEM), else null. */
export const modelDownloadUrl = (source: SourceSpec): string | null =>
  source.strategy === 'repo'
    ? apiUrl(`/models/${source.model_id}/download`)
    : source.strategy === 'custom'
      ? apiUrl(`/gems/${encodeURIComponent(source.gem_token)}/download`)
      : null
