import type { ApiErrorBody } from './types'

const BASE = (import.meta.env.VITE_API_BASE as string | undefined) ?? '/api/v1'

/** Absolute URL for a backend path — for direct downloads (GET with a file body). */
export const apiUrl = (path: string) => `${BASE}${path}`

/** Every backend failure arrives as {error, detail, hint}; this carries it intact. */
export class ApiError extends Error {
  readonly status: number
  readonly type: string
  readonly hint?: string | null

  constructor(status: number, body: ApiErrorBody) {
    super(body.detail ?? body.error ?? `Request failed (${status})`)
    this.status = status
    this.type = body.error ?? 'error'
    this.hint = body.hint
  }
}

function toBody(raw: unknown, status: number): ApiErrorBody {
  if (raw && typeof raw === 'object') {
    const o = raw as Record<string, unknown>
    // FastAPI validation errors put a list in `detail`.
    const detail = Array.isArray(o.detail)
      ? o.detail.map((d) => (d as { msg?: string }).msg ?? String(d)).join('; ')
      : (o.detail as string | undefined)
    return { error: o.error as string | undefined, detail, hint: o.hint as string | undefined }
  }
  return { detail: `Request failed (${status})` }
}

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  let res: Response
  try {
    res = await fetch(`${BASE}${path}`, init)
  } catch {
    throw new ApiError(0, {
      error: 'network',
      detail: 'Could not reach the NetFlux server.',
      hint: 'Check that the backend is running and reachable.',
    })
  }
  const text = await res.text()
  let data: unknown = null
  try {
    data = text ? JSON.parse(text) : null
  } catch {
    /* non-JSON body */
  }
  if (!res.ok) throw new ApiError(res.status, toBody(data, res.status))
  return data as T
}

export const get = <T>(path: string, signal?: AbortSignal) => request<T>(path, { signal })

export const post = <T>(path: string, body: unknown) =>
  request<T>(path, {
    method: 'POST',
    headers: { 'content-type': 'application/json' },
    body: JSON.stringify(body),
  })

export const upload = <T>(path: string, form: FormData) =>
  request<T>(path, { method: 'POST', body: form })
