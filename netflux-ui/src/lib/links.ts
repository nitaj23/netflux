// External database links for genes and reactions. Models hosted in BiGG resolve to the exact
// BiGG page; models from other sources (and uploaded GEMs) fall back to an NCBI gene search, and
// reactions get no link there since NCBI has no reaction namespace.

// Repo models NOT in BiGG (their ids/genes won't resolve there). Matched by prefix so variants
// like "iCTH669_w_GLGC" are covered.
const NON_BIGG_MODELS = ['iCN1361', 'iCTH669', 'iJN1026', 'iLB1027']

/** The BiGG model id to link against, or undefined when the model isn't in BiGG. */
function biggModel(biggId?: string): string | undefined {
  if (!biggId) return undefined
  return NON_BIGG_MODELS.some((m) => biggId.startsWith(m)) ? undefined : biggId
}

export function geneUrl(gene: string, biggId?: string, organism?: string): string {
  const model = biggModel(biggId)
  return model
    ? `https://bigg.bio/models/${encodeURIComponent(model)}/genes/${encodeURIComponent(gene)}`
    : `https://www.ncbi.nlm.nih.gov/gene/?term=${encodeURIComponent(`${gene} ${organism ?? ''}`.trim())}`
}

export function reactionUrl(reaction: string, biggId?: string): string | null {
  const model = biggModel(biggId)
  return model ? `https://bigg.bio/models/${encodeURIComponent(model)}/reactions/${encodeURIComponent(reaction)}` : null
}
