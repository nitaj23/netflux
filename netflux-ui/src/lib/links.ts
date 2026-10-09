// External database links for genes and reactions. Repo (BiGG) models resolve to the exact
// BiGG page; uploaded GEMs (unknown id scheme) fall back to an NCBI search, and reactions get
// no link since NCBI has no reaction namespace.

export function geneUrl(gene: string, biggId?: string, organism?: string): string {
  return biggId
    ? `https://bigg.bio/models/${encodeURIComponent(biggId)}/genes/${encodeURIComponent(gene)}`
    : `https://www.ncbi.nlm.nih.gov/gene/?term=${encodeURIComponent(`${gene} ${organism ?? ''}`.trim())}`
}

export function reactionUrl(reaction: string, biggId?: string): string | null {
  return biggId ? `https://bigg.bio/models/${encodeURIComponent(biggId)}/reactions/${encodeURIComponent(reaction)}` : null
}
