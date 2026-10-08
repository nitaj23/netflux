import { useState } from 'react'
import { Modal } from './Modal'

export function Header() {
  const [modal, setModal] = useState<'about' | 'help' | null>(null)
  return (
    <header className="topbar">
      <h1>NetFlux</h1>
      <nav>
        <button className="nav-link" onClick={() => setModal('about')}>
          About
        </button>
        <button className="nav-link" onClick={() => setModal('help')}>
          Help
        </button>
        <a className="nav-link" href="https://2026.igem.wiki/mit-mahe/" target="_blank" rel="noreferrer">
          MIT-MAHE 2026
        </a>
      </nav>
      {modal === 'about' && (
        <Modal title="About NetFlux" onClose={() => setModal(null)}>
          <p>
            NetFlux extracts bioproduction pathways from genome-scale metabolic models and ranks
            strain-design interventions that couple a target compound to growth.
          </p>
          <p>
            Pick a compound and an organism from the repository, or upload your own COBRA JSON / SBML
            model. All computation is ephemeral: nothing you upload is stored.
          </p>
        </Modal>
      )}
      {modal === 'help' && (
        <Modal title="Using NetFlux" onClose={() => setModal(null)}>
          <ol>
            <li>Type a compound and choose it from the list.</li>
            <li>Choose an organism. Models that carry flow for the compound are listed first.</li>
            <li>Optionally adjust the substrate, its flux bounds, or the objective reaction.</li>
            <li>
              Press <b>Generate</b> for a Subgraph tab. Each Generate adds a new tab.
            </li>
            <li>
              Open <b>Interventions</b> to rank strain designs for the selected Subgraph. This runs in the
              background and can take a few minutes.
            </li>
            <li>
              Use <b>download as</b> to save the Subgraph (.json, .graphml) or the designs (.csv).
            </li>
          </ol>
          <p>
            To use your own model, press <b>Upload</b> beside Generate (.json or .xml, up to 50 MB). Uploaded
            models are kept in memory for 30 minutes.
          </p>
        </Modal>
      )}
    </header>
  )
}
