import { useId, useRef, useState } from 'react'

export interface Option {
  value: string
  label: string
  sub?: string
}

interface Props {
  value: string
  onChange: (text: string) => void
  onSelect: (o: Option) => void
  options: Option[]
  loading?: boolean
  placeholder?: string
  disabled?: boolean
  id?: string
  invalid?: boolean
}

/** Text input with a filterable suggestion list. Free text is allowed; options only suggest. */
export function Combobox({ value, onChange, onSelect, options, loading, placeholder, disabled, id, invalid }: Props) {
  const [open, setOpen] = useState(false)
  const [hi, setHi] = useState(0)
  const listId = useId()
  const wrap = useRef<HTMLDivElement>(null)

  const pick = (o: Option) => {
    onSelect(o)
    setOpen(false)
  }

  return (
    <div
      className="combo"
      ref={wrap}
      onBlur={(e) => {
        if (!wrap.current?.contains(e.relatedTarget as Node)) setOpen(false)
      }}
    >
      <input
        id={id}
        className="field"
        role="combobox"
        aria-expanded={open}
        aria-controls={listId}
        aria-autocomplete="list"
        aria-invalid={invalid || undefined}
        autoComplete="off"
        spellCheck={false}
        disabled={disabled}
        placeholder={placeholder}
        value={value}
        onFocus={() => setOpen(true)}
        onChange={(e) => {
          onChange(e.target.value)
          setOpen(true)
          setHi(0)
        }}
        onKeyDown={(e) => {
          if (e.key === 'ArrowDown') {
            e.preventDefault()
            setOpen(true)
            setHi((h) => Math.min(h + 1, options.length - 1))
          } else if (e.key === 'ArrowUp') {
            e.preventDefault()
            setHi((h) => Math.max(h - 1, 0))
          } else if (e.key === 'Enter' && open && options[hi]) {
            e.preventDefault()
            pick(options[hi])
          } else if (e.key === 'Escape') {
            setOpen(false)
          }
        }}
      />
      {open && (loading || options.length > 0) && (
        <ul className="combo-list" id={listId} role="listbox">
          {loading && options.length === 0 && <li className="combo-empty">Searching…</li>}
          {options.map((o, i) => (
            <li
              key={o.value}
              role="option"
              aria-selected={i === hi}
              className={i === hi ? 'combo-opt active' : 'combo-opt'}
              onMouseEnter={() => setHi(i)}
              onMouseDown={(e) => e.preventDefault()}
              onClick={() => pick(o)}
            >
              <span>{o.label}</span>
              {o.sub && <small>{o.sub}</small>}
            </li>
          ))}
        </ul>
      )}
    </div>
  )
}
