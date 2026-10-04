# EUCAST 16.1 Breakpoint Transcription Audit

Status: candidate extraction only; independent verification required  
Artifact SHA-256: `ee5709eb5b7c9beb1c9ba8182a4e7032e71963aaca30e618b9f76dc53f6bab73`

The official workbook was downloaded on 2026-08-15 from the EUCAST clinical
breakpoint page. Its XLSX archive has 179 members and passes ZIP integrity
testing. The source artifact and acquisition manifest live under
`data/reference/` and are intentionally ignored by Git.

## Structured Cell Extraction

The following values were read from the `Enterobacterales` worksheet through
OOXML cell references. They are evidence for review, not an approved rules file.

| Agent and indication | MIC S | MIC R | MIC ATU | Disk content | Zone S | Zone R | Zone ATU | Cells |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Ceftriaxone, indications other than meningitis | <=1 mg/L | >2 mg/L | — | 30 ug | >=27 mm | <24 mm | — | A55:H55 |
| Ceftriaxone, meningitis | <=1 mg/L | >1 mg/L | — | 30 ug | >=27 mm | <27 mm | — | A56:H56 |
| Ciprofloxacin, indications other than meningitis | <=0.25 mg/L | >0.5 mg/L | 0.5 mg/L | 5 ug | >=25 mm | <22 mm | 22–24 mm | A82:H82 |
| Ciprofloxacin, meningitis | <=0.125 mg/L | >0.125 mg/L | — | method note | method note | method note | — | A83:H83 |
| Gentamicin, systemic infections | bracketed 2 mg/L | bracketed 2 mg/L | — | 10 ug | bracketed 17 mm | bracketed 17 mm | — | A97:H97 |
| Gentamicin, urinary-tract infections | <=2 mg/L | >2 mg/L | — | 10 ug | >=17 mm | <17 mm | — | A98:H98 |

Headers at B/C and F/G explicitly define `S <=`, `R >`, `S >=`, and `R <`.
Bracketed systemic gentamicin values invoke separate EUCAST guidance and must
not be converted into ordinary breakpoints without that policy. Ciprofloxacin
and ceftriaxone also have indication-specific rows. Therefore:

- `clinical_indication` is a required, exact breakpoint-rule constraint;
- specimen source is not accepted as a clinical-indication proxy;
- generic drug rules are prohibited when multiple indication rows exist;
- disk rules require exact disk content and relevant method notes;
- two reviewers must compare candidate rows, notes, and styles against the
  frozen workbook before `data/reference/EUCAST_v16.1_rules.csv` is approved.

## Deterministic OOXML cross-check

`scripts/audit_eucast_workbook.py` independently authenticates the frozen
workbook and acquisition manifest, then records 135 targeted cells from the
`Enterobacterales` sheet with their OOXML style IDs, merged ranges, rich-text
runs, superscript markers, comments, ATU columns, and applicable section-note
anchors. The generated audit is
`data/interim/eucast_v16_1_ooxml_audit.json` (SHA-256
`d2c8eab23a321ec7b20e10f0d128c0f0b5482c720b4be5281241992688a36675`).

The cross-check confirmed that:

- ciprofloxacin for indications other than meningitis has MIC ATU `0.5 mg/L`
  at `D82` and zone ATU `22–24 mm` at `H82`;
- the `2` marker on the meningitis row and the `B` markers in `F83:G83` are
  superscripts linked to note `2/B`, which requires an MIC test or inference
  from a pefloxacin 5 µg screen rather than an ordinary ciprofloxacin disk
  breakpoint;
- the systemic gentamicin values are truly bracketed, and their superscript
  `1/A` markers link to EUCAST bracket guidance in `I95`;
- the urinary gentamicin row is a separate, unbracketed indication row.

This deterministic audit reduces transcription risk but is not the required
independent human review. Approval remains blocked until a second reviewer
checks the cells, section notes, ATU handling, rich-text markers, and visual
styles in the frozen workbook.

Official source page:
https://www.eucast.org/bacteria/clinical-breakpoints-and-interpretation/clinical-breakpoint-tables/
