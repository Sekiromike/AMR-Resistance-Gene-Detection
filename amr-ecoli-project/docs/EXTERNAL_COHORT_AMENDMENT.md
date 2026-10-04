# External Cohort Amendment

Amendment: 002
Date: 2026-09-20
Status: ADOPTED
Decision owner: Rushath Rajeev
Depends on: [`ENDPOINT_AMENDMENT.md`](ENDPOINT_AMENDMENT.md) (amendment 001)

> **Partly superseded by amendment 003**,
> [`EXTERNAL_THRESHOLD_AMENDMENT.md`](EXTERNAL_THRESHOLD_AMENDMENT.md).
> The exclusion rule below uses the near-neighbour threshold, which removes
> lineage membership rather than clone-level leakage. The locked external set
> is now the 1,888 isolates surviving the putative-duplicate rule (1,555 with
> all three MICs); the 536/454 set described here is retained as a nested
> lineage-novel stratum. Everything else in this amendment stands.

## Decision

The external evaluation design changes from "all 3,159 JARBS *E. coli*" to two
prespecified, sequence-only strata plus a new internal holdout:

| Set | Definition | n | Role |
|---|---|---:|---|
| **External-locked** | JARBS isolates with zero genomic collision against development, with all three target MICs | **454** | the one locked external evaluation |
| **External-stress** | JARBS isolates that collide with development | 2,623 | reported separately, never a headline claim |
| **Internal-holdout** | genomically disjoint components carved from the 10,334 development isolates | pending compute | second OOD axis |

## Why

Unity comparison `63772662` completed 2026-08-30 with exit `0:0` after 8h16m,
running all 3,159 frozen external queries against 32,616,080 candidate
comparisons. It found **578,339 collision pairs across 166 components**, every
pair unresolved:

| Collision class | Pairs |
|---|---:|
| `NEAR_NEIGHBOR` (ANI ≥ 99.9%, min reciprocal AF ≥ 0.9) | 556,021 |
| `PUTATIVE_GENOMIC_DUPLICATE` (ANI ≥ 99.99%, min reciprocal AF ≥ 0.95) | 22,318 |

**2,623 of 3,159 JARBS isolates (83.0%) collide with development**; 1,271 have
at least one putative genomic duplicate and 1,352 are near-neighbour only.
3,675 development isolates are involved. Observed ANI spans 99.90–100.00%.

Exact identifier overlap had been **zero** across BioSample, SRA run,
registered assembly, and public isolate ID. Identifier-based disjointness is
therefore demonstrably insufficient for *E. coli*, where globally disseminated
lineages place near-twins on different continents and in different
surveillance programmes. Evaluating on the full JARBS set would have reported
"external validation" on a cohort that is mostly not external.

Provenance: requests SHA-256
`8b2d2689ad3be1b8698783a94afcf3a0b9ec33f199c98353a1e90e28cee9408f`, contract
SHA-256 `8dfe0265a76ba9e75c920e0037773fd4197e10d993decf85398582b1d6392386`,
comparisons SHA-256
`9ce30f8f29e2f5b8a0227b64c9c2ba2c1006e57345b8c1d3ac72556022ab72ff`,
collisions SHA-256
`bd50d587ac88d3a4588a5287186ad5782e4f4f059a1bc0bcb98aeb976eecce3a`.

## Selection integrity

The 536 collision-free isolates were selected by a **sequence-only,
phenotype-blind** criterion computed before any external phenotype was read.
The 454 are those with all three target MICs; the remaining 82 have none. No
model, feature, threshold, or epoch was involved, and `external_phenotypes_read`
remains false in the frozen contract.

## Limitations that must be reported

1. **Selection-induced distinctness.** The 454 are by construction the JARBS
   isolates *least* similar to the development pool. They are a harder
   out-of-distribution test than a random JARBS sample, and they are not
   representative of JARBS, of Japanese surveillance, or of any routine
   population. Performance on them is a lower bound under strong shift, not an
   estimate of deployed accuracy.
2. **Resistance enrichment persists.** JARBS sampling deliberately enriched for
   third-generation-cephalosporin resistance or reduced meropenem
   susceptibility. Prevalence in the 454 does not estimate routine prevalence,
   so predictive values must be reported prevalence-adjusted.
3. **Reduced size.** 454 isolates × 3 drugs bounds achievable confidence-interval
   width. A formal precision calculation is required per drug before the locked
   evaluation, and any drug failing it is reported as underpowered rather than
   quietly included.
4. **Absent stratifiers remain absent.** JARBS publishes no AST testing date,
   hospital, or clinical indication. Leave-site-out evaluation is impossible on
   this cohort and must not be simulated from geography or collection date.

## Internal holdout

`genomic_cluster` is already a required cohort column, so within-development
clustering is a prerequisite for grouped cross-validation regardless of this
amendment. The same frozen near-neighbour rule (ANI ≥ 99.9%, min reciprocal
AF ≥ 0.9) is reused; **no new threshold is introduced**. Connected components
over development-versus-development comparisons define `genomic_cluster`, and
the internal holdout is drawn as whole components so no cluster crosses a split
boundary.

Both the external-locked set and the internal holdout stay sealed from feature
selection, model choice, calibration, thresholding, and early stopping.

## Not resolved by this amendment

- The 166 cross-cohort components remain `UNRESOLVED`. This amendment excludes
  them from the locked external set rather than adjudicating them.
- JARBS still supplies no AST date, site, or indication, so the categorical
  secondary endpoint stays blocked under amendment 001.
- The internal holdout size is unknown until development clustering completes.
