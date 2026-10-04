# Development/JARBS genomic-disjointness protocol

## Evidence boundary

This gate uses sequence identity only. It does not read JARBS MIC values, tune a
model, choose a decision threshold from external outcomes, or create clinical or
causal evidence. JARBS remains a conditional external challenge until every
sequence query is complete and every cross-cohort collision is resolved.

## Frozen population and sequence precedence

- Development: all 10,334 NCBI AST parent isolates from enrichment snapshot
  `20260821T012604Z-v1`.
- External candidate: all 3,159 *E. coli* isolates in JARBS snapshot
  `20260825T021722Z-v1`.
- A registered NCBI assembly is used when available (9,986 development
  isolates). Otherwise, the sole linked short-read run is assembled by the
  frozen read workflow (348 development and 3,159 external isolates).
- No assembly is substituted across BioSamples. Downloaded reads and assemblies
  must be checksumed; failed identity or assembly-QC checks are explicit
  exclusions, never silent replacements.

## Prespecified relatedness rule

All eligible assemblies are compared with skani 0.3.1. Candidate discovery uses
an approximate ANI screen of 95% and a deliberately broader minimum reciprocal
aligned-fraction floor of 0.50. The final leakage rule uses the minimum of the
reciprocal aligned fractions:

- putative genomic duplicate: reported ANI >= 99.99% and minimum reciprocal
  aligned fraction >= 0.95;
- near-neighbor collision: ANI >= 99.9% and minimum reciprocal aligned fraction
  >= 0.90.

The higher category is explicitly putative because rounded ANI is not proof of
identical genomes or transmission. These are conservative study-governance
thresholds for leakage control, not a claim that they universally define
strains, transmission, or species. A
sensitivity ledger at 99.5%, 99.9%, and 99.99% ANI will be retained without
consulting external phenotypes.

The comparison uses skani's default learned-ANI adjustment, compression factor
125, and marker compression factor 1000. Candidate search is executed with
`-s 95 --both-min-af 50`; reported percentage values are normalized to fractions
for the ledger. Thresholds are applied to the two-decimal values emitted by the
pinned skani release, which makes the high category deliberately conservative.

Any cross-cohort edge joins both genomes into a collision component. Each
component must receive a documented, sequence-only resolution. Until the
collision ledger is empty or every edge is resolved under a prespecified cohort
policy, the external lock fails closed.

## Completion evidence

The gate requires all of the following:

1. frozen request table and input/output SHA-256 hashes;
2. raw-read and assembly accession completeness;
3. read/assembly checksums, assembly commands, environment lock, and QC;
4. a completion marker for every one of the 3,159 external queries;
5. normalized ANI/aligned-fraction comparisons and the collision-component
   ledger;
6. tool versions, Slurm job metadata, exit status, and restart-safe logs.

The planning and evaluation contract is implemented in
`scripts/build_genomic_disjointness.py`. After acquisition arrays leave the
queue, submit `hpc/unity/disjointness-audit.sbatch`; it fails closed unless
every frozen request has a completion marker and every declared artifact hash
verifies. Retained failed-task partial directories are reported separately and
are never deleted automatically.

After and only after that audit exits successfully,
`hpc/unity/disjointness-compare.sbatch` builds a skani database from all
development assemblies and searches every external assembly with the frozen
`-s 95 --both-min-af 50` screen. The job maps paths back to frozen sequence IDs,
normalizes skani percentage aligned fractions to fractions, authenticates every
comparison artifact, and writes the unresolved collision-component ledger.
At comparison start it also requires the successful audit to name the same
assembly root and frozen request hash, then re-hashes every requested assembly
and declared provenance artifact. Any post-audit change fails before sketching.

Before submitting that production comparison, run
`hpc/unity/disjointness-search-canary.sbatch` on one already completed frozen
development/external pair. This validates the exact skani 0.3.1 sketch/search
path orientation and output schema only; it is not cohort evidence and does not
unlock the external cohort.

## Bounded Unity pilot

A one-development/one-external pilot completed on 2026-08-25 with NCBI
Datasets 18.29.1, SRA Toolkit 3.2.1, Shovill 1.1.0, and skani 0.3.1. The
external SRA archive MD5 exactly matched its frozen source value, paired reads
were extracted, and Shovill produced a checksumed assembly. The job completed
in 6 minutes 47 seconds with 4.75 GB peak resident memory.

The pilot pair (`GCA_000770275.1` and `DRR385905`) produced reported ANI
99.99%, development aligned fraction 98.21%, and external aligned fraction
99.68%. It is therefore an unresolved putative-genomic-duplicate collision
under the prespecified leakage rule. This immediately confirms that exact
identifier disjointness alone was insufficient; it does not estimate the total
number of collisions. JARBS remains conditional until the complete 13,493-
isolate comparison and resolution ledger are finished.
