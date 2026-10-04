"""Cohort-driven assembly staging, QC, MLST, and AMRFinderPlus rules.

Integration contract:
  1. Include this file from the main Snakefile after loading study config.
  2. Add ``GENOMICS_FINAL_MANIFEST`` to the main ``rule all`` input.
  3. Define the documented ``genomics`` config keys before executing the rule.

The checkpoint is intentional: assembly accessions come from independently
interpreted AST rows rather than from a hand-maintained list.
"""
import csv


GENOMICS_CONFIG = config.get("genomics", {})
GENOMICS_ROOT = GENOMICS_CONFIG.get("root", "data/genomics")
GENOMICS_RESULTS = GENOMICS_CONFIG.get("results_root", "results/genomics")
GENOMICS_REFERENCES = GENOMICS_CONFIG.get("reference_root", "data/reference/genomics")
GENOMICS_COHORT = config.get("outputs", {}).get(
    "interpreted_ast", "data/interim/ast_interpreted.csv"
)
GENOMICS_AMR_DB_VERSION = GENOMICS_CONFIG.get(
    "amrfinder_database_version", "REQUIRED_BEFORE_COHORT_FREEZE"
)
GENOMICS_AMR_DB = (
    f"{GENOMICS_REFERENCES}/amrfinderplus/{GENOMICS_AMR_DB_VERSION}"
)
GENOMICS_QC = GENOMICS_CONFIG.get("assembly_qc", {})

GENOMICS_ACCESSIONS = f"{GENOMICS_ROOT}/requests/assembly_accessions.txt"
GENOMICS_REQUESTS = f"{GENOMICS_ROOT}/requests/assembly_mapping.tsv"
GENOMICS_REQUEST_EXCLUSIONS = f"{GENOMICS_ROOT}/requests/exclusions.tsv"
GENOMICS_REQUEST_MANIFEST = f"{GENOMICS_ROOT}/requests/manifest.json"
GENOMICS_PACKAGE_ZIP = f"{GENOMICS_ROOT}/source/ncbi_genomes.zip"
GENOMICS_PACKAGE_DIR = f"{GENOMICS_ROOT}/source/ncbi_genomes"
GENOMICS_PACKAGE_SHA256 = f"{GENOMICS_ROOT}/source/ncbi_genomes.zip.sha256"
GENOMICS_STAGED_DIR = f"{GENOMICS_ROOT}/staged"
GENOMICS_STAGE_TABLE = f"{GENOMICS_ROOT}/manifests/staging.tsv"
GENOMICS_STAGE_EXCLUSIONS = f"{GENOMICS_ROOT}/manifests/staging_exclusions.tsv"
GENOMICS_STAGE_MANIFEST = f"{GENOMICS_ROOT}/manifests/staging.json"
GENOMICS_QUAST_ROOT = f"{GENOMICS_RESULTS}/quast"
GENOMICS_MLST_ROOT = f"{GENOMICS_RESULTS}/mlst"
GENOMICS_AMR_ROOT = f"{GENOMICS_RESULTS}/amrfinderplus"
GENOMICS_TOOL_VERSIONS = f"{GENOMICS_RESULTS}/provenance/tool_versions.txt"
# quast, mlst and amrfinder live in separate environments because their perl,
# zlib and python_abi requirements are mutually unsatisfiable, so each records
# its own version fragment and a combine step concatenates them.
GENOMICS_TOOL_VERSION_FRAGMENTS = {
    "qc": f"{GENOMICS_RESULTS}/provenance/tool_versions.qc.txt",
    "mlst": f"{GENOMICS_RESULTS}/provenance/tool_versions.mlst.txt",
    "amrfinder": f"{GENOMICS_RESULTS}/provenance/tool_versions.amrfinder.txt",
}
GENOMICS_MLST_INFO = f"{GENOMICS_RESULTS}/provenance/mlst_scheme_info.tsv"
GENOMICS_AMR_DB_METADATA = f"{GENOMICS_RESULTS}/provenance/amrfinder_database.txt"
GENOMICS_AMR_DB_CHECKSUMS = f"{GENOMICS_RESULTS}/provenance/amrfinder_database.sha256"
GENOMICS_FINAL_TABLE = f"{GENOMICS_RESULTS}/genome_manifest.tsv"
GENOMICS_FINAL_EXCLUSIONS = f"{GENOMICS_RESULTS}/genome_exclusions.tsv"
GENOMICS_FINAL_MANIFEST = f"{GENOMICS_RESULTS}/genome_manifest.json"


rule genomics_build_requests:
    input:
        cohort=GENOMICS_COHORT,
    output:
        accessions=GENOMICS_ACCESSIONS,
        mapping=GENOMICS_REQUESTS,
        exclusions=GENOMICS_REQUEST_EXCLUSIONS,
        manifest=GENOMICS_REQUEST_MANIFEST,
    threads: 1
    resources:
        mem_mb=2000,
        runtime=15,
    conda:
        "../envs/controller.yaml"
    shell:
        "python scripts/build_genome_manifest.py request "
        "--cohort {input.cohort:q} --accessions {output.accessions:q} "
        "--mapping {output.mapping:q} --exclusions {output.exclusions:q} "
        "--manifest {output.manifest:q}"


rule genomics_download_ncbi_package:
    input:
        accessions=GENOMICS_ACCESSIONS,
        request_manifest=GENOMICS_REQUEST_MANIFEST,
    output:
        archive=GENOMICS_PACKAGE_ZIP,
        package=directory(GENOMICS_PACKAGE_DIR),
        checksum=GENOMICS_PACKAGE_SHA256,
    threads: 1
    resources:
        mem_mb=4000,
        runtime=360,
    conda:
        "../envs/genomics-qc.yaml"
    shell:
        r"""
        datasets download genome accession \
          --inputfile {input.accessions:q} \
          --include genome \
          --dehydrated \
          --filename {output.archive:q}
        mkdir -p {output.package:q}
        unzip -q {output.archive:q} -d {output.package:q}
        datasets rehydrate --directory {output.package:q}
        sha256sum {output.archive:q} > {output.checksum:q}
        """


checkpoint genomics_stage_ncbi_package:
    input:
        requests=GENOMICS_REQUESTS,
        package=GENOMICS_PACKAGE_DIR,
        package_checksum=GENOMICS_PACKAGE_SHA256,
    output:
        staged=directory(GENOMICS_STAGED_DIR),
        table=GENOMICS_STAGE_TABLE,
        exclusions=GENOMICS_STAGE_EXCLUSIONS,
        manifest=GENOMICS_STAGE_MANIFEST,
    threads: 1
    resources:
        mem_mb=4000,
        runtime=180,
    conda:
        "../envs/genomics-qc.yaml"
    shell:
        "python scripts/build_genome_manifest.py stage "
        "--requests {input.requests:q} --package-dir {input.package:q} "
        "--stage-dir {output.staged:q} --table {output.table:q} "
        "--exclusions {output.exclusions:q} --manifest {output.manifest:q}"


def _genomics_staged_rows(wildcards):
    checkpoint_output = checkpoints.genomics_stage_ncbi_package.get().output.table
    with open(checkpoint_output, encoding="utf-8", newline="") as handle:
        return [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["stage_status"] == "STAGED"
        ]


def _genomics_staged_accessions(wildcards):
    return [row["assembly_accession"] for row in _genomics_staged_rows(wildcards)]


def _genomics_staged_fasta(wildcards):
    matches = [
        row["staged_fasta"]
        for row in _genomics_staged_rows(wildcards)
        if row["assembly_accession"] == wildcards.assembly_accession
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one staged FASTA for {wildcards.assembly_accession}; found {len(matches)}"
        )
    return matches[0]


def _genomics_analysis_dirs(wildcards, root):
    return expand(
        f"{root}/{{assembly_accession}}",
        assembly_accession=_genomics_staged_accessions(wildcards),
    )


rule genomics_freeze_amrfinder_database:
    output:
        database=directory(GENOMICS_AMR_DB),
        metadata=GENOMICS_AMR_DB_METADATA,
        checksums=GENOMICS_AMR_DB_CHECKSUMS,
    params:
        expected_version=GENOMICS_AMR_DB_VERSION,
    threads: 2
    resources:
        mem_mb=8000,
        runtime=180,
    conda:
        "../envs/genomics-amrfinder.yaml"
    shell:
        r"""
        if [ {params.expected_version:q} = REQUIRED_BEFORE_COHORT_FREEZE ]; then
          echo "Set genomics.amrfinder_database_version to an archived NCBI database release." >&2
          exit 2
        fi
        mkdir -p {output.database:q}
        amrfinder_update --force_update --database {output.database:q}
        source_db=$(readlink -f {output.database:q}/latest)
        observed_version=$(head -n 1 "$source_db/version.txt" | tr -d '\r\n')
        if [ "$observed_version" != {params.expected_version:q} ]; then
          echo "AMRFinderPlus DB mismatch: expected {params.expected_version}, observed $observed_version" >&2
          exit 2
        fi
        amrfinder --database "$source_db" --database_version > {output.metadata:q}
        (cd {output.database:q} && find . -type f -print0 | sort -z | xargs -0 sha256sum) \
          > {output.checksums:q}
        """


rule genomics_tool_versions_qc:
    input:
        environment="workflow/envs/genomics-qc.yaml",
    output:
        fragment=GENOMICS_TOOL_VERSION_FRAGMENTS["qc"],
    threads: 1
    resources:
        mem_mb=2000,
        runtime=15,
    conda:
        "../envs/genomics-qc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.fragment:q})
        printf 'ncbi-datasets-cli\t' > {output.fragment:q}
        datasets --version >> {output.fragment:q} 2>&1
        printf 'quast\t' >> {output.fragment:q}
        quast.py --version >> {output.fragment:q} 2>&1
        printf '\n[conda-explicit:qc]\n' >> {output.fragment:q}
        conda list --explicit >> {output.fragment:q} 2>&1 || \
          printf 'unavailable; retain the Snakemake Conda archive\n' >> {output.fragment:q}
        """


rule genomics_tool_versions_mlst:
    input:
        environment="workflow/envs/genomics-mlst.yaml",
    output:
        fragment=GENOMICS_TOOL_VERSION_FRAGMENTS["mlst"],
        mlst_info=GENOMICS_MLST_INFO,
    threads: 1
    resources:
        mem_mb=2000,
        runtime=15,
    conda:
        "../envs/genomics-mlst.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.fragment:q})
        printf 'mlst\t' > {output.fragment:q}
        mlst --version >> {output.fragment:q} 2>&1
        printf '\n[conda-explicit:mlst]\n' >> {output.fragment:q}
        conda list --explicit >> {output.fragment:q} 2>&1 || \
          printf 'unavailable; retain the Snakemake Conda archive\n' >> {output.fragment:q}
        mlst --info > {output.mlst_info:q}
        """


rule genomics_tool_versions_amrfinder:
    input:
        environment="workflow/envs/genomics-amrfinder.yaml",
    output:
        fragment=GENOMICS_TOOL_VERSION_FRAGMENTS["amrfinder"],
    threads: 1
    resources:
        mem_mb=2000,
        runtime=15,
    conda:
        "../envs/genomics-amrfinder.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.fragment:q})
        printf 'amrfinderplus\t' > {output.fragment:q}
        amrfinder --version >> {output.fragment:q} 2>&1
        printf '\n[conda-explicit:amrfinder]\n' >> {output.fragment:q}
        conda list --explicit >> {output.fragment:q} 2>&1 || \
          printf 'unavailable; retain the Snakemake Conda archive\n' >> {output.fragment:q}
        """


rule genomics_tool_metadata:
    input:
        qc=GENOMICS_TOOL_VERSION_FRAGMENTS["qc"],
        mlst=GENOMICS_TOOL_VERSION_FRAGMENTS["mlst"],
        amrfinder=GENOMICS_TOOL_VERSION_FRAGMENTS["amrfinder"],
    output:
        versions=GENOMICS_TOOL_VERSIONS,
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,
    shell:
        "cat {input.qc:q} {input.mlst:q} {input.amrfinder:q} > {output.versions:q}"


rule genomics_quast_assembly:
    input:
        fasta=_genomics_staged_fasta,
    output:
        result=directory(f"{GENOMICS_QUAST_ROOT}/{{assembly_accession}}"),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
    conda:
        "../envs/genomics-qc.yaml"
    shell:
        r"""
        mkdir -p {output.result:q}
        set +e
        quast.py --threads {threads} --min-contig 0 --output-dir {output.result:q} \
          {input.fasta:q} > {output.result:q}/command.stdout 2> {output.result:q}/command.stderr
        rc=$?
        set -e
        printf '%s\n' "$rc" > {output.result:q}/exit_code.txt
        """


rule genomics_mlst_assembly:
    input:
        fasta=_genomics_staged_fasta,
    output:
        result=directory(f"{GENOMICS_MLST_ROOT}/{{assembly_accession}}"),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60,
    conda:
        "../envs/genomics-mlst.yaml"
    shell:
        r"""
        mkdir -p {output.result:q}
        set +e
        mlst --scheme ecoli_achtman_4 --full {input.fasta:q} \
          > {output.result:q}/mlst.tsv 2> {output.result:q}/command.stderr
        rc=$?
        set -e
        printf '%s\n' "$rc" > {output.result:q}/exit_code.txt
        """


rule genomics_amrfinder_assembly:
    input:
        fasta=_genomics_staged_fasta,
        database=GENOMICS_AMR_DB,
        database_metadata=GENOMICS_AMR_DB_METADATA,
    output:
        result=directory(f"{GENOMICS_AMR_ROOT}/{{assembly_accession}}"),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
    conda:
        "../envs/genomics-amrfinder.yaml"
    shell:
        r"""
        mkdir -p {output.result:q}
        set +e
        amrfinder --nucleotide {input.fasta:q} \
          --organism Escherichia \
          --plus \
          --name {wildcards.assembly_accession:q} \
          --threads {threads} \
          --database {input.database:q}/latest \
          --output {output.result:q}/amrfinder.tsv \
          > {output.result:q}/command.stdout 2> {output.result:q}/command.stderr
        rc=$?
        set -e
        printf '%s\n' "$rc" > {output.result:q}/exit_code.txt
        """


rule genomics_finalize_manifest:
    input:
        stage_table=lambda wildcards: checkpoints.genomics_stage_ncbi_package.get().output.table,
        quast=lambda wildcards: _genomics_analysis_dirs(wildcards, GENOMICS_QUAST_ROOT),
        mlst=lambda wildcards: _genomics_analysis_dirs(wildcards, GENOMICS_MLST_ROOT),
        amrfinder=lambda wildcards: _genomics_analysis_dirs(wildcards, GENOMICS_AMR_ROOT),
        versions=GENOMICS_TOOL_VERSIONS,
        mlst_info=GENOMICS_MLST_INFO,
        amrfinder_database_metadata=GENOMICS_AMR_DB_METADATA,
        amrfinder_database_checksums=GENOMICS_AMR_DB_CHECKSUMS,
    output:
        table=GENOMICS_FINAL_TABLE,
        exclusions=GENOMICS_FINAL_EXCLUSIONS,
        manifest=GENOMICS_FINAL_MANIFEST,
    params:
        min_total_length=GENOMICS_QC.get("min_total_length", 4000000),
        max_total_length=GENOMICS_QC.get("max_total_length", 6500000),
        max_contigs=GENOMICS_QC.get("max_contigs", 500),
        min_n50=GENOMICS_QC.get("min_n50", 20000),
        max_ns_per_100kb=GENOMICS_QC.get("max_ns_per_100kb", 1000.0),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    conda:
        "../envs/controller.yaml"
    shell:
        "python scripts/build_genome_manifest.py finalize "
        "--stage-table {input.stage_table:q} --quast-root {GENOMICS_QUAST_ROOT:q} "
        "--mlst-root {GENOMICS_MLST_ROOT:q} --amrfinder-root {GENOMICS_AMR_ROOT:q} "
        "--tool-versions {input.versions:q} --mlst-scheme-info {input.mlst_info:q} "
        "--amrfinder-database-metadata {input.amrfinder_database_metadata:q} "
        "--amrfinder-database-checksums {input.amrfinder_database_checksums:q} "
        "--output {output.table:q} --exclusions {output.exclusions:q} "
        "--manifest {output.manifest:q} "
        "--min-total-length {params.min_total_length} "
        "--max-total-length {params.max_total_length} "
        "--max-contigs {params.max_contigs} --min-n50 {params.min_n50} "
        "--max-ns-per-100kb {params.max_ns_per_100kb}"
