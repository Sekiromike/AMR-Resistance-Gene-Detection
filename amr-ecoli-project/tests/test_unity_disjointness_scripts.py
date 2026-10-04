from __future__ import annotations

import json
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]


class UnityDisjointnessScriptTests(unittest.TestCase):
    def test_raw_prefetch_retry_is_bounded_and_attempt_scoped(self) -> None:
        script = (ROOT / "hpc" / "unity" / "disjointness-raw-array.sbatch").read_text(
            encoding="utf-8"
        )
        self.assertIn("prefetch_attempt_limit=3", script)
        self.assertIn('attempt_dir="${work}/sra/attempt-${prefetch_attempt}"', script)
        self.assertIn('retry_delay=$((prefetch_attempt * 30))', script)
        self.assertIn('exit "${prefetch_exit}"', script)
        self.assertIn('prefetch_attempts=%s', script)

    def test_disjointness_scripts_use_only_parameterized_project_storage(self) -> None:
        for name in (
            "disjointness-raw-array.sbatch",
            "disjointness-audit.sbatch",
            "disjointness-compare.sbatch",
            "development-clusters.sbatch",
            "development-clusters-resume.sbatch",
            "external-clusters.sbatch",
            "development-duplicate-pairs.sbatch",
            "genomics-setup.sbatch",
            "genome-analysis-array.sbatch",
            "cohort-build.sbatch",
            "mlst-array.sbatch",
            "development-baseline.sbatch",
            "development-models.sbatch",
            "model-comparison.sbatch",
            "models-env.sbatch",
            "model-round.sbatch",
            "loci-reference.sbatch",
            "loci-array.sbatch",
            "loci-collect.sbatch",
            "embeddings-units.sbatch",
            "pretraining-overlap.sbatch",
            "source-robustness.sbatch",
            "external-predict.sbatch",
            "external-evaluate.sbatch",
            "posthoc-family-predict.sbatch",
            "posthoc-family-score.sbatch",
            "leave-variant-out.sbatch",
            "leave-variant-readout.sbatch",
            "embeddings-canary.sbatch",
            "embeddings-env.sbatch",
        ):
            script = (ROOT / "hpc" / "unity" / name).read_text(encoding="utf-8")
            self.assertIn("AMR_WORK_ROOT", script, name)
            self.assertNotIn("/tmp/", script, name)
            self.assertNotIn("/home/", script, name)
            self.assertNotIn("/scratch/", script, name)
            self.assertNotIn("/project/", script, name)

    def test_development_clusters_job_fails_closed_on_an_incomplete_source(self) -> None:
        """Clustering must never reuse a partial sketch from an unfinished run."""
        script = (ROOT / "hpc" / "unity" / "development-clusters.sbatch").read_text(
            encoding="utf-8"
        )
        self.assertIn('AMR_SOURCE_COMPARISON:?', script)
        self.assertIn('if [[ ! -f "${SOURCE_ROOT}/COMPLETE" ]]; then', script)
        self.assertIn("sha256sum --check", script)
        self.assertIn("Refusing to overwrite cluster result", script)
        # Reuses the frozen sketch rather than rebuilding it.
        self.assertIn("sketch_reused_from", script)
        self.assertNotIn("skani sketch", script)


class UnityPolicyComplianceTests(unittest.TestCase):
    """Guards for the login-node and storage rules in docs/UNITY_HPC.md."""

    def test_every_conda_env_job_redirects_the_package_cache_off_home(self) -> None:
        for name in ("disjointness-env.sbatch", "embeddings-env.sbatch", "models-env.sbatch"):
            script = (ROOT / "hpc" / "unity" / name).read_text(encoding="utf-8")
            self.assertIn("CONDA_PKGS_DIRS", script, name)
            self.assertIn("AMR_CONDA_PKGS", script, name)
            # A mention inside a comment is fine; an executable line that
            # actually points the cache at home is not.
            executable = [
                line for line in script.splitlines() if not line.lstrip().startswith("#")
            ]
            self.assertFalse(
                [line for line in executable if "~/.conda" in line or "$HOME/.conda" in line],
                name,
            )

    def test_gpu_environment_build_requests_a_gpu(self) -> None:
        """conda's __cuda virtual package only resolves where a GPU is present."""
        script = (ROOT / "hpc" / "unity" / "embeddings-env.sbatch").read_text(
            encoding="utf-8"
        )
        self.assertIn("#SBATCH --gpus=1", script)
        self.assertIn("conda env create", script)

    def test_preempt_partition_jobs_stay_under_the_two_hour_kill_window(self) -> None:
        import re

        for path in sorted((ROOT / "hpc" / "unity").glob("*.sbatch")):
            script = path.read_text(encoding="utf-8")
            if "preempt" not in script:
                continue
            match = re.search(r"#SBATCH --time=(\d+):(\d+):(\d+)", script)
            self.assertIsNotNone(match, path.name)
            hours, minutes, seconds = (int(value) for value in match.groups())
            total = hours * 3600 + minutes * 60 + seconds
            self.assertLessEqual(
                total, 2 * 3600, f"{path.name} exceeds the preempt kill window"
            )

    def test_declared_gpu_constraints_are_valid_unity_names(self) -> None:
        import re

        # Taken from the live cluster (`sinfo -o %f`), not the published list,
        # which omits h100/h200/vram143/sm_87/sm_120.
        valid = {
            f"vram{size}" for size in (8, 11, 12, 16, 23, 32, 40, 48, 80, 102, 143)
        } | {
            f"sm_{cap}" for cap in (52, 61, 70, 75, 80, 86, 87, 89, 90, 120)
        } | {
            "titanx", "m40", "1080ti", "v100", "2080", "2080ti", "rtx8000",
            "a100", "a100-40g", "a100-80g", "a4000", "a16", "a40", "gh200",
            "l40s", "l4", "h100", "h200", "h200_nvl",
        }
        for path in sorted((ROOT / "hpc" / "unity").glob("*.sbatch")):
            script = path.read_text(encoding="utf-8")
            for match in re.findall(r"#SBATCH --constraint=(\S+)", script):
                for token in re.split(r"[&|]", match):
                    self.assertIn(token, valid, f"{path.name}: {token}")


class ArchitectureCompatibilityTests(unittest.TestCase):
    def test_no_job_requests_the_arm_gh200(self) -> None:
        """gh200 is Grace Hopper (aarch64); the x86_64 conda stack cannot run there."""
        for path in sorted((ROOT / "hpc" / "unity").glob("*.sbatch")):
            script = path.read_text(encoding="utf-8")
            executable = [
                line for line in script.splitlines() if not line.lstrip().startswith("# ")
            ]
            self.assertFalse(
                [line for line in executable if "gh200" in line or "arm-gpu" in line],
                f"{path.name} targets the ARM GPU partition",
            )


class ResumeGateTests(unittest.TestCase):
    def test_resume_never_reruns_search_and_verifies_before_writing(self) -> None:
        script = (ROOT / "hpc" / "unity" / "development-clusters-resume.sbatch").read_text(
            encoding="utf-8"
        )
        executable = "\n".join(
            line for line in script.splitlines() if not line.lstrip().startswith("#")
        )
        self.assertNotIn("/bin/skani", executable, "resume must not invoke skani at all")
        gate = script.index("verification=PASSED")
        resume = script.index("--output \"${RESUME_ROOT}/comparisons.tsv\"")
        self.assertLess(gate, resume, "development output written before verification gate")
        self.assertIn("Refusing to overwrite existing", script)
        self.assertIn("resume_distinct_queries", script)


class EmbeddingCanaryTests(unittest.TestCase):
    def test_canary_is_bounded_pins_a_revision_and_keeps_cache_in_allocation(self) -> None:
        script = (ROOT / "hpc" / "unity" / "embeddings-canary.sbatch").read_text(
            encoding="utf-8"
        )
        # One GPU, short wall time, preemptible partition per docs/UNITY_HPC.md.
        self.assertIn("#SBATCH --gpus=1", script)
        self.assertIn("--partition=gpu-preempt", script)
        # A moving model tag is not reproducible provenance.
        self.assertIn("AMR_MODEL_REVISION:?", script)
        self.assertIn("AMR_MODEL_ID:?", script)
        # The HuggingFace cache must not land in the 100 GB home directory.
        self.assertIn('HF_HOME="${WORK_ROOT}/huggingface"', script)
        self.assertIn("Refusing to overwrite canary result", script)


class GenomeAnalysisScriptTests(unittest.TestCase):
    """Regressions from genomics-setup job 64925569."""

    def test_interpreted_tools_run_with_their_environment_first_on_path(self) -> None:
        # mlst is '#!/usr/bin/env perl'; without its env on PATH it used the
        # system perl and failed to find List::MoreUtils.
        script = (ROOT / "hpc" / "unity" / "genome-analysis-array.sbatch").read_text(encoding="utf-8")
        self.assertIn('PATH="${ENVS}/amr-genomics-mlst/bin:${PATH}" "${MLST}"', script)
        self.assertIn('PATH="${ENVS}/amr-genomics-qc/bin:${PATH}" "${QUAST}"', script)

    def test_quast_environment_pins_a_python_that_ships_distutils(self) -> None:
        spec = (ROOT / "workflow" / "envs" / "genomics-qc.yaml").read_text(encoding="utf-8")
        self.assertIn("- python=3.11", spec)

    def test_setup_checks_reference_species_by_taxid_and_tool_versions(self) -> None:
        script = (ROOT / "hpc" / "unity" / "genomics-setup.sbatch").read_text(encoding="utf-8")
        self.assertIn('species.get("id") != 562', script)
        self.assertIn("summary taxonomy taxon", script)
        self.assertNotIn('organism != "Escherichia coli"', script)
        self.assertIn("QUAST unusable", script)
        self.assertIn("mlst unusable", script)


class CohortBuildScriptTests(unittest.TestCase):
    def test_cohort_build_verifies_sources_and_fails_closed(self) -> None:
        script = (ROOT / "hpc" / "unity" / "cohort-build.sbatch").read_text(encoding="utf-8")
        self.assertIn("sha256sum --check --quiet artifact.sha256", script)
        self.assertIn("Genome analysis incomplete", script)
        self.assertIn("Refusing to overwrite cohort result", script)
        # Adequacy runs before metadata, and a non-keep outcome stops the job.
        self.assertLess(script.index("assess_lineage_adequacy.py"), script.index("build_isolate_metadata.py"))
        self.assertIn('exit "${adequacy_exit}"', script)
        # Both MIC sources and both cluster files; primary endpoint profile.
        self.assertEqual(script.count("--interpreted-ast"), 2)
        self.assertEqual(script.count("--clusters "), 2)
        self.assertEqual(script.count("--endpoint-profile mic_regression"), 2)
        # COMPLETE only after validation passes.
        self.assertLess(script.index("validation_exit=$?"), script.index('touch "${RESULT_ROOT}/COMPLETE"'))


class MlstSchemeTests(unittest.TestCase):
    """Regression: 'ecoli' is the Pasteur scheme in mlst >= 2.23."""

    def test_no_script_types_with_the_bare_ecoli_scheme(self) -> None:
        for path in list((ROOT / "hpc" / "unity").glob("*.sbatch")) + [ROOT / "workflow" / "rules" / "genomics.smk"]:
            text = path.read_text(encoding="utf-8")
            self.assertNotRegex(text, r"--scheme ecoli(\s|$)", path.name)

    def test_arrays_read_the_scheme_from_config_and_cohort_build_uses_the_retyping(self) -> None:
        for name in ("mlst-array.sbatch", "genome-analysis-array.sbatch"):
            text = (ROOT / "hpc" / "unity" / name).read_text(encoding="utf-8")
            self.assertIn('["genomics"]["mlst"]', text, name)
        build = (ROOT / "hpc" / "unity" / "cohort-build.sbatch").read_text(encoding="utf-8")
        self.assertIn('--mlst-root "${MLST_ROOT}"', build)
        self.assertIn("MLST re-typing incomplete", build)
        config = json.loads((ROOT / "config" / "study.json").read_text(encoding="utf-8"))
        self.assertEqual(config["genomics"]["mlst"]["scheme"], "ecoli_achtman_4")
        self.assertEqual(sorted(config["genomics"]["mlst"]["loci"]),
                         sorted(["adk", "fumC", "gyrB", "icd", "mdh", "purA", "recA"]))


class DevelopmentBaselineScriptTests(unittest.TestCase):
    def test_baseline_job_requires_a_complete_verified_cohort_and_never_locks_external(self) -> None:
        script = (ROOT / "hpc" / "unity" / "development-baseline.sbatch").read_text(encoding="utf-8")
        self.assertIn('[[ -f "${COHORT_ROOT}/COMPLETE" ]]', script)
        self.assertIn("sha256sum --check --quiet artifact.sha256", script)
        self.assertNotIn("--locked-external", script)
        self.assertIn("for population in human_clinical all_sources", script)


class DevelopmentModelsScriptTests(unittest.TestCase):
    def test_models_job_is_development_only_and_verifies_inputs(self) -> None:
        script = (ROOT / "hpc" / "unity" / "development-models.sbatch").read_text(encoding="utf-8")
        self.assertNotIn("--locked-external", script)
        self.assertIn("comparisons.tsv changed since clustering", script)
        self.assertIn("sha256sum --check --quiet artifact.sha256", script)
        self.assertIn("for model in constant neighbor amrfinder", script)
        self.assertLess(script.index("build_split_manifest.py"), script.index("run_censored_regression.py"))


class LmodTests(unittest.TestCase):
    def test_every_job_that_calls_module_loads_lmod_first(self) -> None:
        """Regression: job 64980720 failed with 'module: command not found'."""
        for path in sorted((ROOT / "hpc" / "unity").glob("*.sbatch")):
            text = path.read_text(encoding="utf-8")
            if "module load" not in text and "module purge" not in text:
                continue
            first_use = min(i for i in (text.find("module load"), text.find("module purge")) if i >= 0)
            self.assertIn("z00-lmod-profile.sh", text[:first_use], path.name)


class ExternalSealScriptTests(unittest.TestCase):
    def test_prediction_job_never_scores_and_evaluation_requires_the_lock(self) -> None:
        predict = (ROOT / "hpc" / "unity" / "external-predict.sbatch").read_text(encoding="utf-8")
        evaluate = (ROOT / "hpc" / "unity" / "external-evaluate.sbatch").read_text(encoding="utf-8")
        self.assertNotIn("evaluate_mic_predictions", predict)
        self.assertNotIn("attach_external_references", predict)
        self.assertIn("config/external_prediction_lock.tsv", evaluate)
        self.assertIn("--locked-external", evaluate)
        self.assertLess(evaluate.index("attach_external_references.py"), evaluate.index("evaluate_mic_predictions.py"))


if __name__ == "__main__":
    unittest.main()
