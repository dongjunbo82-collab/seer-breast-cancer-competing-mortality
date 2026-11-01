# SEER breast cancer competing mortality analysis

Statistical source code for five-year sex differences in breast cancer and competing-cause mortality, with detailed cause decomposition and continuous-age analyses.

This repository contains source code, an input-column schema, and execution instructions only. It contains no manuscript, patient records, completed datasets, fitted models, numerical results, figures, credentials, or author contact details. Obtain SEER records independently through the authorized access process: https://seer.cancer.gov/data/access.html. A data-use agreement is required; do not upload records or outputs to this repository.

## Independent recent cardiovascular age checks

Run `Rscript work/jjco_recent_cv/analyse.R` after authorized input reconstruction. It fits three competing hazards in 2015–2018 alone (breast cancer, cardiovascular and all remaining deaths), with sex-specific log-linear age, complete-case and unknown-category specifications; it also estimates sex/age-group Aalen–Johansen risks. Outputs include private fitted objects that must remain local. This checks direction and magnitude under a simpler age form; it does not validate an exact historical spline or a screening threshold. Table S22 reports observed case mix; Table S23 reports recent male non-age standardization. `clinical_other.R` then uses the original 40 saved fits to supply noncardiovascular other-cause risks at ages 50/65/75 for Table S21, retaining prediction covariance. These scripts require locally generated authorized inputs and original fits, not files supplied by the repository.

## Current JJCO additions

`work/jjco_value_revision/background_matched.R` supplies national background context. `work/jjco_robustness_final/analyse.R` runs independent temporal and boundary checks; `shared_reference.R` predicts at a common recent reference. `work/jjco_clinical_final/nonCOVID.R` sums non-COVID other-cause cumulative incidence while retaining COVID as a competing event and preserving covariance. Ages 65/70 are in the original grid. These scripts require authorized inputs and saved principal fits. Set SEER_ANALYSIS_ROOT locally and run from the repository root. No results or patient records are distributed.

## Current analysis

The principal workflow is `work/revision_20261019/`: seven competing hazards, four imputation chains, 40 completed cohorts, sex-specific age terms where specified, and standardization to a pooled male case mix. `work/cbc_robustness20/` contains alternative-knot sensitivity analyses and scalar chain diagnostics. Directory dates identify source versions, not distinct study populations.

The analysis is exploratory. Standardized comparisons describe a common measured case mix; reference restrictions do not eliminate model extrapolation. Ten retained completions per chain limit ESS and split-Rhat assessment. Alternative-knot models reuse the primary completions and do not re-impute under each alternative formula.

## Requirements and inputs

Developed with R 4.6.1 and Python 3. R packages: data.table, splines, nnet, MASS, survival, ggplot2, patchwork, svglite, ragg, posterior (diagnostics developed with 1.7.0). Python preparation uses the standard library; the optional life-table parser also requires xlrd. `Rscript setup.R` checks R dependencies and prepares local output directories; it does not install software or packages.

Export formatted SEER case listings with SEER*Stat 9.0.43.0. `export_schema.json` specifies the two input schemas and exact column order. These are column definitions, not patient data. Set `SEER_BASE_EXPORT` and `SEER_COD_EXPORT` to your local 30-column and 34-column CSV files, with corresponding records in identical order. Existing database and selection versions must match the study; the code does not download SEER data or authenticate to SEER.

Run commands from the repository root. The default data root is `.`; `SEER_ANALYSIS_ROOT` may point to a different reconstruction directory. The primary result directory defaults to `outputs/revision_20261019`; `SEER_REVISION_OUT` can override it and must be consistent across primary scripts. Preparation and contextual scripts use repository-relative `work/` paths. New spline diagnostics write to `outputs/cbc_robustness20`.

## Execution order

```sh
Rscript setup.R
python work/audit_cod_pair.py
python work/build_extended.py
python work/build_cod_extended.py
Rscript work/cbc_robustness20/prepare_current_cohort.R
```

Run each imputation chain once with chain argument 1, 2, 3, or 4:

```sh
Rscript work/revision_20261019/impute.R 1 40
Rscript work/revision_20261019/impute.R 2 40
Rscript work/revision_20261019/impute.R 3 40
Rscript work/revision_20261019/impute.R 4 40
Rscript work/revision_20261019/analyze.R baseline
```

For each chain argument 1–4, run both commands below after its imputations exist. The chain scripts poll for incomplete jobs and can run for substantial time; launch each chain only once in the same output directory.

```sh
Rscript work/revision_20261019/analyze.R mi 1
Rscript work/revision_20261019/compare_models.R 1
```

After all four chains and comparisons finish:

```sh
Rscript work/revision_20261019/pool.R
Rscript work/revision_20261019/support_grid.R
Rscript work/revision_20261019/validate_likelihood.R
Rscript work/revision_20261019/audit.R
Rscript work/revision_20261019/figures.R
```

The comparison script reads only function definitions from `work/models.R`; do not run that historical engine before the current seven-cause analysis.

For alternative knots, run `Rscript work/cbc_robustness20/robustness20.R CHAIN 10` separately for CHAIN 1, 2, 3, and 4, then:

```sh
Rscript work/cbc_robustness20/pool_robustness20.R
Rscript work/cbc_robustness20/chain_diagnostics20.R
```

Scalar diagnostics require locally generated per-imputation estimates. No aggregate inputs are bundled in this code-only release.

## Analysis-to-output map

| Source | Local output/purpose |
| --- | --- |
| `analyze.R`, `pool.R` | `revised_MI.csv`: standardized risks, cause differences, continuous age, subgroup, restricted references; principal Tables 2–3, supplemental S4–S8/S12 |
| `analyze.R baseline` | `observed_AJ.csv`: observed risks; `revised_sensitivity.csv`: complete-case and unknown-category checks |
| `compare_models.R`, `pool.R` | `model_comparison.csv`: shared-age/sex-specific-age three-cause comparisons |
| `impute.R`, `pool.R` | category traces, FMI and Monte Carlo diagnostics |
| `support_grid.R` | `support_grid.csv`, `support_summary.csv`: prediction support |
| `figures.R` | current principal and diagnostic figure generation |
| `robustness20.R`, `pool_robustness20.R` | two alternative age-knot specifications, complete age-grid estimates, fit checks; supplemental S16 |
| `chain_diagnostics20.R` | lag-one correlation, basic unsplit ESS, split Rhat and correlation-aware MCSE; supplemental S17 |
| `validate_likelihood.R`, `audit.R`, `baseline_audit.R`, `benchmark_likelihood.R` | numerical checks and optional fit audits; require locally generated data/fits |

## Contextual analyses and supporting utilities

These scripts use separate three-cause specifications and must not replace the current principal workflow. They are retained only for context and necessary function dependencies.

- `models.R`: three-cause functions and contextual base workflow. `LOAD_ONLY=1` loads its data/functions without executing the final analysis block; the current comparison script instead parses function definitions only.
- `value_add_analysis.R`: common three-year pre-pandemic window and age-stratum checks, corresponding to supplemental S13 context.
- `stage_trends.R`: stage-distribution models, supplemental S14 context.
- `sensitivity.R`: contextual three-cause checks, including the under-85 comparison needed by `background.R`.
- `get_lifetables.py`, `life_parse.py`, `background.R`: download public CDC life tables, parse annual death probabilities, and calculate descriptive population benchmarks, supplemental S15 context. Run the fetch/parser before `background.R`; run `sensitivity.R` first to generate its under-85 observed comparison. These are descriptive benchmarks, not formal excess-mortality models.
- `audit_seer.py`, `summarize_followup.py`: local input/cohort and follow-up summaries. Any generated audit JSON or CSV may contain sensitive source values and must remain local.

For contextual workflows run `Rscript work/models.R` once to generate contextual base fits, then the needed contextual script. This creates local records, models and outputs covered by `.gitignore`; none are repository content. Existing output files can be reused by some scripts, so use a clean reconstruction directory for a full rerun. Output filenames in source retain historical numbering; the table map above identifies their current use.

## Publication boundary and validation

Only explicitly allowlisted source files, this README, `.gitignore`, and `export_schema.json` are intended for publication. The input schema includes record-key column names solely to define matching; no identifier values are provided. Do not add local inputs, fitted objects, traces, aggregates, plots, logs, personal paths, or contact details. `.gitignore` is a convenience and does not replace review of staged files.

The code-only release is checked for Python/R syntax, local source dependencies, forbidden file types and obvious credentials/contact paths. The underlying analyses were executed in the private research workspace. A full rerun cannot be performed from this repository alone without authorized SEER input files.


## Extended-chain stability check

After the current primary analysis is complete, continue each of the four chains from its saved sweep-40 completed cohort:

```sh
Rscript work/chain_stability/extend_chains.R 1 120
Rscript work/chain_stability/extend_chains.R 2 120
Rscript work/chain_stability/extend_chains.R 3 120
Rscript work/chain_stability/extend_chains.R 4 120
Rscript work/chain_stability/pool_extension.R 120
Rscript work/chain_stability/plot_extension.R
```

The original random-number states were not available. The saved completed covariates at sweep 40 provide the starting chain states; continuation uses independent seeds `20261020 + 100000 * chain`. This continues the same transition kernel with fresh random streams rather than reproducing the original streams. Each continuation saves its completed covariate state and RNG state locally after each sweep; rerunning the same command resumes that checkpoint. Checkpoints contain patient-level covariates and must never be shared.

The first 20 continuation sweeps are excluded. Key overall and age-50/75 effects are evaluated at sweeps 62,64,...,120, yielding 30 retained completions per chain (120 total). They are pooled separately and compared with the principal 40-completion analysis; they do not replace the main analysis or regenerate the alternative-knot models. The continuation uses the existing vectorized exact category likelihood. `validate_kernel.R` compares it with explicit category-by-category likelihood evaluation using locally available primary fits.

The output directory defaults to `outputs/chain_stability` and can be changed with `SEER_EXTENSION_OUT`. It contains only local runtime products, including private checkpoints. The code repository supplies no diagnostic numbers or supplementary CSV files: numerical Supplementary Results accompany the article separately. Low MCSE is not a convergence criterion; examine split Rhat, effective sample size, chain traces, and changes in key estimates and standard errors together. The scalar summaries cover overall sex differences, ages 50 and 75, and their contrasts; this check does not refit every subgroup or restricted-reference analysis.

The extended-chain comparison corresponds to supplementary Table S18 and Figure S9. The full 36-estimand comparison and category/trace CSVs are submitted as independent Supplementary Results attachments with the article, not as files in this repository.
