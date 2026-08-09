# AD-GMRS

**AD-GMRS** stands for **Genotype-based Molecular Risk Score for Atopic Dermatitis**.

AD-GMRS is a lightweight command-line toolkit for calculating gene and protein feature scores from genotype data, standardizing these scores with a precomputed reference distribution, and outputting an atopic dermatitis risk score with a released lasso model.

## Authors

Hao Wu, Jia-Hao Wang, Wei Shi, Xin Ke

## Files

`AD-GMRS` is the main software entry. `AD-GMRS.script/` stores the internal workflow files used by the software, including `AD-GMRS.R`, `AD-GMRS.profile.txt`, and `AD-GMRS.path`. `README.md` is the software manual.

## Dependencies

The table below lists the runtime environment for the released prediction workflow.

| Program or Software | Version |
|---|---|
| Bash | `5.2.37` |
| PLINK 1.9 | `1.90 beta 7.11` |
| Python 3 | `3.11.15` |
| PRScs | `v1.1.0` |
| R | `4.5.3` |
| data.table | `1.18.4` |
| glmnet | `5.0` |

## Public reference data

Before running the software for the first time, edit `AD-GMRS.script/AD-GMRS.path` and replace the data paths in that file, including the 1000 Genomes path used by `AD_GMRS_BIM_TEMPLATE`.

| Public item | Path variable in `AD-GMRS.path` | Source |
|---|---|---|
| `PRScs.py` from the PRScs package | `AD_GMRS_PRSCS_PY` | [PRScs GitHub](https://github.com/getian107/PRScs) |
| PRScs LD reference directory such as `ldblk_1kg_EUR` | `AD_GMRS_REF_DIR` | [PRScs GitHub](https://github.com/getian107/PRScs) |
| 1000 Genomes Phase 3 PLINK BIM prefix such as `eur_chr{chr}` | `AD_GMRS_BIM_TEMPLATE` | [IGSR / 1000 Genomes](https://www.internationalgenome.org/data/) |

## Usage

Users need to prepare genotype data in PLINK binary format (`.bed/.bim/.fam`). If `AD-GMRS` is not in your `PATH`, run it in the project root as `./AD-GMRS`.

Predict risk for a new genotype dataset:

```bash
AD-GMRS --geno /path/to/new_sample_chr{chr}
```

This command calculates gene or protein scores from the genotype input, standardizes them with `AD-GMRS.profile.txt`, and writes risk results to `AD-GMRS.result/`. If `weight_manifest.tsv` is stored inside `AD_GMRS_WEIGHTS_DIR`, no extra manifest setting is required.

Optional overrides can still be passed when needed, for example `--outdir`, `--keep`, `--weights`, `--model`, `--manifest`, and `--profile-reference`.