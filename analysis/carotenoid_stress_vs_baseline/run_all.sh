#!/usr/bin/env bash
# Full rerun of the a* analysis on the curated database (D-53, D-56, D-57): Cr, Cu, Pb only; groups from strain_curation; no population analyses.
# Run from the repo root as a SLURM job:
#   sbatch -p batch -c 8 --mem=48G -t 08:00:00 --wrap="cd $PWD && bash analysis/carotenoid_stress_vs_baseline/run_all.sh"
set -euo pipefail
A=analysis/carotenoid_stress_vs_baseline; S=$A/scripts; R=$A/results; P="pixi run"
step() { echo "=== $(date +%T) $*"; }
# 0. well tables and strain groups
step wells; $P python $S/prepare_wells.py --min-area 2000; $P python $S/prepare_wells.py --min-area 0
for a in 1000 3000 5000; do $P python $S/prepare_wells.py --min-area $a; done
$P python $S/prepare_wells.py --all-objects
step s0; $P python $S/strat/s0_inputs.py
# 1. sensitivity chain on the unfiltered wells (copied into report/sensitivity_unfiltered)
step unfiltered; cp $R/wells_nofilter.csv $R/wells.csv
$P Rscript $S/fit_mixed_models.R; $P python $S/plot_figures.py
$P Rscript $S/strat/s1_species_models.R; $P Rscript $S/strat/s346_models.R; $P python $S/strat/s5_split_half.py; $P python $S/strat/s_plots.py
U=$A/report/sensitivity_unfiltered; rm -rf $U; mkdir -p $U/tables $U/figures
cp $A/report/tables/* $U/tables/; cp $R/mixed_model_summary.csv $R/dose_factor_effects_M2.csv $U/tables/
for f in fig1_astar_vs_size_by_dose fig2_baseline_vs_stressed fig4_repeatability; do cp $R/figures/$f.png $U/figures/; done
for f in s1_species_baseline_forest s3_regimes s5_split_half_summary; do cp $A/report/figures/$f.png $U/figures/; done
mkdir -p $R/unfiltered; cp $R/mixed_model_summary.csv $R/dose_factor_effects_M2.csv $R/unfiltered/; cp $R/wells_nofilter.csv $R/unfiltered/wells.csv
# 2. minimum-area sensitivity
step minarea
for a in 1000 2000 3000 5000; do mkdir -p $R/minarea_$a; $P Rscript $S/fit_mixed_models.R $R/wells_min$a.csv $R/minarea_$a; done
# 3. main chain
step main; cp $R/wells_min2000.csv $R/wells.csv
$P Rscript $S/fit_mixed_models.R; $P python $S/plot_figures.py; $P python $S/traits_species.py
$P Rscript $S/strat/s1_species_models.R
$P python $S/strat/s2_phylo.py; $P python $S/strat/s2_sensitivity.py
$P Rscript $S/strat/s346_models.R; $P python $S/strat/s5_split_half.py; $P python $S/strat/s8_minarea_summary.py
$P Rscript $S/strat/s9_species_strata.R; $P python $S/strat/s9_species_plots.py
$P python $S/strat/s10_build_traits.py; $P Rscript $S/strat/s10_trait_models.R; $P python $S/strat/s10_b_morphology_plots.py
$P python $S/strat/s12_area_floor.py; $P python $S/strat/s_plots.py
step done
