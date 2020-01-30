source ~/.bashrc
conda activate hotspot
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
#snakemake spatialDE/timing_pairs/log_500_1000.json
snakemake timingsPairs -j16
