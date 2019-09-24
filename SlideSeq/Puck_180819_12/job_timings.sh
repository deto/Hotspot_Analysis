source ~/.bashrc
conda activate hotspot
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
snakemake timings
