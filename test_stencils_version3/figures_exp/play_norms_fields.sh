# Devito Config
export DEVITO_ARCH=gcc
export DEVITO_PLATFORM=intel64
export DEVITO_LANGUAGE=openmp
export OMP_NUM_THREADS=40

nohup python rcs_files_norms_fields.py > compile_results_norms_fields.out 2>&1 &