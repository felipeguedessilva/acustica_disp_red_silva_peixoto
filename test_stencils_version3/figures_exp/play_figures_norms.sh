# Devito Config
export DEVITO_ARCH=gcc
export DEVITO_PLATFORM=intel64
export DEVITO_LANGUAGE=openmp
export OMP_NUM_THREADS=40

nohup python test_figures_norms.py > compile_results_figures_norms.out 2>&1 &