# Devito Config
export DEVITO_ARCH=gcc
export DEVITO_PLATFORM=intel64
export DEVITO_LANGUAGE=openmp
export OMP_NUM_THREADS=40

nohup python signal_analysis_cr.py > compile_results_signal_analysis_cr.out 2>&1 &