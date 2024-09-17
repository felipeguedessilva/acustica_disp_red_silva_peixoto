export DEVITO_ARCH=gcc
export DEVITO_PLATFORM=intel64
export DEVITO_LOGGING=DEBUG
export DEVITO_LANGUAGE=openmp
export DEVITO_AUTOTUNING=aggressive
export OMP_NUM_THREADS=24

nohup python test_exec.py &
#nohup python test_exec_multiprocessing.py &
#nohup python test_exec_threading.py &
