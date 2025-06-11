#Limpa o cache de módulos
module purge
module load apptainer

# Devito Config
export DEVITO_ARCH=gcc
export DEVITO_PLATFORM=intel64
export DEVITO_LANGUAGE=openmp
export OMP_NUM_THREADS=40

#Execut
nohup /home/public/devito_cont_mintrop/devito_develop.sif python3 test_exec.py > compile_results.out 2>&1 &