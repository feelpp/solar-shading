# Prerequisites

The submodules may have to be initialized and updated:

```bash
git submodule update --init --recursive
```

More, the Intel oneAPI MKL library has to be installed and, either the environment variable `MKLROOT` has to be set, or the include and library paths have to be added to the `C_INCLUDE_PATH` and `LD_LIBRARY_PATH` environment variables, or directly in the 'runbenchmark.sh' script by changing the commented part of the line:

```bash
-I/home/u4/csmi/2022/devora/intel/oneapi/mkl/2023.1.0/include # -L/home/u4/csmi/2022/devora/intel/oneapi/mkl/2023.1.0/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
```

# Warnings 

Intel's oneAPI MKL library isn't installed and hence isn't tested in the script.

# Compilation and execution

When the prerequisites are met, the benchmarks can be compiled by running the following command:

```bash
cd benchmark
cd distributions
chmod +x runbenchmark.sh
./runbenchmark.sh
```
And then either type `uniform_int` or `uniform_real` when asked to choose the distribution to benchmark.

After running this, the benchmarks can be found in the benchmark/distributions directory next to the scripts. The benchmarks are composed of two files, `uniform_int.csv` (or `uniform_real.csv`) containing the results of the benchmarks and `uniform_int.png` (or `uniform_real.png`) containing the plots of the benchmarks in ascending order of execution time.


After each execution, the results are saved in these two files, but if you want to rerun the benchmark, you will have to respond 'Y' to the question 'Do you want to erase the previous results ?' in the script.

The `-march=native` flag is used to compile the benchmarks with the best possible optimization for the machine on which the benchmarks are run. If you want to compile the benchmarks for an execution on another machine, you can change this flag in the script by using `-mavx2` or `-msse4.2` for example. 