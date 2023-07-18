# Prerequisites

The submodules may have to be initialized and updated:

```bash
git submodule update --init --recursive
```

More, the Intel oneAPI MKL library has to be installed and, either the environment variable `MKLROOT` has to be set, or the include and library paths have to be added to the `C_INCLUDE_PATH` and `LD_LIBRARY_PATH` environment variables, or directly in the 'runbenchmark.sh' script by changing the commented part of the line:

```bash
-I/home/u4/csmi/2022/devora/intel/oneapi/mkl/2023.1.0/include # -L/home/u4/csmi/2022/devora/intel/oneapi/mkl/2023.1.0/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
```

# Benchmarking

## With the Bash script

The `-march=native` flag is used to compile the benchmarks with the best possible optimization for the machine on which the benchmarks are run. If you want to compile the benchmarks for an execution on another machine, you can change this flag in the script by using `-mavx2` or `-msse4.2` for example. 

To run the unique benchmarks for the int and real distributions, one can launch the `benchmark.sh` script:

```bash
cd benchmark/distributions
chmod +x benchmark.sh
./benchmark.sh
int # or real after asked to choose the distribution to benchmark
```
After running this, the benchmark's results can be found in the results/csv directory at the root of this repo. The benchmarks are composed of two files, `uniform_int.csv` (or `uniform_real.csv`) containing the results of the benchmarks and `uniform_int.png` (or `uniform_real.png`) containing the plots of the benchmarks in ascending order of execution time (available inthe results/images folder).

The `-march=native` flag is used to compile the benchmarks with the best possible optimization for the machine on which the benchmarks are run. If you want to compile the benchmarks for an execution on another machine, you can change this flag in the script by using `-mavx2` or `-msse4.2` for example. 
## With Reframe

After installing reframe on your machine and adding the `reframe` command to your path, you can run the benchmarks with the following command when located at the root of the repository:

```bash
cd benchmark/distributions
reframe -c benchmarking_int.py -r --performance-report
```

The results will be displayed in the terminal and individual performance logs are generated in the `benchmark/distributions/output/*/rfm_job.out` files.

One can also run the benchmarks for the real distributions with the following command:

```bash
cd benchmark/distributions
reframe -c benchmarking_real.py -r --performance-report
```

The generated files are located at the same place as for the integer distributions.

# Warnings 

Intel's oneAPI MKL library isn't installed and hence isn't tested in the script.
