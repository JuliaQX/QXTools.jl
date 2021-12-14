# Running Distributed Simulations

## Generating input files

To use QXContexts we first generate a set of input files using the `generate_simulation_files`
function.
This function accepts any circuit type for which the `convert_to_tnc` function has been implemented.
This includes QXZoo and YaoBlocks circuits, see the [Circuits with QXZoo and YaoQX](circuits.md)
tutorial for more details on creating such circuits.
The following example will use a 49 qubit Random Quantum Circuit (RQC) of depth 16.

```
using QXTools

circ = create_rqc_circuit(7, 7, 16, 42, final_h=true)
output_args = output_params_dict(49, 100)
generate_simulation_files(circ, "rqc_7_7_16", time=10, output_args=output_args)
```

The above example will run the contraction finding algorithm for 10 seconds and write input files
`rqc_7_7_16.qx`, `rqc_7_7_16.yaml` and `rqc_7_7_16.jld2`.
By default two bonds will be sliced and the simulation files will specify that simulation should be run for 100 random bitstrings. Further details on the format and contents of these files is described in the [User Guide](users_guide.md) and further details of the slicing and output sampling options of the `generate_simulation_files` function are described in the [Miscellaneous Features](features.md) tutorial.

## Running the simulation

Once the input files have been generated, simulations can be run using the `bin/qxrun.jl` script which can be run directly from the command line (not the Julia REPL).
To run a simulation using the input files generated in the example above one would use

```
julia --project bin/qxrun.jl -d rqc_7_7_16.qx -o rqc_7_7_16_output.jld2
```

This will run the simulation and write the output to the `rqc_7_7_16_output.jld2` file.

To run simulations in parallel the MPI.jl package must first be installed and the `mpiexecjl` utility
installed.
This can be done with the following commands:

```
]add MPI
import MPI
MPI.install_mpiexecjl()
```

On compute clusters it is advised to use the system provided MPI installation. For more details on this see the official MPI.jl documentation [here](https://juliaparallel.github.io/MPI.jl/stable/configuration/#Using-a-system-provided-MPI).

Once MPI.jl has been installed and configured the simulation example described above can be run in parallel on two processes with

```
mpiexecjl --project -n 2 julia  bin/qxrun.jl -d rqc_7_7_16.qx -o rqc_7_7_16_output.jld2 -m
```

where the `-n 2` specifies the number of processes to use and the `-m` option enabled MPI.

## Measuring the time and speedup

When running the above example we would expect to see a reduction in the time taken when using more processes, however this is not what we observe in practice.
This is due to the fact that with Julia when running code for the first time the majority of the time is taken up by compilation.
In the above example since we are only calculating amplitudes for 10 bitstrings there are not enough computations after this compilation to observe a speedup.
In practice for circuits of this size one would be calculating amplitudes for millions of bitstrings.
It is also possible to reduce this startup time by compiling a custom system image as decribed in [this documentation](https://juliaqx.github.io/QXContexts.jl/dev/#Custom-system-image) for QXContexts.

To see more clearly the time taken for compilation vs computation one can run the examples above using the `-t` flag. When this flag is used the simulation will be performed twice and timing information collected from each of these runs.
(We also add `-b 1` which sets the number of BLAS threads to 1 to see clearly the effect of increasing the number of processes without changing the contention between cores that would result when threading is also used.)
Using a single process with the following command

```
mpiexecjl --project -n 1 julia bin/qxrun.jl -d rqc_7_7_16.qx -o rqc_7_7_16_output.jld2 -m  -t -b 1
```

results in the following timing output

```
 ──────────────────────────────────────────────────────────────────────────────
                                       Time                   Allocations
                               ──────────────────────   ───────────────────────
       Tot / % measured:          1127914s / 0.01%          15.5GiB / 90.5%

 Section               ncalls     time   %tot     avg     alloc   %tot      avg
 ──────────────────────────────────────────────────────────────────────────────
 Simulation                 1    52.6s  80.5%   52.6s   13.0GiB  93.0%  13.0GiB
 Init sampler               1    10.7s  16.3%   10.7s    759MiB  5.29%   759MiB
   Parse input files        1    7.80s  11.9%   7.80s    562MiB  3.91%   562MiB
   Create Context           1    1.34s  2.05%   1.34s    147MiB  1.03%   147MiB
   Create sampler           1   73.5ms  0.11%  73.5ms   3.88MiB  0.03%  3.88MiB
 Write results              1    2.12s  3.24%   2.12s    251MiB  1.75%   251MiB
 ──────────────────────────────────────────────────────────────────────────────
 ──────────────────────────────────────────────────────────────────────────────
                                       Time                   Allocations
                               ──────────────────────   ───────────────────────
       Tot / % measured:            16.1s / 100%            5.95GiB / 100%

 Section               ncalls     time   %tot     avg     alloc   %tot      avg
 ──────────────────────────────────────────────────────────────────────────────
 Simulation                 1    16.0s   100%   16.0s   5.94GiB  100%   5.94GiB
 Init sampler               1   60.2ms  0.37%  60.2ms   11.8MiB  0.19%  11.8MiB
   Parse input files        1   41.5ms  0.26%  41.5ms   4.80MiB  0.08%  4.80MiB
   Create Context           1   15.6ms  0.10%  15.6ms   6.78MiB  0.11%  6.78MiB
   Create sampler           1   69.3μs  0.00%  69.3μs   2.08KiB  0.00%  2.08KiB
 Write results              1   1.60ms  0.01%  1.60ms   72.1KiB  0.00%  72.1KiB
 ──────────────────────────────────────────────────────────────────────────────
```

 The first set of timings include the precompilation whereas the second set are for the computations alone.
 From this we see that the simulation part of the code took 16 seconds when using a single process.
 Running on two processes with

```
mpiexecjl --project -n 2 julia bin/qxrun.jl -d rqc_7_7_16.qx -o rqc_7_7_16_output.jld2 -m  -t -b 1
```

we get

```
 ──────────────────────────────────────────────────────────────────────────────
                                       Time                   Allocations
                               ──────────────────────   ───────────────────────
       Tot / % measured:          1129782s / 0.01%          12.5GiB / 88.2%

 Section               ncalls     time   %tot     avg     alloc   %tot      avg
 ──────────────────────────────────────────────────────────────────────────────
 Simulation                 1    46.5s  78.7%   46.5s   10.1GiB  91.1%  10.1GiB
 Init sampler               1    10.5s  17.7%   10.5s    759MiB  6.72%   759MiB
   Parse input files        1    7.62s  12.9%   7.62s    562MiB  4.97%   562MiB
   Create Context           1    1.36s  2.31%   1.36s    147MiB  1.30%   147MiB
   Create sampler           1   74.0ms  0.13%  74.0ms   3.88MiB  0.03%  3.88MiB
 Write results              1    2.10s  3.57%   2.10s    251MiB  2.22%   251MiB
 ──────────────────────────────────────────────────────────────────────────────
 ──────────────────────────────────────────────────────────────────────────────
                                       Time                   Allocations
                               ──────────────────────   ───────────────────────
       Tot / % measured:            8.53s / 100%            2.98GiB / 100%

 Section               ncalls     time   %tot     avg     alloc   %tot      avg
 ──────────────────────────────────────────────────────────────────────────────
 Simulation                 1    8.45s  99.0%   8.45s   2.97GiB  100%   2.97GiB
 Init sampler               1   79.9ms  0.94%  79.9ms   11.8MiB  0.39%  11.8MiB
   Parse input files        1   38.2ms  0.45%  38.2ms   4.80MiB  0.16%  4.80MiB
   Create Context           1   17.0ms  0.20%  17.0ms   6.78MiB  0.22%  6.78MiB
   Create sampler           1   27.5μs  0.00%  27.5μs   2.08KiB  0.00%  2.08KiB
 Write results              1   2.43ms  0.03%  2.43ms   72.1KiB  0.00%  72.1KiB
 ──────────────────────────────────────────────────────────────────────────────
```

where we see that the simulation only took 8.45s, almost half the time.
Running for four processes with

```
mpiexecjl --project -n 4 julia bin/qxrun.jl -d rqc_7_7_16.qx -o rqc_7_7_16_output.jld2 -m  -t -b 1
```

results in

```
 ──────────────────────────────────────────────────────────────────────────────
                                       Time                   Allocations
                               ──────────────────────   ───────────────────────
       Tot / % measured:          1128856s / 0.01%          11.0GiB / 86.6%

 Section               ncalls     time   %tot     avg     alloc   %tot      avg
 ──────────────────────────────────────────────────────────────────────────────
 Simulation                 1    55.6s  79.7%   55.6s   8.57GiB  89.7%  8.57GiB
 Init sampler               1    12.0s  17.2%   12.0s    759MiB  7.76%   759MiB
   Parse input files        1    8.58s  12.3%   8.58s    562MiB  5.74%   562MiB
   Create Context           1    1.60s  2.30%   1.60s    147MiB  1.50%   147MiB
   Create sampler           1   88.1ms  0.13%  88.1ms   3.88MiB  0.04%  3.88MiB
 Write results              1    2.11s  3.03%   2.11s    251MiB  2.57%   251MiB
 ──────────────────────────────────────────────────────────────────────────────
 ──────────────────────────────────────────────────────────────────────────────
                                       Time                   Allocations
                               ──────────────────────   ───────────────────────
       Tot / % measured:            4.63s / 100%            1.50GiB / 100%

 Section               ncalls     time   %tot     avg     alloc   %tot      avg
 ──────────────────────────────────────────────────────────────────────────────
 Simulation                 1    4.56s  98.5%   4.56s   1.48GiB  99.2%  1.48GiB
 Init sampler               1   68.9ms  1.49%  68.9ms   11.8MiB  0.77%  11.8MiB
   Parse input files        1   46.3ms  1.00%  46.3ms   4.80MiB  0.31%  4.80MiB
   Create Context           1   17.9ms  0.39%  17.9ms   6.78MiB  0.44%  6.78MiB
   Create sampler           1   22.6μs  0.00%  22.6μs   2.08KiB  0.00%  2.08KiB
 Write results              1   2.01ms  0.04%  2.01ms   72.1KiB  0.00%  72.1KiB
  ──────────────────────────────────────────────────────────────────────────────
```

which shows the simulation taking 4.56s, almost half again.
This shows that we achieve approximately linear speedup in the time taken for the actual simulation with increasing numbers of processes.
For some scaling to larger numbers of processes see the results included in [arXiv:2110.09894](https://arxiv.org/pdf/2110.09894.pdf).