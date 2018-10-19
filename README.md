# CNV Trees

## Requirements

* C++ compiler that supports C++14 standards (e.g. `gcc>=6.3.0`, `clang>=5.0.0)`)
* CMake version > 3.8

## Installation

1. Open your terminal in your preferred directory and clone this project
```bash
$ git clone https://github.com/anilbey/sc-dna
```
2. Enter the project directory
```bash
$ cd sc-dna
```
3. Create and enter the build directory
```bash
$ mkdir build && cd build
```
4. Compile the program with cmake
```bash
$ cmake ..
```
5. Build the executable
```bash
$ make
```
That's it! :octocat:

Three executables (simulation, sc-dna and test) will be generated after the installation.


## Simulation
Simulates the count matrix. Outputs the count matrix, region sizes, ground truth and the tree that generated the data.

#### Simulation parameters

| Parameter name | Description | Default value |
| ---- | -------- | --- | 
| **n_bins** | Number of bins of the input matrix | 10000 |
| **n_cells** | Number of cells | 500 |
| **n_nodes** | Number of nodes of the tree | 50 |
| **n_regions** | Number of regions | 50 |
| **n_iters** | Number of iterations | 10000 |
| **n_reads** | Number of reads per cell | 10000 |
| **ploidy** | The ploidy information | 2 (diploid, human) |
| **verbosity** | Verbosity of the programme, 1 provides standard output, 2 also writes several files | 0 |
| **seed** | Seed | - |
| **postfix** | Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system | "" |
| **print_precision** | The precision points of the score values to be printed | 16 |

#### *Sample run* :
```shell
$ ./simulations --print_precision 32 --n_bins 100 --n_regions 10 --n_nodes 10 --n_reads 100000 --verbosity 2 --ploidy 2 --n_cells 400 --postfix 9 --seed 42   
```

## Inference

## Test


