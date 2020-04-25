<div align="center">
  <img src="logo.png" >
</div>


[![CircleCI](https://circleci.com/gh/anilbey/SCICoNE.svg?style=svg&circle-token=ccc7efcfe007311847eb0d7260382125518d68ca)](https://circleci.com/gh/anilbey/SCICoNE)
[![C++ Standard](https://img.shields.io/badge/c++-14-blueviolet?style=flat&logo=c%2B%2B)](https://en.wikipedia.org/wiki/C%2B%2B14)
[![License](https://img.shields.io/:license-GPLv3-blueviolet.svg?style=flat&logo=gnu)](http://www.gnu.org/licenses/gpl-3.0.html)


## About

Single-cell copy number calling and event history reconstruction.
A statistical model and MCMC algorithm tailored to single-cell copy
number profiling from shallow whole-genome DNA sequencing data. SCICoNE reconstructs the
history of copy number events in the tumour and uses these evolutionary relationships to identify
the copy number profiles of the individual cells. 

## Requirements

* C++ compiler that supports C++14 standards (e.g. `gcc>=5.2.0`, `clang>=5.0.0)`)
* CMake version > 3.8
* Boost >= 1.6.x

## Installation

1. Clone the repository

2. Enter the project directory
```bash
$ cd SCICoNE
```
3. Create and enter the build directory
```bash
$ mkdir build && cd build
```
4. Compile the program with cmake
```bash
$ cmake ..
```

or more specifically
```bash
$ cmake -DCMAKE_BUILD_TYPE=Release ..
```

5. Build the executable
```bash
$ make
```
That's it! :octocat:

Multiple executables (such as simulation, inference, test, score) will be created after the installation.

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
| **verbosity** | Verbosity of the programme, 0 is non-verbose setting, 1 creates the debug files, 2 writes the inference logs as well, 3 writes the tree logs on top | 0 |
| **seed** | Seed | - |
| **postfix** | Postfix to be added to the output files, this is useful when you are running multiple simulations through a work flow management system | "" |
| **print_precision** | The precision points of the score values to be printed | 16 |

#### *Sample run* :
```shell
$ ./simulation --print_precision 32 --n_bins 100 --n_regions 10 --n_nodes 10 --n_reads 100000 --verbosity 2 --ploidy 2 --n_cells 400 --postfix 9 --seed 42   
```

## Inference
Finds the maximum likelihood tree given cellsxregions matrix or the simulated matrix with params specified.

#### Inference parameters

| Parameter name | Description | Default value |
| ---- | -------- | --- | 
| **region_sizes_file** | Path to the file containing the region sizes, each line contains one region size. Segmentation is performed if the region sizes file is not specified | "" |
| **d_matrix_file** | Path to the counts matrix file, delimiter: ' ', line separator: '\n'  | "" |
| **n_bins** | Number of bins in the input matrix | - |
| **n_iters** | Number of iterations | 10000 |
| **n_cells** | Number of cells in the input matrix | - |
| **n_reads** | The number of reads per cell. This has no functionality besides being appended to the output files | -1 |
| **ploidy** | The ploidy information | 2 (diploid, human) |
| **verbosity** | Verbosity of the programme, 0 is non-verbose setting, 1 creates the debug files, 2 writes the inference logs as well, 3 writes the tree logs on top | 0 |
| **seed** | Seed | - |
| **postfix** | Postfix to be added to the output files, this is useful when you are running multiple simulations through a workflow management system | "" |
| **print_precision** | The precision points of the score values to be printed | 16 |
| ---- | parameters for the random initialised tree | --- | 
| **n_nodes** | the number of nodes in the random initialised tree | 50 |
| **lambda_r** | lambda param for the poisson that generates the number of regions | 0.1 |
| **lambda_c** | lambda param for the poisson that generates the copy number state of a region | 0.2 |
| **max_region_size** | the maximum size that a region can have | 10 |


#### *Sample run* :
```shell
$ ./inference --n_reads 100000 --print_precision 16 --n_nodes 10 --n_bins 10000 --n_iters 100 --n_cells 500 --verbosity 2 --ploidy 2 --seed 42 --postfix d31 --d_matrix_file ./30_d_mat.txt --region_sizes_file ./30_region_sizes.txt 
```

## Test
Runs the validation tests and writes the results to the standard output and error streams.
#### *Sample run* :
```shell
$ ./tests
```

Logo url: https://logomakr.com/9tSL7m
