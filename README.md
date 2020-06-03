<div align="center">
  <img src="logo.png" >
</div>

[![circleci](https://img.shields.io/circleci/project/github/cbg-ethz/SCICoNE/master.svg?&logo=circleci&color=blueviolet)](https://circleci.com/gh/cbg-ethz/SCICoNE/)
[![C++ Standard](https://img.shields.io/badge/c++-14-blueviolet?style=flat&logo=c%2B%2B)](https://en.wikipedia.org/wiki/C%2B%2B14)
[![License](https://img.shields.io/:license-GPLv3-blueviolet.svg?style=flat&logo=gnu)](http://www.gnu.org/licenses/gpl-3.0.html)


## About
Single-cell copy number calling and event history reconstruction.

A statistical model and MCMC algorithm tailored to single-cell copy
number profiling from shallow whole-genome DNA sequencing data. SCICoNE reconstructs the history of copy number events in the tumour and uses these evolutionary relationships to identify the copy number profiles of the individual cells.

## Requirements
* C++ compiler that supports C++14 standards (e.g. `gcc>=5.2.0`, `clang>=5.0.0)`)
* CMake >= 3.9
* Boost >= 1.6.x
* OpenMP >= 4.5
* NLopt >= 2.6.2

## Installation
Once the requirements are in place, downloading and installing SCICoNE takes about 5 minutes.
```bash
git clone git@github.com:cbg-ethz/SCICoNE.git # Clone the repository
cd SCICoNE
mkdir build && cd build                       # Create and enter the build directory
cmake ..                                      # Compile the program with cmake
make                                          # Build the executables
```

### Python wrapper
We also provide a Python 3 wrapper to facilitate plotting and easily integrate SCICoNE with other data analysis tools. This interface runs the C++ binaries and reads the outputs into `numpy` arrays.
```bash
cd pyscicone
pip install .
```

## Quick start
SCICoNE takes a read counts matrix of cells by genomic bins and outputs the copy number profile of each cell and the underlying event history. An introductory notebook showcasing the SCICoNE workflow and API is available [here](https://github.com/cbg-ethz/SCICoNE/blob/master/notebooks/tutorial.ipynb).

## C++ reference
### Breakpoint detection
Find the breakpoints defining copy number regions from a cells by bins matrix.

| Parameter name | Description | Default value |
| ---- | -------- | --- |
| **d_matrix_file** | Path to the counts matrix file, delimiter: ' ', line separator: '\n'  | "" |
| **n_cells** | Number of cells in the input matrix | - |
| **n_bins** | Number of bins in the input matrix | - |
| **min_cells** | Minimum number of cells to consider for a bin being a breakpoint | 4 |
| **window_size** | Size of the window used in breakpoint detection | 10 |
| **threshold** | Breakpoint calling sensitivity. If -1, stop after computing the LR | 3 |
| **bp_limit** | Maximum number of breakpoints | 300 |
| **compute_lr** | Boolean indicator of wether the per bin cell-wise breakpoint evidence should be computed | True |
| **lr_file** | Path to a matrix containing the evidence for breakpoint at each bin in each cell | "" |
| **sp_file** | Path to a vector containing the combined evidence for breakpoint at each bin across all cells | "" |
| **compute_sp** | Boolean indicator of wether the per bin breakpoint evidence should be computed | True |
| **evaluate_peaks** | Boolean indicator of wether to evaluate peaks and call breakpoints | True |
| **input_breakpoints_file** | Path to file indicating bins which correspond to known breakpoints that must be included. These may be, for example, the chromosome bin coordinates. | ""|
| **verbosity** | Verbosity of the programme, 0 is non-verbose setting, 1 creates the debug files, 2 writes the logs as well | 0 |
| **postfix** | Postfix to be added to the output files, this is useful when you are running multiple simulations through a workflow management system | "" |
| **print_precision** | The precision points of the score values to be printed | 16 |


#### *Sample run* :
```shell
$ ./breakpoint_detection --d_matrix_file ./d_mat.txt --n_bins 1000 --n_cells 400 --window_size 10 --threshold 3.0 --bp_limit 300 --verbosity=1 --evaluate_peaks=True
```

### Inference
Finds the maximum a posteriori tree given a cells by regions matrix.

| Parameter name | Description | Default value |
| ---- | -------- | --- |
| **region_sizes_file** | Path to the file containing the region sizes, each line contains one region size. Segmentation is performed if the region sizes file is not specified | "" |
| **d_matrix_file** | Path to the counts matrix file, delimiter: ' ', line separator: '\n'  | "" |
| **n_iters** | Number of iterations | 10000 |
| **n_cells** | Number of cells in the input matrix | - |
| **n_regions** | Number of regions in the input matrix | - |
| **ploidy** | The ploidy information | 2 (diploid, human) |
| **verbosity** | Verbosity of the programme, 0 is non-verbose setting, 1 creates the debug files, 2 writes the inference logs as well, 3 writes the tree logs on top | 0 |
| **seed** | Seed | - |
| **postfix** | Postfix to be added to the output files, this is useful when you are running multiple simulations through a workflow management system | "" |
| **print_precision** | The precision points of the score values to be printed | 16 |
| ---- | parameters for the random initialised tree | --- |
| **n_nodes** | the number of nodes in the random initialised tree | 50 |
| **lambda_r** | lambda param for the poisson that generates the number of regions | 0.1 |
| **lambda_c** | lambda param for the poisson that generates the copy number state of a region | 0.2 |

#### *Sample run* :
```shell
$ ./inference --n_cells 400 --n_regions 10 --n_iters 100 --n_nodes 10 --ploidy 2 --verbosity 2 --seed 42 --d_matrix_file ./30_d_mat.txt --region_sizes_file ./30_region_sizes.txt
```

### Simulation
Simulates the count matrix. Outputs the count matrix, region sizes, ground truth CNVs and the tree that generated the data.

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
$ ./simulation --n_cells 400 --n_bins 100 --n_regions 10 --n_nodes 10 --n_reads 100000 --ploidy 2 --verbosity 2  --seed 42 --print_precision 32
```

### Test
Runs the validation tests and writes the results to the standard output and error streams.
#### *Sample run* :
```shell
$ ./tests
```
