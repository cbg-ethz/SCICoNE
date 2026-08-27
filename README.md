<div align="center">
  <img src="logo.png" >
</div>

[![Build](https://github.com/cbg-ethz/SCICoNE/actions/workflows/build.yml/badge.svg)](https://github.com/cbg-ethz/SCICoNE/actions/workflows/build.yml)
[![C++ Standard](https://img.shields.io/badge/c++-14-blueviolet?style=flat&logo=c%2B%2B)](https://en.wikipedia.org/wiki/C%2B%2B14)
[![License](https://img.shields.io/:license-GPLv3-blueviolet.svg?style=flat&logo=gnu)](http://www.gnu.org/licenses/gpl-3.0.html)


Single-cell copy number calling and event history reconstruction.

A statistical model and MCMC algorithm tailored to single-cell copy
number profiling from shallow whole-genome DNA sequencing data. SCICoNE reconstructs the history of copy number events in the tumour and uses these evolutionary relationships to identify the copy number profiles of the individual cells.

## Quick start
SCICoNE takes a read counts matrix of cells by genomic bins and outputs the copy number profile of each cell and the underlying event history.
* [Command line interface tutorial](https://github.com/cbg-ethz/SCICoNE/blob/master/docs/tutorial.md)
* [Python package tutorial](https://github.com/cbg-ethz/SCICoNE/blob/master/docs/notebooks/tutorial.ipynb)
* [Example run on 10X Genomics data using the Python package](https://github.com/cbg-ethz/SCICoNE/blob/master/docs/notebooks/10X_example.ipynb)

## Installation
### conda (recommended)
```bash
conda install -c bioconda -c conda-forge scicone
```
This installs the command line executables (`scicone-simulation`,
`scicone-breakpoint_detection`, `scicone-inference`, `scicone-score`) and the
`scicone` Python package, which finds them on your `PATH`.

### From source
Requirements:
* C++ compiler that supports C++14 standards (e.g. `gcc>=5.2.0`, `clang>=5.0.0)`)
* CMake >= 3.16
* Boost >= 1.6.x (headers only)
* OpenMP >= 4.5
* NLopt >= 2.6.2

Once the requirements are in place, downloading and installing SCICoNE takes about 5 minutes.
```bash
git clone https://github.com/cbg-ethz/SCICoNE.git # Clone the repository
cd SCICoNE
cmake -S scicone -B build                     # Configure the build
cmake --build build                           # Build the executables
cmake --install build --prefix <install_dir>  # Optional: install them onto your PATH
```
The executables are named `scicone-simulation`, `scicone-breakpoint_detection`,
`scicone-inference` and `scicone-score`, and are left in `build/` if you do not
install them.

### Python package
We also provide a Python 3 package to facilitate plotting and easily integrate SCICoNE with other data analysis tools. This interface runs the C++ binaries and reads the outputs into `numpy` arrays. Even if you don't want to use the complete package, we recommend you install it to facilitate usage of the C++ command line interface. It comes with the conda package; to install it separately:
```bash
pip install pyscicone/
```
It looks for the executables on your `PATH` first. If you did not install them,
point it at your build directory, either with
`SCICoNE(binary_path="/path/to/SCICoNE/build")` or by setting `$SCICONE_BINARY_PATH`.
