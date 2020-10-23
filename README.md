# Pair-Wise Conditional analysis and Colocalisation analysis (PWCoCo)
A C++ implementation of the PWCoCo algorithm described by Zheng, et al in their paper, [Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases](https://doi.org/10.1038/s41588-020-0682-6). 

This tool integrates methods from [GCTA-COJO](https://cnsgenomics.com/software/gcta/#Overview) and the [coloc](https://chr1swallace.github.io/coloc/index.html) R package.

## Requirements
- Eigen >= 3.3.7
- spdlog >= 1.8.1

These libraries are bundled in `/include/` at the required versions for ease of the user.

## How to Build
Currently, only Unix and Windows are supported.

- Unix
To build on Unix systems, follow the code below:
```
git clone https://github.com/jwr-git/pwcoco

mkdir build
cd build
cmake ..
make
```

- Windows

A .sln file is provided for Visual Studio 2019.
