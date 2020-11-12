# Pair-Wise Conditional analysis and Colocalisation analysis (PWCoCo)
A C++ implementation of the PWCoCo algorithm described by Zheng, et al in their paper, [Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases](https://doi.org/10.1038/s41588-020-0682-6). 

This tool integrates methods from [GCTA-COJO](https://cnsgenomics.com/software/gcta/#Overview) and the [coloc](https://chr1swallace.github.io/coloc/index.html) R package.

## Requirements
- Cmake >= 3.10
- Eigen >= 3.3.7
- spdlog >= 1.8.1

These libraries are bundled in `/include/` at the required versions for ease of the user.

## How to Build
Currently, only Unix and Windows are supported.

### Unix
To build on Unix systems, follow the code below:
```
git clone https://github.com/jwr-git/pwcoco

mkdir build
cd build
cmake ..
make
```

If building on the University of Bristol's HPC, load the module `languages/gcc-9.1.0` and ensure this is the **only** gcc module loaded. This should be all you need to build the program.

### Windows

A .sln file is provided for Visual Studio 2019. Be aware that there may be issues using the Eigen library, as packaged, on Windows machines.

## How to Use
PWCoCo is a command-line program. Here is a list of accepted flags with a description of each one:

**Required**
- `--bfile` - specifies the location of the reference dataset, normally from Plink, in the bed/bim/fam formats. Each of the bed/bim/fam files should have the same name and in the same directory.
- `--phen1_file` - first phenotype file, file ending does not matter.
- `--phen2_file` - second phenotype file, file ending does not matter.

For acceptable formats for these files, please see below.

**Optional**
- `--log` - specifies log name, default is "pwcoco_log.txt" and will save in the same folder from where the program is run.
- `--out` - prefix for the result files, default is "pwcoco_out".
- `--p_cutoff` - P value cutoff for SNPs to be selected by the stepwise selection process, default is 5e-8. 
- `--chr` - when reading the reference files, the program will limit the analysis to those SNPs on this chromosome. 
- `--top_snp` - maximum number of SNPs that may be selected by the stepwise selection process, default is 1e10, i.e. a lot.
- `--ld_window` - distance (in kb) that, when exceeded, is assumed for SNPs to be in total LE, default is 1e7.
- `--collinear` - threshold that, when exceeded, determines if SNPs are collinear, default is 0.9.
- `--maf` - filters SNPs from the reference dataset according to this threshold, default is 0.1.
- `--freq_threshold` - SNPs in the phenotype datasets which differ by more than this amount in the reference dataset will be excluded, default is 0.2.
- `--init_h4` - PWCoCo will run an initial colocalisation on the unconditioned dataset. If the H4 for this analysis reaches this threshold, the program will terminate early. Default is 80 (i.e. 80%). Set to 0 if you would like the program to always continue regardless of the initial colocalisation result.
- `--out_cond` - true/false: would you like for the conditioned data to be saved as text files as well?

## Example
Example files will be provided soon so that a full analysis can be run. Instead, here is an example command to run the analysis:

`pwcoco --bfile "../../1kg_plink/chr5" --phen1_file "tgfbi/TGFBI_exposure.txt" --phen2_file "tgfbi/TGFBI_outcome.txt" --out "tgfbi/res" --chr 5 --maf 0.01 --out_cond 1`

### Input File Formats
The bfile, or reference files, must be in Plink format. This means a .bed, .bim and .fam file in the same directory with the same name.

Phenotype files do not require a certain file format. Instead, they _must_ follow this structure:

`SNP	effect_allele	other_allele	effect_allele_freq	beta	se	p	n`

Column names do not matter, only the order of the data.

### Output File Formats
The program by default will output a file with the ending `.coloc` which contains the results for each of the colocalisation analyses run:

`SNP1	SNP2	H0	H1	H2	H3	H4	log_abf_all`

If the data has been unconditioned, then the SNP column will contain "unconditioned" instead of a SNP name. Please note that output files **are not** deleted or overwritten between runs. That means if you run the program twice with the same output file name, results will be appended to the output file. SNPs correspond to the same numbered phenotype file, e.g. SNP1 comes from phen1_file.

The program will also output two extra files (ending in `.included`) listing the SNPs included in the analysis from each dataset. Finally, if the `out_cond` flag is set to true, the program will output the conditioned results after the conditional analysis in files ending with `.cojo`:

`Chr	SNP	bp	refA	freq	b	se	p	n	freq_geno	bC	bC_se	pC`

The columns `freq`, `b`, `se` and `p` should be unaltered from the original dataset. The final columns will be post-conditional analysis. 

## To Do
- Allow more functions to work asynchronously.
- Rewrite the .bed reading function for efficiency.
- Provide example files.
- Clean up input flags.
