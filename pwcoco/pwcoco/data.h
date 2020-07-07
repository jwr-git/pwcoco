#pragma once

#include <algorithm>
#include <bitset>
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/algorithm/string.hpp>

#include "helper_funcs.h"

using namespace std;

class phenotype {
public:
	phenotype(string name);
	phenotype();

	void read_phenofile(string filename);
	void phenotype_clear();
	double calc_variance();

	// From phenotype file
	vector<string> snp_name;
	vector<string> allele1;
	vector<string> allele2;
	vector<double> freq;
	vector<double> beta;
	vector<double> se;
	vector<double> pval;
	vector<double> n;
	vector<double> n_case;
	vector<double> Vp_v;
	vector<double> mu;

	vector<size_t> matched_idx; /// Indicies of SNPs that have been matched

private:
	string pheno_name;
	double pheno_variance; /// Estimated phenotypic variance from summary stats
};

phenotype *init_exposure(string filename, string pheno_name);

class mdata {
public:
	mdata(phenotype *ph1, phenotype *ph2);
	mdata();

	// Data from datasets
	// These are matched!
	vector<string> snps1, snps2;
	vector<double> betas1, betas2;
	vector<double> ses1, ses2;
	vector<double> pvals1, pvals2;
	vector<double> mafs1, mafs2;
	vector<double> ns1, ns2; // Must these be double? Can they be int?

private:
	map<size_t, size_t> snp_map; // Positions of SNPs in original datasets
};

class reference {
public:
	reference(string out, unsigned short chr, bool verbose);
	reference();
	void reference_clear();

	void read_bimfile(string bimfile, phenotype *pheno);
	void read_famfile(string famfile);
	void read_bedfile(string bedfile);
	void bim_clear();
	void fam_clear();

	void filter_snp_maf(double maf);
	void sanitise_list();
	void pair_fam();
	void calculate_allele_freq();
	void get_read_individuals(vector<int> &read_individuals);
	void get_read_snps(vector<int> &read_snps);

	// From .bim
	vector<string> bim_snp_name; /// SNP names
	vector<string> ref_A; /// Reference allele
	vector<size_t> to_include; /// SNP list to include in analysis after sanitising
	map<string, size_t> snp_map; /// Maps rsID/SNP identifer to vector position
	vector<string> bim_allele1; /// A1
	vector<string> bim_allele2; /// A2
	vector<int> bim_chr; /// Chromosome
	vector<int> bim_bp; /// BP position

	// From .fam
	vector<size_t> fam_ids_inc; /// Family IDs that are included in the analysis
	vector<double> mu; /// Calculated allele frequencies using fam data

	// From .bed file
	vector<vector<bool>> bed_snp_1;
	vector<vector<bool>> bed_snp_2;

private:
	string a_out;
	bool a_verbose;
	unsigned short a_chr;

	// From .bim file
	vector<double> bim_genet_dst; /// Distance 
	// Extra helper info
	vector<size_t> og_indices; /// Original indices for the SNP names in bim_snp_name - use this to get information of a SNP from the other, unsorted vectors
	size_t num_snps; /// Number of SNPs in analysis
	vector<string> other_A; /// Other allele

	// From .fam file
	vector<string> fam_fid; /// FID or family ID
	vector<string> fam_iid; /// IID or within-family ID
	vector<string> fam_fa_id; /// Father ID
	vector<string> fam_mo_id; /// Mother ID
	vector<unsigned short> fam_sex; /// Sex '1' = male, '2' = female, '0' = unknown
	vector<unsigned short> fam_pheno; /// Phenotype value '1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control
	size_t individuals; /// Number of individuals read from the .fam file
	map<string, int> fam_map; /// Mapping between FIDs and IIDs
};
