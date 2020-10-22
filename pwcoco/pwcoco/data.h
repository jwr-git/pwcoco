#pragma once

#include <algorithm>
#include <bitset>
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <omp.h>
#include <string>
#include <sstream>
#include <vector>

#include "helper_funcs.h"

using namespace std;

class cond_analysis;

class phenotype {
public:
	phenotype(string name);
	phenotype();

	void read_phenofile(string filename);
	void phenotype_clear();
	double calc_variance(vector<size_t> idx);

	string get_phenoname() {
		return pheno_name;
	}

	double get_variance() {
		return pheno_variance;
	}

	vector<string> &get_snp_names() {
		return snp_name;
	}

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

phenotype *init_pheno(string filename, string pheno_name);

class mdata {
public:
	mdata(phenotype *ph1, phenotype *ph2);
	mdata(cond_analysis *ca1, cond_analysis *ca2);
	mdata(cond_analysis *ca, phenotype *ph);
	mdata();

	vector<string> &get_snp_list() {
		return snps1;
	}

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

	void read_bimfile(string bimfile);
	void read_famfile(string famfile);
	void read_bedfile(string bedfile);
	void bim_clear();
	void fam_clear();
	void match_bim(vector<string> &names, vector<string> &names2);

	void filter_snp_maf(double maf);
	void sanitise_list();
	void pair_fam();
	void calculate_allele_freq();
	void get_read_individuals(vector<int> &read_individuals);
	void update_inclusion(const vector<size_t> idx, const vector<string> snps);

	void includes_clear() {
		to_include.clear();
	}

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
	vector<size_t> bim_og_pos; /// Original position in the .bim file
	vector<double> bim_genet_dst; /// Distance 
	// Extra helper info
	size_t num_snps; /// Number of SNPs in analysis
	size_t num_snps_matched; /// Number of SNPs in analysis after matching
	vector<string> other_A; /// Other allele
	vector<bool> read_snps; /// SNPs that are read and included by the bim file after phenotype matching
	size_t tot_read_snps; /// Number of these SNPs

	// From .fam file
	vector<string> fam_fid; /// FID or family ID
	vector<string> fam_iid; /// IID or within-family ID
	vector<string> fam_fa_id; /// Father ID
	vector<string> fam_mo_id; /// Mother ID
	vector<unsigned short> fam_sex; /// Sex '1' = male, '2' = female, '0' = unknown
	vector<unsigned short> fam_pheno; /// Phenotype value '1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control
	size_t individuals; /// Number of individuals read from the .fam file
	map<string, size_t> fam_map; /// Mapping between FIDs and IIDs
};
