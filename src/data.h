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
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

#include "helper_funcs.h"

using namespace std;

#if defined(_MSC_VER)
#define fseek64 _fseeki64
#else
#define _FILE_OFFSET_BITS 64
#define fseek64 fseeko
#endif

enum class coloc_type : int {
	COLOC_NONE = 0,
	COLOC_QUANT,
	COLOC_CC,
};

class cond_analysis;

class phenotype {
public:
	phenotype(string name, double n, double n_case, double pve, bool no_freq, string pve_file);
	phenotype();

	void read_phenofile(string filename);
	void phenotype_clear();

	string get_phenoname() {
		return pheno_name;
	}

	double get_variance() {
		return pheno_variance;
	}

	vector<string> &get_snp_names() {
		return snp_name;
	}

	bool has_failed() {
		return failed;
	}

	enum coloc_type get_coloc_type() {
		return ctype;
	}

	void copy_frequencies(const vector<double> freqs);

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
	void read_phenotypic_variance_file(string pve_file);
	void calc_pheno_variance(string filename);
	void calc_pheno_variance(string filename, vector<double> freq, vector<double> beta, vector<double> se, vector<double> n);

	string pheno_name;
	double pheno_variance; /// Estimated phenotypic variance from summary stats
	double n_from_cmd; /// N passed from command line
	double n_case_from_cmd; /// N_cases passed from command line

	bool no_freq; // Dataset does not have frequencies; ignore freq-related functions
	bool failed; // Phenotype reading failed in some way
	coloc_type ctype; // Type of coloc to use: cc or quant
};

phenotype *init_pheno(string filename, string pheno_name, double n, double n_case, double pve, bool no_freq, string pve_file);

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
	vector<double> s1, s2; // Case to total sample size proportion
	coloc_type type1, type2;

private:
	map<size_t, size_t> snp_map; // Positions of SNPs in original datasets
};

class reference {
public:
	reference(string out, unsigned short chr);
	reference();
	void reference_clear();

	int read_bimfile(string bimfile);
	int read_famfile(string famfile);
	int read_bedfile(string bedfile);
	void parse_bed_data(char *buf, size_t i, vector<int> read_individuals);
	void bim_clear();
	void fam_clear();
	void match_bim(vector<string> &names, vector<string> &names2, bool keep_frequencies);
	void whole_bim();
	void reset_vectors();

	int filter_snp_maf(double maf);
	void sanitise_list();
	void pair_fam();
	void get_read_individuals(vector<int> &read_individuals);
	void update_inclusion(const vector<size_t> idx, const vector<string> snps);

	void includes_clear() {
		to_include.clear();
	}

	bool has_failed() {
		return failed;
	}

	bool is_ready() {
		return read;
	}

	vector<double> get_freqs() {
		return mu_m;
	}

	// From .bim
	vector<string> bim_snp_name; /// SNP names
	vector<string> ref_A; /// Reference allele
	vector<size_t> to_include; /// SNP list to include in analysis after sanitising
	vector<size_t> to_include_bim; /// Positions in the original bim file
	map<string, size_t> snp_map; /// Maps rsID/SNP identifer to vector position
	vector<string> bim_allele1; /// A1
	vector<string> bim_allele2; /// A2
	vector<unsigned short> bim_chr; /// Chromosome
	vector<int> bim_bp; /// BP position
	vector<size_t> bim_read_pos; /// Read position in the .bim file

	// Unaltered vectors
	vector<string> bim_snp_name_m; /// Unaltered SNP names
	vector<string> bim_allele1_m; /// A1
	vector<string> bim_allele2_m; /// A2
	vector<unsigned short> bim_chr_m; /// Chromosome
	vector<int> bim_bp_m; /// BP position

	// From .fam
	vector<size_t> fam_ids_inc; /// Family IDs that are included in the analysis
	vector<double> mu; /// Calculated allele frequencies using fam data

	// From .bed file
	vector<vector<bool>> bed_snp_1;
	vector<vector<bool>> bed_snp_2;

private:
	string a_out;
	unsigned short a_chr;
	bool failed; // Reference files failed to read in some way
	bool read; // Reference files have already been read and cleaned, if true.

	// From .bim file
	vector<size_t> bim_og_pos; /// Position in the .bim file
	vector<double> bim_genet_dst; /// Distance 
	// Extra helper info
	size_t start_snps, end_snps; /// Location of first read and last read SNP in the reference panel
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

	// Unaltered vectors for frequencies
	vector<vector<bool>> bed_snp_1_m;
	vector<vector<bool>> bed_snp_2_m;
	vector<double> mu_m; /// Calculated allele frequencies using fam data

	vector<vector<bool>>snp_1;
	vector<vector<bool>>snp_2;
};
