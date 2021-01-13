#pragma once

#include "data.h"
#include "conditional.h"
#include "helper_funcs.h"

using namespace std;

enum hypothesis {
	H0 = 0,
	H1,
	H2,
	H3,
	H4,
};

class coloc_analysis {
public:
	coloc_analysis(mdata *mdat, string out, double pval1, double pval2, double pval3);
	coloc_analysis();
	~coloc_analysis();

	void init_coloc();
	void init_coloc(string snp1, string snp2);

	vector<double> pp_abf; // Results from colocalisation
	
private:
	bool perform_coloc();
	bool estimate_bf(const vector<double> beta, const vector<double> se, const vector<double> freq, const vector<double> n, coloc_type type, vector<double> *ABF);
	void combine_abf(size_t abf_size);
	void results_to_file(string s1, string s2);

	mdata *matched;
	string outfile;

	// Coloc stuff
	double p1; // Prior probability of SNP associated with trait 1
	double p2; // Prior probability of SNP associated with trait 2
	double p3; // Prior probability of SNP associated with both traits
	map<size_t, size_t> snp_map; // 1st refers to snps1, 2nd refers to snps2

	vector<double> *ABF_1, *ABF_2;
	vector<double> ABF_sum;
	double h0, h1, h2, h3, h4;
	double log_ABF_sum;
	double log_abf_all;
};
