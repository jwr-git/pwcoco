#pragma once

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "data.h"
#include "helper_funcs.h"

using namespace Eigen;
using namespace std;

typedef VectorXd eigenVector;

#ifdef SINGLE_PRECISION
typedef DiagonalMatrix<float, Dynamic, Dynamic> eigenDiagMat;
typedef MatrixXf eigenMatrix;
typedef VectorXf eigenVector;
typedef SparseMatrix<float> eigenSparseMat;
typedef DynamicSparseMatrix<float> eigenDynSparseMat;
#else
typedef DiagonalMatrix<double, Dynamic, Dynamic> eigenDiagMat;
typedef MatrixXd eigenMatrix;
typedef VectorXd eigenVector;
typedef SparseMatrix<double> eigenSparseMat;
typedef DynamicSparseMatrix<double> eigenDynSparseMat;
#endif

enum cond_type {
	CO_COND = 0,
	CO_JOINT,
};

class cond_analysis {
public:
	cond_analysis(double p_cutoff, double collinear, double ld_window, string out, bool verbose, double top_snp, bool actual_geno, double freq_thres);
	cond_analysis();

	void init_conditional(phenotype *pheno, reference *ref);
	void match_gwas_phenotype(phenotype *pheno, reference *ref);
	void massoc(reference *ref, string snplist);

	vector<size_t> read_snplist(string snplist, vector<size_t> &remain, reference *ref);
	void makex_eigenVector(size_t j, eigenVector &x, bool resize, reference *ref);
	bool init_b(const vector<size_t> &idx, reference *ref);
	void init_z(const vector<size_t> &idx, reference *ref);
	bool insert_B_Z(const vector<size_t> &idx, size_t pos, reference *ref);
	void erase_B_and_Z(const vector<size_t> &idx, size_t erase, reference *ref);
	void stepwise_select(vector<size_t> &selected, vector<size_t> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref);

	bool select_entry(vector<size_t> &selected, vector<size_t> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref);
	void selected_stay(vector<size_t> &select, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, reference *ref);
	void massoc_conditional(const vector<size_t> &selected, vector<size_t> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref);
	void massoc_joint(const vector<size_t> &idx, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, reference *ref);

	double massoc_calcu_Ve(const vector<size_t> &selected, eigenVector &bJ, eigenVector &b);
	void LD_rval(const vector<size_t> &idx, eigenMatrix &rval);
	void sanitise_output(vector<size_t> &selected, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, eigenMatrix &rval, cond_type ctype, reference *ref);

	double a_ld_window; // Distance in kb after which SNPs are considered to be in LD

	size_t sw_snps; // Amount of SNPs selected after stepwise selection

	// Joint analysis related
	double a_collinear; // Collinearity check between SNPs
	int jma_snpnum_collinear;
	int jma_snpnum_backward;
	eigenVector ja_freq;
	eigenVector ja_beta;
	eigenVector ja_beta_se;
	eigenVector ja_pval;
	eigenVector ja_N_outcome;

	eigenVector msx; /// What am I?
	eigenVector msx_b; /// What am I?
	eigenVector nD; /// What am I?

	// For coloc
	vector<string> snps_cond; /// SNP names
	vector<double> b_cond; /// Beta
	vector<double> se_cond; /// SE(beta)
	vector<double> maf_cond; /// Minor allele frequency
	vector<double> p_cond; /// P values
	vector<double> n_cond; /// Sample sizes

	eigenSparseMat B;
	eigenMatrix B_i; // Identity matrix of B
	eigenSparseMat B_N;
	eigenMatrix B_N_i; // Identity matrix of B_N
	eigenVector D_N;
	eigenSparseMat Z;
	eigenSparseMat Z_N;

private:
	string a_out;
	double a_top_snp;
	double a_p_cutoff;
	bool a_verbose;
	bool a_actual_geno;
	double a_freq_threshold;
	int num_snps;

	// Joint analysis related
	double jma_Ve; /// What am I?
	double jma_Vp; /// Phenotypic variance
	double GC_val; /// What am I? TODO
};
