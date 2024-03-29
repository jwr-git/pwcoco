#pragma once

#include <algorithm>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <sys/stat.h>

#include <Eigen/StdVector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "data.h"
#include "helper_funcs.h"

#ifdef PYTHON_INC
#include <Python.h>
#endif

#ifdef SINGLE_PRECISION
typedef Eigen::SparseMatrix<float, Eigen::ColMajor, long long> eigenSparseMat;
#else
typedef Eigen::SparseMatrix<double, Eigen::ColMajor, long long> eigenSparseMat;
#endif

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(eigenSparseMat);
using namespace Eigen;
using namespace std;

#ifdef SINGLE_PRECISION
typedef DiagonalMatrix<float, Dynamic, Dynamic> eigenDiagMat;
typedef MatrixXf eigenMatrix;
typedef VectorXf eigenVector;
typedef DynamicSparseMatrix<float> eigenDynSparseMat;
#else
typedef DiagonalMatrix<double, Dynamic, Dynamic> eigenDiagMat;
typedef MatrixXd eigenMatrix;
typedef VectorXd eigenVector;
typedef DynamicSparseMatrix<double> eigenDynSparseMat;
#endif


enum cond_type {
	CO_COND = 0,
	CO_JOINT,
};

struct conditional_dat {
	eigenSparseMat B;
	eigenSparseMat B_i; // Identity matrix of B
	eigenSparseMat B_N;
	eigenSparseMat B_N_i; // Identity matrix of B_N
	eigenVector D_N;
	eigenSparseMat Z;
	eigenSparseMat Z_N;

	// Initialiser
	conditional_dat() : B{ 0, 0 }, B_i{ 0, 0 }, B_N{ 0, 0 }, B_N_i{ 0, 0 }, D_N(0), Z{ 0, 0 }, Z_N{ 0, 0 } {}; // TODO All will be resized later which may be inefficient
};

class cond_analysis {
public:
	cond_analysis(double p_cutoff, double collinear, double ld_window, string out, double top_snp, double freq_thres, string name, bool cond_ssize, bool verbose);
	cond_analysis();

	bool coloc_ready() {
		return cond_passed;
	}

	string get_cond_name() {
		return cname;
	}

	size_t get_num_ind() {
		return num_ind_snps;
	}

	string get_ind_snp_name(size_t pos) {
		try {
			return ja_snp_name[ind_snps[pos]] + std::string(num_ind_snps == 1 ? "*" : "");
		}
		catch (...) {
			cout << "Independent SNP index out of bound." << endl;
			return "";
		}
	}

	enum coloc_type get_coloc_type() {
		return ctype;
	}

	void init_conditional(phenotype *pheno, reference *ref);
	void find_independent_snps(conditional_dat *cdat, reference *ref);
	void pw_conditional(int pos, bool out_cond, conditional_dat *cdat, reference *ref);

	// For coloc
	vector<string> snps_cond; /// SNP names
	vector<double> b_cond; /// Beta
	vector<double> se_cond; /// SE(beta)
	vector<double> maf_cond; /// Minor allele frequency
	vector<double> p_cond; /// P values
	vector<double> n_cond; /// Sample sizes
	vector<double> s_cond; /// Cases for case-control (TODO probably needs conditioned)

private:
	void match_gwas_phenotype(phenotype *pheno, reference *ref);

	void makex_eigenVector(size_t j, eigenVector &x, bool resize, reference *ref);
	bool init_b(const vector<size_t> &idx, conditional_dat *cdat, reference *ref);
	void init_z(const vector<size_t> &idx, conditional_dat *cdat, reference *ref);
	bool insert_B_Z(const vector<size_t> &idx, size_t pos, conditional_dat *cdat, reference *ref);
	void erase_B_and_Z(const vector<size_t> &idx, size_t erase, conditional_dat *cdat);
	void stepwise_select(vector<size_t> &selected, vector<size_t> &remain, conditional_dat *cdat, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref);

	bool select_entry(vector<size_t> &selected, vector<size_t> &remain, conditional_dat *cdat, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref);
	void selected_stay(vector<size_t> &select, conditional_dat *cdat, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, reference *ref);
	void massoc_conditional(const vector<size_t> &selected, vector<size_t> &remain, conditional_dat *cdat, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref);
	void massoc_joint(const vector<size_t> &idx, conditional_dat *cdat, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, reference *ref);

	void LD_rval(const vector<size_t> &idx, eigenMatrix &rval, conditional_dat *cdat);
	void LD_rval(const vector<size_t> &v1, const vector<size_t> &v2, eigenMatrix &rval, reference *ref);
	void sanitise_output(vector<size_t> &selected, vector<size_t> &remain, int pos, conditional_dat *cdat, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, reference *ref);
	void locus_plot(char *filename, char *datafile, char *to_save, char *snpname, double bp, double p, double pC);

	double a_ld_window; // Distance in kb after which SNPs are considered to be in LD

	size_t num_ind_snps; // Amount of SNPs selected after stepwise selection; can be 0 if the unconditional data are to be used or -1 so no SNPs are excluded in the conditional
	vector<size_t> ind_snps; // Position of SNPs selected after stepwise selection

	// Joint analysis related
	double a_collinear; // Collinearity check between SNPs
	int jma_snpnum_collinear;
	int jma_snpnum_backward;
	vector<string> ja_snp_name;
	eigenVector ja_freq;
	eigenVector ja_beta;
	eigenVector ja_beta_se;
	eigenVector ja_pval;
	eigenVector ja_chisq;
	eigenVector ja_N_outcome; // May very well be unused
	eigenVector ja_n_cases;

	eigenVector msx; 
	eigenVector msx_b; 
	eigenVector nD;

	bool cond_ssize; /// Whether to use conditional sample sizes or not
	vector<double> nsample; /// Note that this is not conditioned like nD
	vector<double> ncases; /// Note that this is not conditioned like nD

	string cname;

	string a_out;
	double a_top_snp;
	double a_p_cutoff;
	double a_freq_threshold;
	int num_snps;
	bool verbose;

	vector<size_t> to_include; /// SNP list to include in analysis after sanitising
	vector<size_t> to_include_bim; /// SNP bim locations
	vector<size_t> fam_ids_inc; /// Family IDs that are included in the analysis

	// Joint analysis related
	double jma_Ve;
	double jma_Vp; /// Phenotypic variance

	vector<size_t> remain_snps; // Remainder of SNPs after the stepwise selection process
	bool cond_passed; // Ready for coloc after conditional analysis
	vector<double> mu; // Calculated allele frequencies using fam data

	// Coloc related
	coloc_type ctype; // Type of coloc to use: cc or quant
};
