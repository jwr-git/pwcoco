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


class cond_analysis {
public:
	cond_analysis(double p_cutoff, double collinear, int ld_window, string out, bool verbose, int top_snp);
	cond_analysis();

	void init_conditional(phenotype *pheno, reference *ref);
	void match_gwas_phenotype(phenotype *pheno, reference *ref);
	void massoc(cond_analysis *p_analysis);

	void makex_eigenVector(int j, eigenVector &x, reference *ref);
	bool init_b(const vector<int> &idx, reference *ref);
	void init_z(const vector<int> &idx, reference *ref);
	bool insert_B_Z(const vector<int> &idx, int pos, reference *ref);
	void stepwise_select(vector<int> &selected, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref);

	bool select_entry(vector<int> &selected, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref);
	void massoc_conditional(const vector<int> &selected, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref);

	int a_ld_window; // Distance in kb after which SNPs are considered to be in LD

	// Joint analysis related
	double jma_Ve; /// What am I?
	double a_collinear; // Collinearity check between SNPs
	int jma_snpnum_collinear;
	eigenVector ja_freq;
	eigenVector ja_beta;
	eigenVector ja_beta_se;
	eigenVector ja_pval;
	eigenVector ja_N_outcome;

	eigenVector msx; /// What am I?
	eigenVector msx_b; /// What am I?
	eigenVector nD; /// What am I?

	eigenSparseMat B;
	eigenMatrix B_i; // Identity matrix of B
	eigenSparseMat B_N;
	eigenMatrix B_N_i; // Identity matrix of B_N
	eigenVector D_N;
	eigenSparseMat Z;
	eigenSparseMat Z_N;

private:
	string a_out;
	int a_top_snp;
	double a_p_cutoff;
	bool a_verbose;
	int num_snps;

	// Joint analysis related
	double jma_Vp; /// What am I?
};
