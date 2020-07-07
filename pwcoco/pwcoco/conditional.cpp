#include "conditional.h"

/*
 * cond_analysis constructor
 */
cond_analysis::cond_analysis(double p_cutoff, double collinear, double ld_window, string out, bool verbose, double top_snp)
{
	a_out = out;
	a_p_cutoff = p_cutoff;
	a_collinear = collinear;
	a_ld_window = ld_window;
	a_verbose = verbose;

	a_top_snp = (top_snp < 0 ? 1e30 : top_snp);

	num_snps = 0;
	GC_val = -1;
}

/*
 * cond_analysis default constructor
 */
cond_analysis::cond_analysis()
{
	a_out = "result";
	a_p_cutoff = 5e-8;
	a_collinear = 0.9;
	a_ld_window = 1e7;
	a_verbose = true;

	a_top_snp = -1;

	num_snps = 0;
	GC_val = -1;
}

/*
 * Initialise the conditional analysis by matching SNPs
 * and calculating frequencies. 
 * From `gcta::init_massoc`.
 */
void cond_analysis::init_conditional(phenotype *pheno, reference *ref)
{
	int i = 0, j = 0;
	size_t n, m;

	// First match the datasets
	match_gwas_phenotype(pheno, ref);
	// Re-caculate variance after matching to reference
	pheno->calc_variance();
	jma_Vp = pheno->pheno_variance;

	n = ref->to_include.size();
	m = ref->fam_ids_inc.size();

	msx_b.resize(n);
	nD.resize(n);

#pragma omp parallel for
	for (i = 0; i < n; i++) {
		eigenVector x;
		makex_eigenVector(i, x, ref);
		msx_b[i] = x.squaredNorm() / (double)m;
	}

	//msx = msx_b;
	//nD = ja_N_outcome;

	msx = 2.0 * ja_freq.array() * (1.0 - ja_freq.array());
	for (i = 0; i < n; i++) {
		nD[i] = (jma_Vp - msx[i] * ja_beta[i] * ja_beta[i]) / (msx[i] * ja_beta_se[i] * ja_beta_se[i]) + 1;
	}
}

/*
 * Matches reference SNPs to phenotype SNPs
 * @ret void
 */
void cond_analysis::match_gwas_phenotype(phenotype *pheno, reference *ref)
{
	size_t i = 0;
	map<string, int> snp_map_buffer(ref->snp_map);
	//map<size_t, size_t> bim_pheno; // Map of positions of SNPs in reference and phenotype, respectively
	map<string, int>::iterator iter;
	map<string, int> id_map;

	// Match the pre-matched SNPs from the exposure/outcome datasets to the reference dataset
	for (i = 0; i < pheno->matched_idx.size(); i++) {
		snp_map_buffer.erase(pheno->snp_name[pheno->matched_idx[i]]);
	}

	for (iter = snp_map_buffer.begin(); iter != snp_map_buffer.end(); iter++) {
		ref->snp_map.erase(iter->first);
	}

	ref->to_include.clear();
	for (iter = ref->snp_map.begin(); iter != ref->snp_map.end(); iter++) {
		ref->to_include.push_back(iter->second);
	}
	stable_sort(ref->to_include.begin(), ref->to_include.end());

	// Debug
	if (true) {
		string filename = "D:/Users/Jamie/Desktop/debug.txt";
		ofstream ofile(filename.c_str());
		size_t st;

		if (!ofile)
			ShowError("Cannot open file \"" + filename + "\" for writing.");

		// Header
		cout << "Ref SNP\tPheno SNP" << endl;
		ofile << "Ref SNP\tPheno SNP" << endl;
		for (st = 0; st < pheno->snp_name.size(); st++) {
			if (ref->snp_map[pheno->snp_name[st]] == 0)
				continue;
			cout << ref->bim_snp_name[ref->snp_map[pheno->snp_name[st]]] << "\t" << ref->snp_map[pheno->snp_name[st]] << "\t" << pheno->snp_name[st] << endl;
			ofile << ref->bim_snp_name[ref->snp_map[pheno->snp_name[st]]] << "\t" << ref->snp_map[pheno->snp_name[st]] << "\t" << pheno->snp_name[st] << endl;
		}
	}

	// Use the matched SNPs to find alleles and calculate mu
	vector<int> idx(ref->to_include.size());

	for (i = 0; i < pheno->snp_name.size(); i++)
		id_map.insert(pair<string, int>(pheno->snp_name[i], i));

	for (i = 0; i < ref->to_include.size(); i++) {
		iter = id_map.find(ref->bim_snp_name[ref->to_include[i]]);
		idx[i] = iter->second;
		ref->ref_A[ref->to_include[i]] = pheno->allele1[iter->second];

		if (!ref->mu.empty() && pheno->allele1[iter->second] == ref->bim_allele2[ref->to_include[i]])
			ref->mu[ref->to_include[i]] = 2.0 - ref->mu[ref->to_include[i]];
	}

	if (ref->to_include.empty()) {
		ShowError("Included list of SNPs is empty - could not match SNPs from phenotype file with reference SNPs.");
	}
	else {
		cout << "Total amount of SNPs matched from phenotype file with reference SNPs are: " << ref->to_include.size() << endl;
	}

	// Resize and get ready for the conditional analysis
	ja_freq.resize(ref->to_include.size());
	ja_beta.resize(ref->to_include.size());
	ja_beta_se.resize(ref->to_include.size());
	ja_pval.resize(ref->to_include.size());
	ja_N_outcome.resize(ref->to_include.size());

	for (i = 0; i < ref->to_include.size(); i++) {
		ja_freq[i] = pheno->freq[idx[i]];
		ja_beta[i] = pheno->beta[idx[i]];
		ja_beta_se[i] = pheno->se[idx[i]];
		ja_pval[i] = pheno->pval[idx[i]];
		ja_N_outcome[i] = pheno->n[idx[i]];
	}
}

void cond_analysis::stepwise_select(vector<int> &selected, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref)
{
	vector<double> p_temp;
	eigenVector2Vector(ja_pval, p_temp);
	size_t i = 0, prev_num = 0,
		m = min_element(p_temp.begin(), p_temp.end()) - p_temp.begin();

	if (p_temp[m] >= a_p_cutoff) {
		return;
	}
	selected.push_back(m);

	for (i = 0; i < ref->to_include.size(); i++) {
		if (i != m)
			remain.push_back(i);
	}

	if (a_p_cutoff > 1e-3) {
		ShowWarning("P value level is too low for stepwise model.", a_verbose);
	}

	while (!remain.empty()) {
		if (select_entry(selected, remain, bC, bC_se, pC, ref)) {
			selected_stay(selected, bC, bC_se, pC, ref);
		}
		else
			break;
		if (selected.size() % 5 == 0 && selected.size() > prev_num)
			cout << selected.size() << " associated SNPs have been selected." << endl;
		if (selected.size() > prev_num)
			prev_num = selected.size();
		if (selected.size() >= a_top_snp)
			break;
	}

	if (a_p_cutoff > 1e-3) {
		cout << "Performing backward elimination..." << endl;
		selected_stay(selected, bC, bC_se, pC, ref);
	}

	sw_snps = selected.size();
	cout << "Finally, " << selected.size() << " associated SNPs have been selected." << endl;
}

bool cond_analysis::insert_B_Z(const vector<int> &idx, int pos, reference *ref)
{
	bool get_ins_col = false, get_ins_row = false;
	int i = 0, j = 0, p,
		n = ref->fam_ids_inc.size(),
		m = ref->to_include.size();
	double d_temp = 0.0;
	vector<int> ix(idx);
	eigenSparseMat B_temp(B), B_N_temp(B_N);

	ix.push_back(pos);
	stable_sort(ix.begin(), ix.end());

	B.resize(ix.size(), ix.size());
	B_N.resize(ix.size(), ix.size());

	p = find(ix.begin(), ix.end(), pos) - ix.begin();
	eigenVector diagB(ix.size());
	eigenVector x_i(n), x_j(n);
	for (j = 0; j < ix.size(); j++) {
		B.startVec(j);
		B_N.startVec(j);
		B.insertBack(j, j) = msx_b[ix[j]];
		B_N.insertBack(j, j) = msx[ix[j]] * nD[ix[j]];

		diagB[j] = msx_b[ix[j]];
		if (pos == ix[j]) {
			get_ins_col = true;
			get_ins_row = true;
		}
		makex_eigenVector(ix[j], x_j, ref);
		
		for (i = j + 1; i < ix.size(); i++) {
			if (pos == ix[i])
				get_ins_row = true;

			if (pos == ix[j] || pos == ix[i]) {
				if (true // _jma_actual_geno
					|| (ref->bim_chr[ref->to_include[ix[i]]] == ref->bim_chr[ref->to_include[ix[j]]]
						&& abs(ref->bim_bp[ref->to_include[ix[i]]] - ref->bim_bp[ref->to_include[ix[j]]]) < a_ld_window)
					)
				{
					makex_eigenVector(ix[i], x_i, ref);
					d_temp = x_i.dot(x_j) / (double)n;
					B.insertBack(i, j) = d_temp;
					B_N.insertBack(i, j) = d_temp
											* min(nD[ix[i]], nD[ix[j]])
											* sqrt(msx[ix[i]] * msx[ix[j]] / (msx_b[ix[i]] * msx_b[ix[j]]));
				}
			}
			else {
				if (B_temp.coeff(i - get_ins_row, j - get_ins_col) != 0) {
					B.insertBack(i, j) = B_temp.coeff(i - get_ins_row, j - get_ins_col);
					B_N.insertBack(i, j) = B_N_temp.coeff(i - get_ins_row, j - get_ins_col);
				}
			}
		}
	}
	B.finalize();
	B_N.finalize();

	LDLT<eigenMatrix> ldlt_B(B);
	B_i = eigenMatrix::Identity(ix.size(), ix.size());
	ldlt_B.solveInPlace(B_i);
	if (ldlt_B.vectorD().minCoeff() < 0 || sqrt(ldlt_B.vectorD().maxCoeff() / ldlt_B.vectorD().minCoeff()) > 30
		|| (1 - eigenVector::Constant(ix.size(), 1).array() / (diagB.array() * B_i.diagonal().array())).maxCoeff() > a_collinear)
	{
		jma_snpnum_collinear++;
		B = B_temp;
		B_N = B_N_temp;
		return false;
	}

	LDLT<eigenMatrix> ldlt_B_N(B_N);
	B_N_i = eigenMatrix::Identity(ix.size(), ix.size());
	ldlt_B_N.solveInPlace(B_N_i);
	D_N.resize(ix.size());
	for (j = 0; j < ix.size(); j++) {
		D_N[j] = msx[ix[j]] * nD[ix[j]];
	}

	if (Z_N.cols() < 1)
		return true;

	eigenSparseMat Z_temp(Z), Z_N_temp(Z_N);
	Z.resize(ix.size(), m);
	Z_N.resize(ix.size(), m);

	for (j = 0; j < m; j++) {
		Z.startVec(j);
		Z_N.startVec(j);

		get_ins_row = false;
		makex_eigenVector(j, x_j, ref);
		for (i = 0; i < ix.size(); i++) {
			if (pos == ix[i]) {
				if (true // _jma_actual_geno
					|| (ix[i] != j 
						&& ref->bim_chr[ref->to_include[ix[i]]] == ref->bim_chr[ref->to_include[j]]
						&& abs(ref->bim_bp[ref->to_include[ix[i]]] - ref->bim_bp[ref->to_include[j]]) < a_ld_window)
					)
				{
					makex_eigenVector(ix[i], x_i, ref);
					d_temp = x_j.dot(x_i) / (double)n;
					Z.insertBack(i, j) = d_temp;
					Z_N.insertBack(i, j) = d_temp
											* min(nD[ix[i]], nD[j])
											* sqrt(msx[ix[i]] * msx[j] / (msx_b[ix[i]] * msx_b[j]));
				}
				get_ins_row = true;
			}
			else {
				if (Z_temp.coeff(i - get_ins_row, j) != 0) {
					Z.insertBack(i, j) = Z_temp.coeff(i - get_ins_row, j);
					Z_N.insertBack(i, j) = Z_N_temp.coeff(i - get_ins_row, j);
				}
			}
		}
	}

	Z.finalize();
	Z_N.finalize();
	return true;
}

void cond_analysis::erase_B_and_Z(const vector<int> &idx, int erase, reference *ref)
{
	bool get_ins_col = false, get_ins_row = false;
	int i = 0, j = 0,
		i_size = idx.size(),
		pos = find(idx.begin(), idx.end(), erase) - idx.begin(),
		m = ref->to_include.size();
	eigenSparseMat B_dense(B), B_N_dense(B_N);

	B.resize(i_size - 1, i_size - 1);
	B_N.resize(i_size - 1, i_size - 1);
	D_N.resize(i_size - 1);

	for (j = 0; j < i_size; j++) {
		if (erase == idx[j]) {
			get_ins_col = true;
			continue;
		}

		B.startVec(j - get_ins_col);
		B_N.startVec(j - get_ins_col);
		D_N[j - get_ins_col] = msx[idx[j]] * nD[idx[j]];
		get_ins_row = get_ins_col;

		for (i = j; i < i_size; i++) {
			if (erase == idx[i]) {
				get_ins_row = true;
				continue;
			}

			if (B_dense.coeff(i, j) != 0) {
				B.insertBack(i - get_ins_row, j - get_ins_col) = B_dense.coeff(i, j);
				B_N.insertBack(i - get_ins_row, j - get_ins_col) = B_N_dense.coeff(i, j);
			}
		}
	}
	B.finalize();
	B_N.finalize();

	if (Z_N.cols() < 1)
		return;

	LDLT<eigenMatrix> ldlt_B(B);
	B_i = eigenMatrix::Identity(i_size - 1, i_size - 1);
	ldlt_B.solveInPlace(B_i);
	LDLT<eigenMatrix> ldlt_B_N(B_N);
	B_N_i = eigenMatrix::Identity(i_size - 1, i_size - 1);
	ldlt_B_N.solveInPlace(B_N_i);

	eigenSparseMat Z_temp(Z), Z_N_temp(Z_N);
	Z.resize(i_size - 1, m);
	Z_N.resize(i_size - 1, m);
	for (j = 0; j < m; j++) {
		Z.startVec(j);
		Z_N.startVec(j);
		get_ins_row = false;
		for (i = 0; i < i_size; i++) {
			if (erase == idx[i]) {
				get_ins_row = true;
				continue;
			}

			if (Z_temp.coeff(i, j) != 0) {
				Z.insertBack(i - get_ins_row, j) = Z_temp.coeff(i, j);
				Z_N.insertBack(i - get_ins_row, j) = Z_N_temp.coeff(i, j);
			}
		}
	}
	Z.finalize();
	Z_N.finalize();
}

bool cond_analysis::select_entry(vector<int> &selected, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref)
{
	size_t m = 0;
	vector<double> pC_temp;

	massoc_conditional(selected, remain, bC, bC_se, pC, ref);

	eigenVector2Vector(pC, pC_temp);

	while (true) {
		m = min_element(pC_temp.begin(), pC_temp.end()) - pC_temp.begin();
		if (pC_temp[m] >= a_p_cutoff)
			return false;

		if (insert_B_Z(selected, remain[m], ref)) {
			selected.push_back(remain[m]);
			stable_sort(selected.begin(), selected.end());
			remain.erase(remain.begin() + m);
			return true;
		}

		pC_temp.erase(pC_temp.begin() + m);
		remain.erase(remain.begin() + m);
	}
}

void cond_analysis::selected_stay(vector<int> &select, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, reference *ref)
{
	if (B_N.cols() < 1) {
		if (!init_b(select, ref)) {
			ShowError("Stepwise Selection Error: There is a collinearity problem with the given list of SNPs.");
		}
	}

	vector<double> pJ_temp;
	while (!select.empty()) {
		massoc_joint(select, bJ, bJ_se, pJ, ref);
		eigenVector2Vector(pJ, pJ_temp);
		size_t m = max_element(pJ_temp.begin(), pJ_temp.end()) - pJ_temp.begin();
		if (pJ[m] > a_p_cutoff) {
			jma_snpnum_backward++;
			erase_B_and_Z(select, select[m], ref);
			select.erase(select.begin() + m);
		}
		else {
			break;
		}
	}
}

void cond_analysis::massoc_conditional(const vector<int> &selected, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref)
{
	int i = 0, j = 0, n = selected.size(), m = remain.size();
	double chisq = 0.0, B2 = 0.0;
	eigenVector b(n), se(n);

	if (B_N.cols() < 1) {
		if (!init_b(selected, ref)) {
			ShowError("Conditional Error: There is a collinearity problem with the SNPs given.\n");
		}
	}

	if (Z_N.cols() < 1) {
		init_z(selected, ref);
	}

	for (i = 0; i < n; i++) {
		b[i] = ja_beta[selected[i]];
		se[i] = ja_beta_se[selected[i]];
	}
	eigenVector bJ1 = B_N_i * D_N.asDiagonal() * b;
	//jma_Ve = massoc_calcu_Ve(selected, bJ1, b);

	eigenVector Z_Bi(n), Z_Bi_temp(n);
	bC = eigenVector::Zero(m);
	bC_se = eigenVector::Zero(m);
	pC = eigenVector::Constant(m, 2);
	for (i = 0; i < m; i++) {
		j = remain[i];
		B2 = msx[j] * nD[j];
		if (!isFloatEqual(B2, 0.0)) {
			Z_Bi = Z_N.col(j).transpose() * B_N_i;
			Z_Bi_temp = Z.col(j).transpose() * B_i;
			if (Z.col(j).dot(Z_Bi_temp) / msx_b[j] < a_collinear) {
				bC[i] = ja_beta[j] - Z_Bi.cwiseProduct(D_N).dot(b) / B2;
				bC_se[i] = (B2 - Z_N.col(j).dot(Z_Bi)) / (B2 * B2);
			}
		}
		bC_se[i] *= jma_Ve;
		if (bC_se[i] > 1.0e-7) {
			bC_se[i] = sqrt(bC_se[i]);
			chisq = bC[i] / bC_se[i];
			// If want to adjust for GC do it here
			// if (GC_val)
			pC[i] = pchisq(chisq * chisq, 1);
		}
	}
}

double cond_analysis::massoc_calcu_Ve(const vector<int> &selected, eigenVector &bJ, eigenVector &b)
{
	double Ve = 0.0, d_temp = 0.0;
	size_t n = bJ.size();
	vector<double> nD_temp(n);

	for (size_t k = 0; k < n; k++) {
		nD_temp[k] = nD[selected[k]];
		Ve += D_N[k] * bJ[k] * b[k];
	}

	d_temp = v_calc_median(nD_temp);
	if (d_temp - n < 1) {
		ShowError("DoF Error: Model is over-fitting due to lack of degree of freedom. Provide a more stringent P-value cutoff.");
	}
	Ve = ((d_temp - 1) * jma_Vp - Ve) / (d_temp - n);
	if (Ve <= 0.0) {
		ShowError("Residual Error: Residual variance is out of bounds meaning the model is over-fitting. Provide a more stringent P-value cutoff.");
	}
	return Ve;
}

void cond_analysis::makex_eigenVector(int j, eigenVector &x, reference *ref)
{
	size_t i = 0,
		n = ref->fam_ids_inc.size(),
		m = ref->to_include.size();
	x.resize(n);
#pragma omp parallel for
	for (i = 0; i < n; i++) {
		if (!ref->bed_snp_1[ref->to_include[j]][ref->fam_ids_inc[i]] || ref->bed_snp_2[ref->to_include[j]][ref->fam_ids_inc[i]])
		{
			if (ref->bim_allele1[ref->to_include[j]] == ref->ref_A[ref->to_include[j]])
				x[i] = (ref->bed_snp_1[ref->to_include[j]][ref->fam_ids_inc[i]] + ref->bed_snp_2[ref->to_include[j]][ref->fam_ids_inc[i]]);
			else
				x[i] = 2.0 - (ref->bed_snp_1[ref->to_include[j]][ref->fam_ids_inc[i]] + ref->bed_snp_2[ref->to_include[j]][ref->fam_ids_inc[i]]);
		}
		else {
			x[i] = ref->mu[ref->to_include[j]];
		}
		x[i] -= ref->mu[ref->to_include[j]];
	}
}

bool cond_analysis::init_b(const vector<int> &idx, reference *ref)
{
	int i = 0, j = 0, k = 0,
		n = ref->fam_ids_inc.size(),
		i_size = idx.size();
	double d_temp = 0.0;
	eigenVector diagB(i_size),
		x_i(n),
		x_j(n);

	B.resize(i_size, i_size);
	B_N.resize(i_size, i_size);
	D_N.resize(i_size);

	for (i = 0; i < i_size; i++) {
		D_N[i] = msx[idx[i]] * nD[idx[i]];
		B.startVec(i);
		B.insertBack(i, i) = msx_b[idx[i]];
		
		B_N.startVec(i);
		B_N.insertBack(i, i) = D_N[i];

		diagB[i] = msx_b[idx[i]];
		makex_eigenVector(idx[i], x_i, ref);

		for (j = i + 1; j < i_size; j++) {
			if (ref->bim_chr[ref->to_include[idx[i]]] == ref->bim_chr[ref->to_include[idx[j]]]
				&& abs(ref->bim_bp[ref->to_include[idx[i]]] - ref->bim_bp[ref->to_include[idx[j]]]) < a_ld_window)
			{
				makex_eigenVector(idx[j], x_j, ref);

				d_temp = x_i.dot(x_j) / (double)n;
				B.insertBack(j, i) = d_temp;
				B_N.insertBack(j, i) = d_temp 
									* min(nD[idx[i]], nD[idx[j]]) 
									* sqrt(msx[idx[i]] * msx[idx[j]] / (msx_b[idx[i]] * msx_b[idx[j]]));
			}
		}
	}

	B.finalize();
	B_N.finalize();

	LDLT<eigenMatrix> ldlt_B(B);
	if (ldlt_B.vectorD().minCoeff() < 0 || sqrt(ldlt_B.vectorD().maxCoeff() / ldlt_B.vectorD().minCoeff()) > 30)
		return false;
	B_i = eigenMatrix::Identity(i_size, i_size);
	ldlt_B.solveInPlace(B_i);
	if ((1 - eigenVector::Constant(i_size, 1).array() / (diagB.array() * B_i.diagonal().array())).maxCoeff() > a_collinear)
		return false;

	LDLT<eigenMatrix> ldlt_B_N(B_N);
	B_N_i = eigenMatrix::Identity(i_size, i_size);
	ldlt_B_N.solveInPlace(B_N_i);
	return true;
}

void cond_analysis::init_z(const vector<int> &idx, reference *ref)
{
	int i = 0, j = 0,
		n = ref->fam_ids_inc.size(),
		m = ref->to_include.size(),
		i_size = idx.size();
	double d_temp = 0.0;
	eigenVector x_i(n), x_j(n);

	Z.resize(i_size, m);
	Z_N.resize(i_size, m);

	for (j = 0; j < m; j++) {
		Z.startVec(j);
		Z_N.startVec(j);

		makex_eigenVector(j, x_j, ref);
		for (i = 0; i < i_size; i++) {
			if (idx[i] != j && ref->bim_chr[ref->to_include[idx[i]]] == ref->bim_chr[ref->to_include[j]]
				&& abs(ref->bim_bp[ref->to_include[idx[i]]] - ref->bim_bp[ref->to_include[j]]) < a_ld_window)
			{
				makex_eigenVector(idx[i], x_i, ref);

				d_temp = x_j.dot(x_i) / (double)n;
				Z.insertBack(i, j) = d_temp;
				Z_N.insertBack(i, j) = d_temp
										* min(nD[idx[i]], nD[j])
										* sqrt(msx[idx[i]] * msx[j] / (msx_b[idx[i]] * msx_b[j]));
			}
		}
	}

	Z.finalize();
	Z_N.finalize();
}

void cond_analysis::massoc_joint(const vector<int> &idx, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, reference *ref)
{
	int i = 0, n = idx.size();
	double chisq = 0.0;
	eigenVector b(n);
	for (i = 0; i < n; i++)
		b[i] = ja_beta[idx[i]];

	if (B_N.cols() < 1) {
		if (!init_b(idx, ref))
			ShowError("`massoc_joint`: There is a collinearity problem with the given list of SNPs.");
	}
	
	bJ.resize(n);
	bJ_se.resize(n);
	pJ.resize(n);
	bJ = B_N_i * D_N.asDiagonal() * b;
	bJ_se = B_N_i.diagonal();
	pJ = eigenVector::Ones(n);
	bJ_se *= jma_Ve;
	for (i = 0; i < n; i++) {
		if (bJ_se[i] > 1.0e-7) {
			bJ_se[i] = sqrt(bJ_se[i]);
			chisq = bJ[i] / bJ_se[i];
			if (GC_val > 0) {
				pJ[i] = pchisq(chisq * chisq / GC_val, 1.0);
			}
			else {
				pJ[i] = pchisq(chisq * chisq, 1.0);
			}
		}
		else {
			bJ[i] = 0.0;
			bJ_se[i] = 0.0;
		}
	}
}

vector<int> cond_analysis::read_snplist(string snplist, vector<int> &remain, reference *ref)
{
	size_t i = 0, n = ref->to_include.size();
	vector<string> givenSNPs;
	vector<int> pgiven;
	string temp;
	ifstream i_snplist(snplist.c_str());
	
	// Read from file
	givenSNPs.clear();
	if (!i_snplist) {
		ShowError("IO Error: Cannot read " + snplist + " to read SNP list.");
	}
	cout << "Reading SNPs upon which to condition from " + snplist + "." << endl;
	while (i_snplist >> temp) {
		givenSNPs.push_back(temp);
		getline(i_snplist, temp);
	}
	i_snplist.close();
	if (givenSNPs.empty()) {
		ShowError("No SNPs were read from the SNP list file - please check the format of this file.");
	}

	map<string, int> m_gSNPs;
	size_t snps_size = givenSNPs.size();
	pgiven.clear();
	remain.clear();
	for (i = 0; i < snps_size; i++) {
		m_gSNPs.insert(pair<string, int>(givenSNPs[i], static_cast<int>(i)));
	}
	for (i = 0; i < n; i++) {
		if (m_gSNPs.find(ref->bim_snp_name[ref->to_include[i]]) != m_gSNPs.end()) {
			pgiven.push_back(i);
		}
		else {
			remain.push_back(i);
		}
	}
	if (pgiven.size() > 0) {
		cout << pgiven.size() << " conditional SNP(s) were matched to the reference dataset." << endl;
	}
	else {
		ShowError("None of the SNPs from the SNP list could be matched. Please double check the datasets.");
	}
	return pgiven;
}

/*
 * Run step-wise selection to find independent association signals/SNP
 */
void cond_analysis::massoc(reference *ref, string snplist)
{
	vector<int> selected, remain, pgiven;
	eigenVector bC, bC_se, pC;

	if (a_top_snp < 0.0)
		a_top_snp = 1e30;

	cout << "Performing stepwise model selection on " << ref->to_include.size() << " SNPs to select association signals (p cutoff = " << a_p_cutoff << ", assuming complete LE between SNPs more than " << a_ld_window / 1e6 << " Mb away)." << endl;
	stepwise_select(selected, remain, bC, bC_se, pC, ref);

	if (selected.empty()) {
		ShowError("Conditional Error: No SNPs have been selected by the step-wise selection algorithm.");
	}
	else if (selected.size() >= ref->fam_ids_inc.size()) {
		ShowError("Conditional Error: Too many SNPs. The number of SNPs should not be larger than the sample size.");
	}

	cout << "(" << jma_snpnum_backward << " SNPs eliminated by backward selection.)" << endl;

	// Perform the joint analysis
	eigenVector bJ, bJ_se, pJ;
	massoc_joint(selected, bJ, bJ_se, pJ, ref);

	eigenMatrix rval(selected.size(), selected.size());
	LD_rval(selected, rval);
	
	sanitise_output(selected, bC, bC_se, pC, rval, CO_COND, ref);
	sanitise_output(selected, bJ, bJ_se, pJ, rval, CO_JOINT, ref);
}

void cond_analysis::LD_rval(const vector<int> &idx, eigenMatrix &rval)
{
	int i = 0, j = 0,
		i_size = idx.size();
	eigenVector sd(i_size);

	for (i = 0; i < i_size; i++)
		sd[i] = sqrt(msx_b[idx[i]]);

	for (j = 0; j < i_size; j++) {
		rval(j, j) = 1.0;
		for (i = j + 1; i < i_size; i++) {
			rval(i, j) = rval(j, i) = B.coeff(i, j) / sd[i] / sd[j];
		}
	}
}

void cond_analysis::sanitise_output(vector<int> &selected, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, eigenMatrix &rval, enum cond_type ctype, reference *ref)
{
	string filename = a_out + (ctype == CO_COND ? ".cma.cojo" : "jma.cojo");
	ofstream ofile(filename.c_str());
	size_t i = 0, j = 0;
	
	if (!ofile)
		ShowError("Cannot open file \"" + filename + "\" for writing.");

	// Header
	ofile << "Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp\tn\tfreq_geno\tbC\tbC_se\tpC";
	if (ctype == CO_JOINT)
		ofile << "\tLD_r";
	ofile << endl;

	snps_cond.resize(selected.size());
	for (i = 0; i < selected.size(); i++) {
		j = selected[i];
		if (ctype == CO_COND)
			snps_cond.push_back(ref->bim_snp_name[ref->to_include[j]]);
		ofile << ref->bim_chr[ref->to_include[j]] << "\t" << ref->bim_snp_name[ref->to_include[j]] << "\t" << ref->bim_bp[ref->to_include[j]] << "\t";
		ofile << ref->ref_A[ref->to_include[j]] << "\t" << ja_freq[j] << "\t" << ja_beta[j] << "\t" << ja_beta_se[j] << "\t";
		ofile << ja_pval[j] << "\t" << nD[j] << "\t" << 0.5 * ref->mu[ref->to_include[j]];

		if (ctype == CO_COND) {
			if (pJ[i] > 1.5)
				ofile << "NA\tNA\tNA" << endl;
			else
				ofile << bJ[i] << "\t" << bJ_se[i] << "\t" << pJ[i] << endl;
		}
		else {
			ofile << bJ[i] << "\t" << bJ_se[i] << "\t" << pJ[i] << endl;
			if (i == selected.size() - 1)
				ofile << 0 << endl;
			else
				ofile << rval(i, i + 1) << endl;
		}
	}
	ofile.close();

	if (ctype == CO_COND) {
		eigenVector2Vector(bJ, b_cond);
		eigenVector2Vector(bJ_se, se_cond);
		eigenVector2Vector(ja_freq, maf_cond); // Coloc requires MAF - this is converted elsewhere
		eigenVector2Vector(pJ, p_cond);
		eigenVector2Vector(nD, n_cond);
	}
}
