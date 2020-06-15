#include "conditional.h"

/*
 * cond_analysis constructor
 */
cond_analysis::cond_analysis(double p_cutoff, double collinear, int ld_window, string out, bool verbose, int top_snp)
{
	a_out = out;
	a_p_cutoff = p_cutoff;
	a_collinear = collinear;
	a_ld_window = ld_window;
	a_verbose = verbose;

	a_top_snp = (top_snp < 0 ? 1e30 : top_snp);

	num_snps = 0;
}

/*
 * cond_analysis default constructor
 */
cond_analysis::cond_analysis()
{
	a_out = "result";
	a_p_cutoff = 1e-5;
	a_collinear = 0.9;
	a_ld_window = 10000 * 1000;
	a_verbose = true;

	a_top_snp = -1;

	num_snps = 0;
}

/*
 * Initialise the conditional analysis by matching SNPs
 * and calculating frequencies. 
 */
void cond_analysis::init_conditional(phenotype *pheno, reference *ref)
{
	int i = 0, j = 0;
	size_t inc_size = ref->to_include.size(),
		keep_size = ref->fam_ids_inc.size();

	// First match the datasets
	match_gwas_phenotype(pheno, ref);

	// Now calculate the allele frequencies
	ref->calculate_allele_freq();
	
	msx_b.resize(keep_size);
	nD.resize(keep_size);

#pragma omp parallel for
	for (i = 0; i < keep_size; i++) {
		eigenVector x;
		x.resize(keep_size);
		for (j = 0; j < inc_size; j++) {
			if (!ref->bed_snp_1[ref->to_include[j]][ref->fam_ids_inc[i]]
				|| ref->bed_snp_2[ref->to_include[j]][ref->fam_ids_inc[i]]) 
			{
				int temp = ref->bed_snp_1[ref->to_include[j]][ref->fam_ids_inc[i]] + ref->bed_snp_2[ref->to_include[j]][ref->fam_ids_inc[i]];
				if (ref->bim_allele1[ref->to_include[j]] == ref->ref_A[ref->to_include[j]])
					x[i] = temp;
				else
					x[i] = 2.0 - temp;
			}
			else {
				x[i] = ref->mu[ref->to_include[j]];
			}
			x[i] -= ref->mu[ref->to_include[j]];
		}

		msx_b[i] = x.squaredNorm() / (double)keep_size;
	}

	if (false) { // _jma_actual_geno flag
		;
	}
	else {
		msx = 2.0 * ja_freq.array() * (1.0 - ja_freq.array());
		for (i = 0; i < inc_size; i++) {
			nD[i] = (jma_Vp - msx[i] * ja_beta[i] * ja_beta[i]) / (msx[i] * ja_beta_se[i] * ja_beta_se[i]) + 1;
		}
	}
}

/*
 * Matches GWAS SNPs to phenotype SNPs
 * @ret void
 */
void cond_analysis::match_gwas_phenotype(phenotype *pheno, reference *ref)
{
	int i = 0;
	map<string, int> snp_map_buffer(ref->snp_map);
	map<string, int>::iterator iter;
	map<string, int> id_map;

	// First check if the SNPs can be found in the reference GWAS
	// Now match the GWAS data to the genotype data
	// We will also handle "bad" SNPs here - i.e. those that do not match
	// alleles with the genotype data.
	for (i = 0; i < pheno->snp_name.size(); i++)
		snp_map_buffer.erase(pheno->snp_name[i]);

	for (iter = snp_map_buffer.begin(); iter != snp_map_buffer.end(); iter++) {
		ref->snp_map.erase(iter->first);
	}

	ref->to_include.clear();
	for (iter = ref->snp_map.begin(); iter != ref->snp_map.end(); iter++)
		ref->to_include.push_back(iter->second);
	stable_sort(ref->to_include.begin(), ref->to_include.end());

	// Use the matched SNPs to find alleles and calculate mu
	vector<int> idx(ref->to_include.size());

	for (i = 0; i < pheno->snp_name.size(); i++)
		id_map.insert(pair<string, int>(pheno->snp_name[i], i));

	for (i = 0; i < ref->to_include.size(); i++) {
		iter = id_map.find(ref->bim_snp_name[ref->to_include[i]]);
		idx[i] = iter->second;
		ref->ref_A[ref->to_include[i]] = pheno->allele1[iter->second];

		if (!ref->mu.empty() && pheno->allele1[iter->second] == pheno->allele2[ref->to_include[i]])
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

	jma_Ve = pheno->pheno_variance;
}

void cond_analysis::stepwise_select(vector<int> &selected, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref)
{
	int i = 0, buf = 0, prev_num = 0, m;
	vector<double> p_temp;
	eigenVector2Vector(ja_pval, p_temp);

	m = min_element(p_temp.begin(), p_temp.end()) - p_temp.begin();
	if (p_temp[m] >= a_p_cutoff) {
		return;
	}

	selected.push_back(m);
	for (i = 0; i < ref->to_include.size(); i++) {
		if (i != m)
			remain.push_back(i);
	}

	while (!remain.empty()) {
		if (select_entry(selected, remain, bC, bC_se, pC, ref)) {
			// Do me
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
		//selected_stay
	}

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
				if (ref->bim_chr[ref->to_include[ix[i]]] == ref->bim_chr[ref->to_include[ix[j]]]
					&& abs(ref->bim_bp[ref->to_include[ix[i]]] - ref->bim_bp[ref->to_include[ix[j]]]) < a_ld_window)
				{
					makex_eigenVector(ix[i], x_i, ref);
					d_temp = x_i.dot(x_j) / (double)n;
					B.insertBack(i, j) = d_temp;
					B_N.insertBack(i, j) = d_temp
											* min(nD[idx[i]], nD[idx[j]])
											* sqrt(msx[idx[i]] * msx[idx[j]] / (msx_b[idx[i]] * msx_b[idx[j]]));
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
				if (ix[i] != j 
					&& ref->bim_chr[ref->to_include[ix[i]]] == ref->bim_chr[ref->to_include[j]]
					&& abs(ref->bim_bp[ref->to_include[ix[i]]] - ref->bim_bp[ref->to_include[j]]) < a_ld_window)
				{
					makex_eigenVector(ix[i], x_i, ref);
					d_temp = x_j.dot(x_i) / (double)n;
					Z.insertBack(i, j) = d_temp;
					Z_N.insertBack(i, j) = d_temp
											* min(nD[idx[i]], nD[idx[j]])
											* sqrt(msx[idx[i]] * msx[idx[j]] / (msx_b[idx[i]] * msx_b[idx[j]]));
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

bool cond_analysis::select_entry(vector<int> &selected, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref)
{
	int i = 0, m = 0;
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

void cond_analysis::massoc_conditional(const vector<int> &selected, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref)
{
	int i = 0, j = 0, n = selected.size(), m = remain.size();
	double chisq = 0.0, B2 = 0.0;
	eigenVector b(n), se(n);

	if (B_N.cols() < 1) {
		if (init_b(selected, ref) == false) {
			ShowError("There is a collinearity problem with the SNPs given.\n");
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

void cond_analysis::makex_eigenVector(int j, eigenVector &x, reference *ref)
{
	int i = 0,
		n = ref->fam_ids_inc.size(),
		m = ref->to_include.size();

#pragma omp parallel for
	for (i = 0; i < n; i++) {
		if (!ref->bed_snp_1[ref->to_include[j]][ref->fam_ids_inc[i]] || ref->bed_snp_2[ref->to_include[j]][ref->fam_ids_inc[i]])
		{
			if (ref->bim_allele1[ref->to_include[j]] == ref->ref_A[ref->to_include[j]])
				x[i] = (ref->bed_snp_1[ref->to_include[j]][ref->fam_ids_inc[i]] + ref->bed_snp_2[ref->to_include[j]][ref->fam_ids_inc[i]]);
			else
				x[i] = 2.0 - (ref->bed_snp_1[ref->to_include[j]][ref->fam_ids_inc[i]] + ref->bed_snp_2[ref->to_include[j]][ref->fam_ids_inc[i]]);

			x[i] -= ref->mu[ref->to_include[j]];
		}
		else {
			x[i] = 0;
		}
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
										* min(nD[idx[i]], nD[idx[j]])
										* sqrt(msx[idx[i]] * msx[idx[j]] / (msx_b[idx[i]] * msx_b[idx[j]]));
			}
		}
	}

	Z.finalize();
	Z_N.finalize();
}

void cond_analysis::massoc(cond_analysis *p_analysis)
{
	
}