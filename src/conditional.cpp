#include "conditional.h"

/*
 * cond_analysis constructor
 */
cond_analysis::cond_analysis(double p_cutoff, double collinear, double ld_window, string out, double top_snp, double freq_thres, string name, bool cond_ssize)
{
	cname = name;
	a_out = out;
	a_p_cutoff = p_cutoff;
	a_collinear = collinear;
	a_ld_window = ld_window;
	a_freq_threshold = freq_thres;

	a_top_snp = top_snp;

	num_snps = 0;
	cond_ssize = cond_ssize;
}

/*
 * cond_analysis default constructor
 */
cond_analysis::cond_analysis()
{
	cname = "conditional-default";
	a_out = "result";
	a_p_cutoff = 5e-8;
	a_collinear = 0.9;
	a_ld_window = 1e7;
	a_freq_threshold = 0.2;

	a_top_snp = 1e10;

	num_snps = 0;
	cond_ssize = false;
}

/*
 * Initialise the conditional analysis by matching SNPs
 * and calculating frequencies. 
 * From `gcta::init_massoc`.
 */
void cond_analysis::init_conditional(phenotype *pheno, reference *ref)
{
	size_t j = 0;
	size_t n, m;

	// First match the datasets
	match_gwas_phenotype(pheno, ref);
	jma_Vp = jma_Ve = pheno->get_variance();

	n = to_include.size();
	m = fam_ids_inc.size();

	msx_b.resize(n);
	nD.resize(n);

#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		eigenVector x;
		makex_eigenVector(i, x, true, ref);
		msx_b[i] = x.squaredNorm() / (double)m;
	}

	msx = 2.0 * ja_freq.array() * (1.0 - ja_freq.array());

	for (int i = 0; i < n; i++) {
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
	map<string, size_t>::iterator iter;
	map<string, size_t> id_map;
	vector<size_t> idx, pheno_idx, bad_idx;
	vector<string> snps, pheno_snps = pheno->get_snp_names();
	unsigned int unmatched = 0;

	to_include.clear();
	ref->includes_clear();

	// Match GWAS data to reference data and initialise the inclusion list
	//ref->update_inclusion(pheno->matched_idx, pheno->snp_name);
	for (i = 0; i < pheno_snps.size(); i++) {
		if ((iter = ref->snp_map.find(pheno_snps[i])) == ref->snp_map.end()
			|| (pheno->allele1[i] != ref->bim_allele1[iter->second] && (pheno->allele1[i] != ref->bim_allele2[iter->second]))
			)
		{
			if (iter != ref->snp_map.end()) {
				unmatched++;
				bad_idx.push_back(iter->second);
			}
			continue;
		}
		pheno_idx.push_back(i);
		id_map.insert(pair<string, size_t>(pheno_snps[i], i));
		snps.push_back(iter->first);
		idx.push_back(iter->second);
	}
	ref->update_inclusion(idx, snps);
	snps.clear();
	idx.clear();

	// Allele matching and swapping
	mu = ref->mu; // Copy across as this will be morphed below
	for (i = 0; i < ref->to_include.size(); i++) {
		bool flip_allele = false;
		iter = id_map.find(ref->bim_snp_name[ref->to_include[i]]);
		if (iter == id_map.end())
			continue;

		ref->ref_A[ref->to_include[i]] = pheno->allele1[iter->second];

		if (!ref->mu.empty() && pheno->allele1[iter->second] == ref->bim_allele2[ref->to_include[i]]) {
			mu[ref->to_include[i]] = 2.0 - ref->mu[ref->to_include[i]];
			flip_allele = true;
		}
		else {
			mu[ref->to_include[i]] = ref->mu[ref->to_include[i]];
		}

		double cur_freq = mu[ref->to_include[i]] / 2.0;
		double freq_diff = abs(cur_freq - pheno->freq[iter->second]);
		if (
			//(ref->bim_allele1[ref->to_include[i]] == pheno->allele1[iter->second]
			//	|| ref->bim_allele2[ref->to_include[i]] == pheno->allele1[iter->second]
			//) &&
			freq_diff < a_freq_threshold) 
		{
			snps.push_back(iter->first);
			idx.push_back(iter->second);
		}
		else {
			unmatched++;
			bad_idx.push_back(ref->to_include[i]);
		}
	}

	if (unmatched) {
		spdlog::info("[{}] There were {} SNPs that had a large difference in the allele frequency to that of the reference sample.", pheno->get_phenoname(), unmatched);
		string filename = a_out + "." + pheno->get_phenoname() + ".badfreq";
		ofstream file(filename.c_str());

		file << "SNP\tAllele1\tAllele2\tRefA\tfreq_ref\tfreq_pheno" << endl;
		for (i = 0; i < bad_idx.size(); i++) {
			double freq = (iter = id_map.find(ref->bim_snp_name[bad_idx[i]])) == id_map.end() ? -1.0 : pheno->freq[iter->second];
			file << ref->bim_snp_name[bad_idx[i]] << "\t" << ref->bim_allele1[bad_idx[i]] << "\t" << ref->bim_allele2[bad_idx[i]] << "\t" << ref->ref_A[bad_idx[i]] << "\t" << mu[bad_idx[i]] / 2.0 << "\t" << freq << endl;
		}
		file.close();
	}

	ref->update_inclusion(idx, snps);
	to_include = ref->to_include;
	to_include_bim = ref->to_include_bim;
	fam_ids_inc = ref->fam_ids_inc;

	if (to_include.empty()) {
		spdlog::critical("Included list of SNPs is empty - could not match SNPs from phenotype file with reference SNPs.");
		return;
	}
	else {
		spdlog::info("[{}] Total amount of SNPs matched from phenotype file with reference SNPs are: {}.", cname, to_include.size());
	}

	ctype = pheno->get_coloc_type();

	// Resize and get ready for the conditional analysis
	ja_snp_name.resize(to_include.size());
	ja_freq.resize(to_include.size());
	ja_beta.resize(to_include.size());
	ja_beta_se.resize(to_include.size());
	ja_pval.resize(to_include.size());
	ja_chisq.resize(to_include.size());
	ja_N_outcome.resize(to_include.size()); // Unused?
	ja_n_cases.resize(to_include.size());
	nsample.resize(to_include.size());
	if (ctype == coloc_type::COLOC_CC)
		ncases.resize(to_include.size());

	for (i = 0; i < to_include.size(); i++) {
		ja_snp_name[i] = pheno->snp_name[idx[i]];
		ja_freq[i] = pheno->freq[idx[i]];
		ja_beta[i] = pheno->beta[idx[i]];
		ja_beta_se[i] = pheno->se[idx[i]];
		//ja_pval[i] = pheno->pval[idx[i]];
		ja_chisq[i] = (ja_beta[i] / ja_beta_se[i]) * (ja_beta[i] / ja_beta_se[i]);
		ja_pval[i] = pchisq(ja_chisq[i], 1.0);
		ja_N_outcome[i] = pheno->n[idx[i]];
		nsample[i] = pheno->n[idx[i]];
		if (ctype == coloc_type::COLOC_CC) {
			ja_n_cases[i] = pheno->n_case[idx[i]];
			ncases[i] = pheno->n_case[idx[i]];
		}
	}

	string filename = a_out + "." + pheno->get_phenoname() + ".included";
	ofstream file(filename.c_str());

	file << "SNP\tChisq\tB\tSE\tPval\tFreq" << endl;
	for (size_t j = 0; j < ja_snp_name.size(); j++) {
		file << ja_snp_name[j] << "\t" << ja_chisq[j] << "\t" << ja_beta[j] << "\t" << ja_beta_se[j] << "\t" << ja_pval[j] << "\t" << ja_freq[j] << endl;
	}
	file.close();
}

void cond_analysis::stepwise_select(vector<size_t> &selected, vector<size_t> &remain, conditional_dat *cdat, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref)
{
	vector<double> p_temp, chisq;
	eigenVector2Vector(ja_pval, p_temp);
	eigenVector2Vector(ja_chisq, chisq);
	size_t i = 0, prev_num = 0,
		//m = min_element(p_temp.begin(), p_temp.end()) - p_temp.begin();
		m = max_element(chisq.begin(), chisq.end()) - chisq.begin();;

	spdlog::info("[{}] Selected SNP {} with chisq {:.2f} and pval {:.2e}.", cname, ja_snp_name[m], ja_chisq[m], ja_pval[m]);
	if (ja_pval[m] >= a_p_cutoff) {
		spdlog::info("[{}] SNP did not meet threshold.", cname);
		return;
	}
	selected.push_back(m);

	for (i = 0; i < to_include.size(); i++) {
		if (i != m)
			remain.push_back(i);
	}

	if (a_p_cutoff > 1e-3) {
		spdlog::warn("P value level is too low for stepwise model.");
	}

	while (!remain.empty()) {
		if (select_entry(selected, remain, cdat, bC, bC_se, pC, ref)) {
			selected_stay(selected, cdat, bC, bC_se, pC, ref);
		}
		else
			break;
		if (selected.size() % 5 == 0 && selected.size() > prev_num)
			spdlog::info("[{}] {} associated SNPs have been selected.", cname, selected.size());
		if (selected.size() > prev_num)
			prev_num = selected.size();
		if (selected.size() >= a_top_snp)
			break;
	}

	if (a_p_cutoff > 1e-3) {
		spdlog::info("Performing backward elimination...");
		selected_stay(selected, cdat, bC, bC_se, pC, ref);
	}

	spdlog::info("[{}] Finally, {} associated SNPs have been selected.", cname, selected.size());
}

bool cond_analysis::insert_B_Z(const vector<size_t> &idx, size_t pos, conditional_dat *cdat, reference *ref)
{
	bool get_ins_col = false, get_ins_row = false;
	size_t i = 0, j = 0, p,
		n = fam_ids_inc.size(),
		m = to_include.size();
	double d_temp = 0.0;
	vector<size_t> ix(idx);
	eigenSparseMat B_temp(cdat->B), B_N_temp(cdat->B_N);

	ix.push_back(pos);
	stable_sort(ix.begin(), ix.end());

	cdat->B.resize(ix.size(), ix.size());
	cdat->B_N.resize(ix.size(), ix.size());

	p = find(ix.begin(), ix.end(), pos) - ix.begin();
	eigenVector diagB(ix.size());
	eigenVector x_i(n), x_j(n);
	for (j = 0; j < ix.size(); j++) {
		cdat->B.startVec(j);
		cdat->B_N.startVec(j);
		cdat->B.insertBack(j, j) = msx_b[ix[j]];
		cdat->B_N.insertBack(j, j) = msx[ix[j]] * nD[ix[j]];

		diagB[j] = msx_b[ix[j]];
		if (pos == ix[j]) {
			get_ins_col = true;
		}
		get_ins_row = get_ins_col;
		makex_eigenVector(ix[j], x_j, false, ref);
		
		for (i = j + 1; i < ix.size(); i++) {
			if (pos == ix[i])
				get_ins_row = true;

			if (pos == ix[j] || pos == ix[i]) {
				if ((ref->bim_chr[to_include[ix[i]]] == ref->bim_chr[to_include[ix[j]]]
						&& abs(ref->bim_bp[to_include[ix[i]]] - ref->bim_bp[to_include[ix[j]]]) < a_ld_window)
					)
				{
					makex_eigenVector(ix[i], x_i, false, ref);
					d_temp = x_i.dot(x_j) / (double)n;
					cdat->B.insertBack(i, j) = d_temp;
					cdat->B_N.insertBack(i, j) = d_temp
											* min(nD[ix[i]], nD[ix[j]])
											* sqrt(msx[ix[i]] * msx[ix[j]] / (msx_b[ix[i]] * msx_b[ix[j]]));
				}
			}
			else {
				size_t ins_row_val = get_ins_row ? 1 : 0,
					ins_col_val = get_ins_col ? 1 : 0;
				if (B_temp.coeff(i - ins_row_val, j - ins_col_val) != 0) {
					cdat->B.insertBack(i, j) = B_temp.coeff(i - ins_row_val, j - ins_col_val);
					cdat->B_N.insertBack(i, j) = B_N_temp.coeff(i - ins_row_val, j - ins_col_val);
				}
			}
		}
	}
	cdat->B.finalize();
	cdat->B_N.finalize();

	SimplicialLDLT<eigenSparseMat> ldlt_B(cdat->B);
	cdat->B_i.resize(ix.size(), ix.size());
	cdat->B_i.setIdentity();
	cdat->B_i = ldlt_B.solve(cdat->B_i).eval();
	if (ldlt_B.vectorD().minCoeff() < 0 || sqrt(ldlt_B.vectorD().maxCoeff() / ldlt_B.vectorD().minCoeff()) > 30
		|| (1 - eigenVector::Constant(ix.size(), 1).array() / (diagB.array() * cdat->B_i.diagonal().array())).maxCoeff() > a_collinear)
	{
		jma_snpnum_collinear++;
		cdat->B = B_temp;
		cdat->B_N = B_N_temp;
		return false;
	}

	SimplicialLDLT<eigenSparseMat> ldlt_B_N(cdat->B_N);
	cdat->B_N_i.resize(ix.size(), ix.size());
	cdat->B_N_i.setIdentity();
	cdat->B_N_i = ldlt_B_N.solve(cdat->B_N_i).eval();
	cdat->D_N.resize(ix.size());
	for (j = 0; j < ix.size(); j++) {
		cdat->D_N[j] = msx[ix[j]] * nD[ix[j]];
	}

	if (cdat->Z_N.cols() < 1)
		return true;

	eigenSparseMat Z_temp(cdat->Z), Z_N_temp(cdat->Z_N);
	cdat->Z.resize(ix.size(), m);
	cdat->Z_N.resize(ix.size(), m);

	for (j = 0; j < m; j++) {
		cdat->Z.startVec(j);
		cdat->Z_N.startVec(j);

		get_ins_row = false;
		makex_eigenVector(j, x_j, false, ref);
		for (i = 0; i < ix.size(); i++) {
			if (pos == ix[i]) {
				if ((ix[i] != j 
						&& ref->bim_chr[to_include[ix[i]]] == ref->bim_chr[to_include[j]]
						&& abs(ref->bim_bp[to_include[ix[i]]] - ref->bim_bp[to_include[j]]) < a_ld_window)
					)
				{
					makex_eigenVector(ix[i], x_i, false, ref);
					d_temp = x_j.dot(x_i) / (double)n;
					cdat->Z.insertBack(i, j) = d_temp;
					cdat->Z_N.insertBack(i, j) = d_temp
											* min(nD[ix[i]], nD[j])
											* sqrt(msx[ix[i]] * msx[j] / (msx_b[ix[i]] * msx_b[j]));
				}
				get_ins_row = true;
			}
			else {
				size_t ins_row_val = get_ins_row ? 1 : 0;
				if (Z_temp.coeff(i - ins_row_val, j) != 0) {
					cdat->Z.insertBack(i, j) = Z_temp.coeff(i - ins_row_val, j);
					cdat->Z_N.insertBack(i, j) = Z_N_temp.coeff(i - ins_row_val, j);
				}
			}
		}
	}

	cdat->Z.finalize();
	cdat->Z_N.finalize();
	return true;
}

void cond_analysis::erase_B_and_Z(const vector<size_t> &idx, size_t erase, conditional_dat *cdat)
{
	bool get_ins_col = false, get_ins_row = false;
	size_t i = 0, j = 0,
		i_size = idx.size(),
		pos = find(idx.begin(), idx.end(), erase) - idx.begin(),
		m = to_include.size();
	eigenSparseMat B_dense(cdat->B), B_N_dense(cdat->B_N);

	cdat->B.resize(i_size - 1, i_size - 1);
	cdat->B_N.resize(i_size - 1, i_size - 1);
	cdat->D_N.resize(i_size - 1);

	for (j = 0; j < i_size; j++) {
		if (erase == idx[j]) {
			get_ins_col = true;
			continue;
		}

		cdat->B.startVec(j - get_ins_col);
		cdat->B_N.startVec(j - get_ins_col);
		cdat->D_N[j - (get_ins_col ? 1 : 0)] = msx[idx[j]] * nD[idx[j]];
		get_ins_row = get_ins_col;

		for (i = j; i < i_size; i++) {
			if (erase == idx[i]) {
				get_ins_row = true;
				continue;
			}

			if (B_dense.coeff(i, j) != 0) {
				size_t ins_row_val = get_ins_row ? 1 : 0,
					ins_col_val = get_ins_col ? 1 : 0;
				cdat->B.insertBack(i - ins_row_val, j - ins_col_val) = B_dense.coeff(i, j);
				cdat->B_N.insertBack(i - ins_row_val, j - ins_col_val) = B_N_dense.coeff(i, j);
			}
		}
	}
	cdat->B.finalize();
	cdat->B_N.finalize();

	if (cdat->Z_N.cols() < 1)
		return;

	SimplicialLDLT<eigenSparseMat> ldlt_B(cdat->B);
	cdat->B_i.resize(i_size - 1, i_size - 1);
	cdat->B_i.setIdentity();

	cdat->B_i = ldlt_B.solve(cdat->B_i).eval();
	SimplicialLDLT<eigenSparseMat> ldlt_B_N(cdat->B_N);
	cdat->B_N_i.resize(i_size - 1, i_size - 1);
	cdat->B_N_i.setIdentity();
	cdat->B_N_i = ldlt_B_N.solve(cdat->B_N_i).eval();

	eigenSparseMat Z_temp(cdat->Z), Z_N_temp(cdat->Z_N);
	cdat->Z.resize(i_size - 1, m);
	cdat->Z_N.resize(i_size - 1, m);
	for (j = 0; j < m; j++) {
		cdat->Z.startVec(j);
		cdat->Z_N.startVec(j);
		get_ins_row = false;
		for (i = 0; i < i_size; i++) {
			if (erase == idx[i]) {
				get_ins_row = true;
				continue;
			}

			if (Z_temp.coeff(i, j) != 0) {
				size_t ins_row_val = get_ins_row ? 1 : 0;
				cdat->Z.insertBack(i - ins_row_val, j) = Z_temp.coeff(i, j);
				cdat->Z_N.insertBack(i - ins_row_val, j) = Z_N_temp.coeff(i, j);
			}
		}
	}
	cdat->Z.finalize();
	cdat->Z_N.finalize();
}

bool cond_analysis::select_entry(vector<size_t> &selected, vector<size_t> &remain, conditional_dat *cdat, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref)
{
	size_t m = 0;
	vector<double> pC_temp;

	massoc_conditional(selected, remain, cdat, bC, bC_se, pC, ref);

	eigenVector2Vector(pC, pC_temp);

	while (true) {
		m = min_element(pC_temp.begin(), pC_temp.end()) - pC_temp.begin();
		spdlog::info("[{}] Selected entry SNP {} with cpval {:.2e}.", cname, ja_snp_name[m], pC_temp[m]);
		if (pC_temp[m] >= a_p_cutoff) {
			spdlog::info("[{}] {} does not meet threshold", cname, ja_snp_name[m]);
			return false;
		}

		if (insert_B_Z(selected, remain[m], cdat, ref)) {
			selected.push_back(remain[m]);
			stable_sort(selected.begin(), selected.end());
			remain.erase(remain.begin() + m);
			return true;
		}

		pC_temp.erase(pC_temp.begin() + m);
		remain.erase(remain.begin() + m);
	}
}

void cond_analysis::selected_stay(vector<size_t> &select, conditional_dat *cdat, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, reference *ref)
{
	if (cdat->B_N.cols() < 1) {
		if (!init_b(select, cdat, ref)) {
			spdlog::critical("There is a collinearity problem with the given list of SNPs.");
			return;
		}
	}

	vector<double> pJ_temp;
	while (!select.empty()) {
		massoc_joint(select, cdat, bJ, bJ_se, pJ, ref);
		eigenVector2Vector(pJ, pJ_temp);
		size_t m = max_element(pJ_temp.begin(), pJ_temp.end()) - pJ_temp.begin();
		if (pJ[m] > a_p_cutoff) {
			jma_snpnum_backward++;
			erase_B_and_Z(select, select[m], cdat);
			select.erase(select.begin() + m);
			spdlog::info("[{}] Erasing SNP {}.", cname, ja_snp_name[m]);
		}
		else {
			break;
		}
	}
}

void cond_analysis::massoc_conditional(const vector<size_t> &selected, vector<size_t> &remain, conditional_dat *cdat, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, reference *ref)
{
	size_t i = 0, j = 0, n = selected.size(), m = remain.size();
	double chisq = 0.0, B2 = 0.0;
	eigenVector b(n), se(n);

	if (cdat->B_N.cols() < 1) {
		if (!init_b(selected, cdat, ref)) {
			spdlog::critical("There is a collinearity problem with the SNPs given.");
			return;
		}
	}

	if (cdat->Z_N.cols() < 1) {
		init_z(selected, cdat, ref);
	}

	for (i = 0; i < n; i++) {
		b[i] = ja_beta[selected[i]];
		se[i] = ja_beta_se[selected[i]];
	}
	eigenVector bJ1 = cdat->B_N_i * cdat->D_N.asDiagonal() * b;

	eigenVector Z_Bi(n), Z_Bi_temp(n);
	bC = eigenVector::Zero(m);
	bC_se = eigenVector::Zero(m);
	pC = eigenVector::Constant(m, 2);
	for (i = 0; i < m; i++) {
		j = remain[i];
		B2 = msx[j] * nD[j];
		if (!isFloatEqual(B2, 0.0)) {
			Z_Bi = cdat->Z_N.col(j).transpose() * cdat->B_N_i;
			Z_Bi_temp = cdat->Z.col(j).transpose() * cdat->B_i;
			if (cdat->Z.col(j).dot(Z_Bi_temp) / msx_b[j] < a_collinear) {
				bC[i] = ja_beta[j] - Z_Bi.cwiseProduct(cdat->D_N).dot(b) / B2;
				bC_se[i] = 1.0 / B2; //(B2 - Z_N.col(j).dot(Z_Bi)) / (B2 * B2);
			}
		}
		bC_se[i] *= jma_Ve;
		if (bC_se[i] > 1e-10 * jma_Vp) {
			bC_se[i] = sqrt(bC_se[i]);
			chisq = bC[i] / bC_se[i];
			pC[i] = pchisq(chisq * chisq, 1);
		}
	}
}

void cond_analysis::makex_eigenVector(size_t j, eigenVector &x, bool resize, reference *ref)
{
	size_t n = fam_ids_inc.size(),
		m = to_include.size();

	if (resize)
		x.resize(n);

#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		if (!ref->bed_snp_1[to_include[j]][fam_ids_inc[i]] || ref->bed_snp_2[to_include[j]][fam_ids_inc[i]])
		{
			double snp1 = ref->bed_snp_1[to_include[j]][fam_ids_inc[i]] ? 1.0 : 0.0,
				snp2 = ref->bed_snp_2[to_include[j]][fam_ids_inc[i]] ? 1.0 : 0.0;
			if (ref->bim_allele1[to_include[j]] == ref->ref_A[to_include[j]])
				x[i] = snp1 + snp2;
			else
				x[i] = 2.0 - (snp1 + snp2);
		}
		else {
			x[i] = mu[to_include[j]];
		}
		x[i] -= mu[to_include[j]];
	}
}

bool cond_analysis::init_b(const vector<size_t> &idx, conditional_dat *cdat, reference *ref)
{
	size_t i = 0, j = 0, k = 0,
		n = fam_ids_inc.size(),
		i_size = idx.size();
	double d_temp = 0.0;
	eigenVector diagB(i_size),
		x_i(n),
		x_j(n);

	cdat->B.resize(i_size, i_size);
	cdat->B_N.resize(i_size, i_size);
	cdat->D_N.resize(i_size);

	for (i = 0; i < i_size; i++) {
		cdat->D_N[i] = msx[idx[i]] * nD[idx[i]];
		cdat->B.startVec(i);
		cdat->B.insertBack(i, i) = msx_b[idx[i]];
		
		cdat->B_N.startVec(i);
		cdat->B_N.insertBack(i, i) = cdat->D_N[i];

		diagB[i] = msx_b[idx[i]];
		makex_eigenVector(idx[i], x_i, false, ref);

		for (j = i + 1; j < i_size; j++) {
			if ((ref->bim_chr[to_include[idx[i]]] == ref->bim_chr[to_include[idx[j]]]
					&& abs(ref->bim_bp[to_include[idx[i]]] - ref->bim_bp[to_include[idx[j]]]) < a_ld_window)
				)
			{
				makex_eigenVector(idx[j], x_j, false, ref);

				d_temp = x_i.dot(x_j) / (double)n;
				cdat->B.insertBack(j, i) = d_temp;
				cdat->B_N.insertBack(j, i) = d_temp 
									* min(nD[idx[i]], nD[idx[j]]) 
									* sqrt(msx[idx[i]] * msx[idx[j]] / (msx_b[idx[i]] * msx_b[idx[j]]));
			}
		}
	}

	cdat->B.finalize();
	cdat->B_N.finalize();

	SimplicialLDLT<eigenSparseMat> ldlt_B(cdat->B);
	if (ldlt_B.vectorD().minCoeff() < 0 || sqrt(ldlt_B.vectorD().maxCoeff() / ldlt_B.vectorD().minCoeff()) > 30)
		return false;
	cdat->B_i.resize(i_size, i_size);
	cdat->B_i.setIdentity();
	cdat->B_i = ldlt_B.solve(cdat->B_i).eval();
	if ((1 - eigenVector::Constant(i_size, 1).array() / (diagB.array() * cdat->B_i.diagonal().array())).maxCoeff() > a_collinear)
		return false;

	SimplicialLDLT<eigenSparseMat> ldlt_B_N(cdat->B_N);
	cdat->B_N_i.resize(i_size, i_size);
	cdat->B_N_i.setIdentity();
	cdat->B_N_i = ldlt_B_N.solve(cdat->B_N_i).eval();
	return true;
}

void cond_analysis::init_z(const vector<size_t> &idx, conditional_dat *cdat, reference *ref)
{
	size_t i = 0, j = 0,
		n = fam_ids_inc.size(),
		m = to_include.size(),
		i_size = idx.size();
	double d_temp = 0.0;
	eigenVector x_i(n), x_j(n);

	cdat->Z.resize(i_size, m);
	cdat->Z_N.resize(i_size, m);

	for (j = 0; j < m; j++) {
		cdat->Z.startVec(j);
		cdat->Z_N.startVec(j);

		makex_eigenVector(j, x_j, false, ref);
		for (i = 0; i < i_size; i++) {
			if ((idx[i] != j 
					&& ref->bim_chr[to_include[idx[i]]] == ref->bim_chr[to_include[j]]
					&& abs(ref->bim_bp[to_include[idx[i]]] - ref->bim_bp[to_include[j]]) < a_ld_window)
				)
			{
				makex_eigenVector(idx[i], x_i, false, ref);

				d_temp = x_j.dot(x_i) / (double)n;
				cdat->Z.insertBack(i, j) = d_temp;
				cdat->Z_N.insertBack(i, j) = d_temp
										* min(nD[idx[i]], nD[j])
										* sqrt(msx[idx[i]] * msx[j] / (msx_b[idx[i]] * msx_b[j]));
			}
		}
	}

	cdat->Z.finalize();
	cdat->Z_N.finalize();
}

void cond_analysis::massoc_joint(const vector<size_t> &idx, conditional_dat *cdat, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, reference *ref)
{
	size_t i = 0, n = idx.size();
	double chisq = 0.0;
	eigenVector b(n);
	for (i = 0; i < n; i++)
		b[i] = ja_beta[idx[i]];

	if (cdat->B_N.cols() < 1) {
		if (!init_b(idx, cdat, ref)) {
			spdlog::critical("There is a collinearity problem with the given list of SNPs.");
			return;
		}
	}
	
	bJ.resize(n);
	bJ_se.resize(n);
	pJ.resize(n);
	bJ = cdat->B_N_i * cdat->D_N.asDiagonal() * b;
	bJ_se = cdat->B_N_i.diagonal();
	pJ = eigenVector::Ones(n);
	bJ_se *= jma_Ve;
	for (i = 0; i < n; i++) {
		if (bJ_se[i] > 1.0e-30) {
			bJ_se[i] = sqrt(bJ_se[i]);
			chisq = bJ[i] / bJ_se[i];
			pJ[i] = pchisq(chisq * chisq, 1.0);
		}
		else {
			bJ[i] = 0.0;
			bJ_se[i] = 0.0;
		}
	}
}

/*
 * Determine number of independent association signals within the region
 * without conducting a conditional analysis.
 */
void cond_analysis::find_independent_snps(conditional_dat *cdat, reference *ref)
{
	vector<size_t> selected, remain;
	eigenVector bC, bC_se, pC;
	jma_snpnum_backward = 0;

	if (a_top_snp <= 0.0)
		a_top_snp = 1e10;

	spdlog::info("[{}] Performing stepwise model selection on {} SNPs; p cutoff = {}, collinearity = {} assuming complete LE between SNPs more than {} Mb away).", cname, to_include.size(), a_p_cutoff, a_collinear, a_ld_window / 1e6);
	stepwise_select(selected, remain, cdat, bC, bC_se, pC, ref);

	if (selected.empty()) {
		spdlog::warn("[{}] No SNPs have been selected by the step-wise selection algorithm. Using the unconditioned dataset.", cname);

		num_ind_snps = 0;
		return;
	}
	else if (selected.size() >= fam_ids_inc.size()) {
		spdlog::warn("[{}] Too many SNPs. The number of SNPs should not be larger than the sample size.", cname);
	}

	spdlog::info("[{}] ({} SNPs eliminated by backward selection.)", cname, jma_snpnum_backward);

	num_ind_snps = selected.size();
	ind_snps = selected;
	remain_snps = remain;
}

/*
 * Run conditional analysis on marginal data for single association.
 */
void cond_analysis::pw_conditional(int pos, bool out_cond, conditional_dat *cdat, reference *ref)
{
	vector<size_t> selected(ind_snps), remain(remain_snps);
	eigenVector bC, bC_se, pC;
	size_t compare = selected[pos];

	// Move SNPs into the remain category
	if (pos >= 0) {
		size_t i = 0;
		for (i = 0; i < selected.size(); i++) {
			if (selected[i] == compare) {
				remain.push_back(selected[i]);
				erase_B_and_Z(selected, selected[i], cdat);
				selected.erase(selected.begin() + i);
				break;
			}
		}
	}

	massoc_conditional(selected, remain, cdat, bC, bC_se, pC, ref);
	if (out_cond) {
		sanitise_output(selected, remain, pos, cdat, bC, bC_se, pC, ref);
	}

	// Save in friendly format for mdata class
	snps_cond.clear();
	b_cond.clear();
	se_cond.clear();
	maf_cond.clear();
	p_cond.clear();
	n_cond.clear();
	if (ctype == coloc_type::COLOC_CC)
		s_cond.clear();

	for (size_t i = 0; i < selected.size(); i++) {
		size_t j = selected[i];
		snps_cond.push_back(ref->bim_snp_name[to_include[j]]);
		b_cond.push_back(ja_beta[j]);
		se_cond.push_back(ja_beta_se[j]);
		maf_cond.push_back(0.5 * mu[to_include[j]]);
		p_cond.push_back(ja_pval[j]);
		if (cond_ssize) {
			n_cond.push_back(nD[j]);
			if (ctype == coloc_type::COLOC_CC)
				s_cond.push_back(ja_n_cases[j]);
		}
		else {
			n_cond.push_back(nsample[j]);
			if (ctype == coloc_type::COLOC_CC)
				s_cond.push_back(ncases[j]);
		}
	}

	for (size_t i = 0; i < remain.size(); i++) {
		size_t j = remain[i];
		snps_cond.push_back(ref->bim_snp_name[to_include[j]]);
		b_cond.push_back(bC[i]);
		se_cond.push_back(bC_se[i]);
		maf_cond.push_back(0.5 * mu[to_include[j]]);
		p_cond.push_back(pC[i]);
		if (cond_ssize) {
			n_cond.push_back(nD[j]);
			if (ctype == coloc_type::COLOC_CC)
				s_cond.push_back(ja_n_cases[j]);
		}
		else {
			n_cond.push_back(nsample[j]);
			if (ctype == coloc_type::COLOC_CC)
				s_cond.push_back(ncases[j]);
		}
	}
	cond_passed = bC.size() > 0;
}

/*
 * Constructs LD matrix for a single vector of SNPs.
 * Will be NxN size, where N is number of SNPs in the vector.
 */
void cond_analysis::LD_rval(const vector<size_t> &idx, eigenMatrix &rval, conditional_dat *cdat)
{
	size_t i = 0, j = 0,
		i_size = idx.size();
	eigenVector sd(i_size);

	for (i = 0; i < i_size; i++)
		sd[i] = sqrt(msx_b[idx[i]]);

	for (j = 0; j < i_size; j++) {
		rval(j, j) = 1.0;
		for (i = j + 1; i < i_size; i++) {
			rval(i, j) = rval(j, i) = cdat->B.coeff(i, j) / sd[i] / sd[j];
		}
	}
}

/*
 * Constructs LD matrix for two vectors of SNPs.
 * Will be NxM size, where N and M are the number of SNPs in either vector.
 */
void cond_analysis::LD_rval(const vector<size_t> &v1, const vector<size_t> &v2, eigenMatrix &rval, reference *ref)
{
	size_t i = 0, j = 0, k = 0,
		n = fam_ids_inc.size(),
		v1_size = v1.size(),
		v2_size = v2.size();
	double d_temp = 0.0;
	eigenVector x_i(n),
		x_j(n);
	eigenMatrix B_ld(v1_size, v2_size);

	for (i = 0; i < v1_size; i++) {
		makex_eigenVector(v1[i], x_i, false, ref);

		for (j = 0; j < v2_size; j++) {
			if (v1[i] == v2[j]) {
				B_ld(i, j) = msx_b[v1[i]];
				continue;
			}

			if ((ref->bim_chr[to_include[v1[i]]] == ref->bim_chr[to_include[v2[j]]]
				&& abs(ref->bim_bp[to_include[v1[i]]] - ref->bim_bp[to_include[v2[j]]]) < a_ld_window)
				)
			{
				makex_eigenVector(v2[j], x_j, false, ref);
				B_ld(i, j) = x_i.dot(x_j) / (double)n;
			}
			else {
				B_ld(i, j) = 0;
			}
		}
	}

	eigenVector sd_v1(v1_size), sd_v2(v2_size);

	for (i = 0; i < v1_size; i++)
		sd_v1[i] = sqrt(msx_b[v1[i]]);

	for (i = 0; i < v2_size; i++)
		sd_v2[i] = sqrt(msx_b[v2[i]]);

	for (i = 0; i < v1_size; i++) {
		for (j = 0; j < v2_size; j++) {
			rval(i, j) = B_ld.coeff(i, j) / sd_v1[i] / sd_v2[j];
		}
	}
}

void cond_analysis::sanitise_output(vector<size_t> &selected, vector<size_t> &remain, int pos, conditional_dat *cdat, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, reference *ref)
{
	string filename;
	size_t i = 0, j = 0, k;

	filename = a_out + "." + get_cond_name();
	filename = filename + "." + ref->bim_snp_name[to_include[remain.back()]] + ".cojo";
	ofstream ofile(filename.c_str());

	if (!ofile) {
		spdlog::warn("Cannot open file {} for writing.", filename);
		return;
	}

	// LD matrix
	eigenMatrix ld(remain.size(), selected.size());
	LD_rval(remain, selected, ld, ref);

	// Header
	ofile << "Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp\tn\tfreq_geno\tbC\tbC_se\tpC";
	for (i = 0; i < selected.size(); i++) {
		ofile << "\t" << ref->bim_snp_name[to_include[selected[i]]];
	}
	ofile << endl;

	for (i = 0; i < remain.size(); i++) {
		j = remain[i];
		ofile << ref->bim_chr[to_include[j]] << "\t" << ref->bim_snp_name[to_include[j]] << "\t" << ref->bim_bp[to_include[j]] << "\t";
		ofile << ref->ref_A[to_include[j]] << "\t" << ja_freq[j] << "\t" << ja_beta[j] << "\t" << ja_beta_se[j] << "\t";
		ofile << ja_pval[j] << "\t" << nD[j] << "\t" << 0.5 * mu[to_include[j]] << "\t";
		ofile << bJ[i] << "\t" << bJ_se[i] << "\t" << pJ[i];

		// LD structure
		for (k = 0; k < selected.size(); k++) {
			ofile << "\t" << ld(i, k) * ld(i, k);
		}
		ofile << endl;
	}
	ofile.close();

#ifdef PYTHON_INC
	string plotname = a_out + "." + get_cond_name() + "." + ref->bim_snp_name[to_include[selected[0]]] + ".png";
	locus_plot(_strdup("../../python/locusplotter.py"), (char *)filename.c_str(), (char *)plotname.c_str(), (char *)(ref->bim_snp_name[to_include[selected[0]]].c_str()), ref->bim_bp[to_include[selected[0]]], ja_pval[selected[0]], 1e-25);
#endif
}

void cond_analysis::locus_plot(char *filename, char *datafile, char *to_save, char *snpname, double bp, double p, double pC)
{
#ifdef PYTHON_INC
	int argc = 8;
	char *argv[8];

	argv[0] = _strdup("../../python/locusplotter.py");
	argv[1] = filename;
	argv[2] = datafile;
	argv[3] = to_save;
	argv[4] = snpname;
	argv[5] = _strdup(to_string(bp).c_str());
	argv[6] = _strdup(to_string(p).c_str());
	argv[7] = _strdup(to_string(pC).c_str());

	wchar_t **_argv = (wchar_t **)PyMem_Malloc(sizeof(wchar_t *) * argc);
	for (int i = 0; i < argc; i++) {
		wchar_t *arg = Py_DecodeLocale(argv[i], NULL);
		_argv[i] = arg;
	}
	Py_SetProgramName(_argv[0]);

	PySys_SetArgv(argc, _argv);

	Py_Main(argc, _argv);

	for (int i = 0; i < argc; i++) {
		PyMem_RawFree(_argv[i]);
		_argv[i] = nullptr;
	}
	_argv = nullptr;
#endif

	return;
}

/*
 * Initialise matched data class from two conditional analyses
 */
mdata::mdata(cond_analysis *ca1, cond_analysis *ca2)
{
	if (!ca1->coloc_ready() || !ca2->coloc_ready()) {
		return; // TODO Handle me
	}

	size_t n = ca1->snps_cond.size();
	vector<string>::iterator it;

	// Match SNPs
	for (it = ca1->snps_cond.begin(); it != ca1->snps_cond.end(); it++) {
		vector<string>::iterator it2;
		if ((it2 = find(ca2->snps_cond.begin(), ca2->snps_cond.end(), *it)) != ca2->snps_cond.end()) 
		{
			size_t dist1 = distance(ca1->snps_cond.begin(), it),
				dist2 = distance(ca2->snps_cond.begin(), it2);
			// Some SNPs after the conditional analysis have 0 beta - need to figure out why that is
			if (ca1->b_cond[dist1] == 0 || ca2->b_cond[dist2] == 0)
				continue;
			snp_map.insert(pair<size_t, size_t>(dist1, dist2));
		}
	}

	// Extract data we want
	map<size_t, size_t>::iterator itmap = snp_map.begin();
	size_t m = snp_map.size();
	while (itmap != snp_map.end()) {
		snps1.push_back(ca1->snps_cond[itmap->first]);
		betas1.push_back(ca1->b_cond[itmap->first]);
		ses1.push_back(ca1->se_cond[itmap->first]);
		pvals1.push_back(ca1->p_cond[itmap->first]);
		mafs1.push_back(ca1->maf_cond[itmap->first]);
		ns1.push_back(ca1->n_cond[itmap->first]);
		if (ca1->get_coloc_type() == coloc_type::COLOC_CC) {
			s1.push_back(ca1->s_cond[itmap->first]);
		}

		snps2.push_back(ca2->snps_cond[itmap->second]);
		betas2.push_back(ca2->b_cond[itmap->second]);
		ses2.push_back(ca2->se_cond[itmap->second]);
		pvals2.push_back(ca2->p_cond[itmap->second]);
		mafs2.push_back(ca2->maf_cond[itmap->second]);
		ns2.push_back(ca2->n_cond[itmap->second]);
		if (ca2->get_coloc_type() == coloc_type::COLOC_CC) {
			s2.push_back(ca2->s_cond[itmap->second]);
		}

		itmap++;
	}

	type1 = ca1->get_coloc_type();
	type2 = ca2->get_coloc_type();
}

/*
 * Initialise matched data class from one conditional analysis and one phenotype
 */
mdata::mdata(cond_analysis *ca, phenotype *ph)
{
	vector<string> pheno_snps = ph->get_snp_names();
	if (!ca->coloc_ready()) {
		return; // TODO Handle me
	}

	// Match SNPs
	for (auto i : ph->matched_idx) {
		vector<string>::iterator it;
		string snp_to_find = pheno_snps[i];

		if ((it = find(ca->snps_cond.begin(), ca->snps_cond.end(), snp_to_find)) != ca->snps_cond.end()) 
		{
			size_t dist = distance(ca->snps_cond.begin(), it);
			// Some SNPs after the conditional analysis have 0 beta - need to figure out why that is
			if (ca->b_cond[dist] == 0)
				continue;
			snp_map.insert(pair<size_t, size_t>(dist, i));
		}
	}

	// Extract data we want
	map<size_t, size_t>::iterator itmap = snp_map.begin();
	size_t m = snp_map.size();
	while (itmap != snp_map.end()) {
		snps1.push_back(ca->snps_cond[itmap->first]);
		betas1.push_back(ca->b_cond[itmap->first]);
		ses1.push_back(ca->se_cond[itmap->first]);
		pvals1.push_back(ca->p_cond[itmap->first]);
		mafs1.push_back(ca->maf_cond[itmap->first]);
		ns1.push_back(ca->n_cond[itmap->first]);
		if (ca->get_coloc_type() == coloc_type::COLOC_CC) {
			s1.push_back(ca->s_cond[itmap->first]);
		}

		snps2.push_back(ph->snp_name[itmap->second]);
		betas2.push_back(ph->beta[itmap->second]);
		ses2.push_back(ph->se[itmap->second]);
		pvals2.push_back(ph->pval[itmap->second]);
		mafs2.push_back(ph->freq[itmap->second]);
		ns2.push_back(ph->n[itmap->second]);
		if (ph->get_coloc_type() == coloc_type::COLOC_CC) {
			s2.push_back(ph->n_case[itmap->second]);
		}

		itmap++;
	}

	type1 = ca->get_coloc_type();
	type2 = ph->get_coloc_type();
}
