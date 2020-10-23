#include "coloc.h"

/*
 * cond_analysis constructor
 */
coloc_analysis::coloc_analysis(mdata *mdat, string out, double pval1, double pval2, double pval3)
{
	p1 = pval1;
	p2 = pval2;
	p3 = pval3;

	ABF_1 = NULL;
	ABF_2 = NULL;
	h0 = h1 = h2 = h3 = h4 = 0.0;
	log_abf_all = log_ABF_sum = 0.0;

	matched = mdat;
	outfile = out;
	spdlog::info("Colocalisation analysis initialised with {} SNPs.", matched->snps1.size());
}

/*
 * cond_analysis default constructor
 */
coloc_analysis::coloc_analysis()
{
	matched = NULL;
	outfile = "pwcoco_coloc.out";
	
	p1 = 1e-4;
	p2 = 1e-4;
	p3 = 1e-5;

	ABF_1 = NULL;
	ABF_2 = NULL;
	h0 = h1 = h2 = h3 = h4 = 0.0;
	log_abf_all = log_ABF_sum = 0.0;
}

void coloc_analysis::estimate_bf(const vector<double> beta, const vector<double> se, const vector<double> freq, 
	const vector<double> n, vector<double> *ABF)
{
	vector<double> varbeta, invbeta,
		nvx = freq, 
		z = beta, 
		sesq = se, 
		ssize = n,
		r, log_temp;
	double sdY, sd_prior;

	// Square standard errors and ensure frequencies are MINOR allele frequencies
	transform(sesq.begin(), sesq.end(), sesq.begin(), [](double x) { return x * x; });
	transform(nvx.begin(), nvx.end(), nvx.begin(), [](double x) { return x > 0.5 ? 1.0 - x : x; });

	// Inverse beta calculation
	varbeta = sesq;
	invbeta = varbeta;
	transform(invbeta.begin(), invbeta.end(), invbeta.begin(), [](double x) { return 1 / x; });
	// Estimate sdY from 
	transform(nvx.begin(), nvx.end(), ssize.begin(), nvx.begin(), [](double x, double n_) { return 2 * n_ * x * (1 - x); });

	// Regress b*var(x) against 1/var(beta)
	sdY = lm_fixed(invbeta, nvx); // same as: lm(nvx ~ invbeta - 1)
	sdY = sqrt(sdY);

	// Calculate z and estimate Bayes factors
	transform(z.begin(), z.end(), se.begin(), z.begin(), [](double b, double se_) { return b / se_; });

	// if type == COLOC_QUANT
	sd_prior = 0.15 * sdY;
	// else sd_prior = 0.2
	
	r = varbeta;
	transform(r.begin(), r.end(), r.begin(), [sd_prior](double v) { return (sd_prior * sd_prior) / (sd_prior * sd_prior + v); });

	log_temp = varbeta;
	transform(log_temp.begin(), log_temp.end(), log_temp.begin(), [sd_prior](double x) { return log(x) - log(sd_prior * sd_prior + x); });

	// Caluclate approximate Bayes factor
	size_t i, size = varbeta.size();
	//ABF->resize(size);
	for (i = 0; i < size; i++) {
		ABF->push_back(0.5 * (log_temp[i] + (r[i] * z[i] * z[i])));
	}
}

/*
 * `combine.abf`
 * l1 = lABF1, l2 = lABF2
 * p1 = 1e-4, p2 = 1e-4, p12 = 1e-5
 */
void coloc_analysis::combine_abf(size_t abf_size)
{
	double log_ABF1 = logsum(*ABF_1),
		log_ABF2 = logsum(*ABF_2);

	h0 = 0.0;
	h1 = log(p1) + log_ABF1;
	h2 = log(p2) + log_ABF2;
	h3 = log(p1) + log(p2) + logdiff(log_ABF1 + log_ABF2, log_ABF_sum);
	h4 = log(p3) + log_ABF_sum;
	
	vector<double> h_temp = { h0, h1, h2, h3, h4 };
	log_abf_all = logsum(h_temp);
	pp_abf = h_temp;
	transform(pp_abf.begin(), pp_abf.end(), pp_abf.begin(), [=](double x) { return exp(x - log_abf_all); });
}

void coloc_analysis::init_coloc()
{
	perform_coloc();
	spdlog::info("Unconditioned colocalisation results.");
	spdlog::info("H0: {:.2f}; H1: {:.2f}; H2: {:.2f}; H3: {:.2f}; H4: {:.2f}; abf_all: {:.2f}.", pp_abf[H0], pp_abf[H1], pp_abf[H2], pp_abf[H3], pp_abf[H4], log_abf_all);
	results_to_file("unconditioned", "unconditioned");
}

void coloc_analysis::init_coloc(string snp1, string snp2)
{
	perform_coloc();
	spdlog::info("Conditioned results for SNP1: {}, SNP2: {}", snp1, snp2);
	spdlog::info("H0: {:.2f}; H1: {:.2f}; H2: {:.2f}; H3: {:.2f}; H4: {:.2f}; abf_all: {:.2f}.", pp_abf[H0], pp_abf[H1], pp_abf[H2], pp_abf[H3], pp_abf[H4], log_abf_all);
	results_to_file(snp1, snp2);
}

void coloc_analysis::perform_coloc()
{
	vector<double> temp, temp2;
	// Estimate sdY and then Bayes factor for the two datasets
	// First the conditional analysis dataset
	estimate_bf(matched->betas1, matched->ses1, matched->mafs1, matched->ns1, &temp);
	ABF_1 = &temp;

	// Now the outcome phenotype dataset
	estimate_bf(matched->betas2, matched->ses2, matched->mafs2, matched->ns2, &temp2);
	ABF_2 = &temp2;

	// "Merge" ABFs for SNPs in both datasets and sum
	size_t i, n = ABF_1->size();
	for (i = 0; i < n; i++) {
		ABF_sum.push_back(ABF_1->at(i) + ABF_2->at(i));
	}
	log_ABF_sum = logsum(ABF_sum);

	// Likely unnecessary - but this is posterior probability of each SNP's H4
	vector<double> pp_h4(n);
	for (i = 0; i < n; i++) {
		pp_h4[i] = exp(ABF_sum[i] - log_ABF_sum);
	}

	// Combine the PPs to find each H
	combine_abf(ABF_sum.size());
}

void coloc_analysis::results_to_file(string s1, string s2)
{
	ifstream ifile(outfile + ".coloc");
	bool write_header = file_is_empty(ifile);
	ofstream file;

	ifile.close();
	file.open(outfile + ".coloc", std::ios::out | std::ios::app);
	if (file.fail()) {
		spdlog::warn("Could not write colocalisation results to file {}. Please check permissions for this folder.", outfile + ".coloc");
		return;
	}
	
	if (write_header) {
		file << "SNP1\tSNP2\tH0\tH1\tH2\tH3\tH4\tlog_abf_all" << endl;
	}
	file << s1 << "\t" << s2 << "\t" << pp_abf[H0] << "\t" << pp_abf[H1] << "\t" << pp_abf[H2] << "\t" << pp_abf[H3] << "\t" << pp_abf[H4] << "\t" << log_abf_all << endl;
	file.close();
}