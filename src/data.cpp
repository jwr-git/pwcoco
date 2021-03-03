#include "data.h"

/*
 * Phenotype constructor
 */
phenotype::phenotype(string name, double n, double n_case)
{
	pheno_name = name;
	pheno_variance = 0;
	failed = false;
	ctype = coloc_type::COLOC_NONE;
	n_from_cmd = n;
	n_case_from_cmd = n_case;
}

/*
 * Phenotype default constructor
 */
phenotype::phenotype()
{
	pheno_name = "";
	pheno_variance = 0;
	failed = false;
	ctype = coloc_type::COLOC_NONE;
	n_from_cmd = 0;
	n_case_from_cmd = 0;
}

void phenotype::phenotype_clear()
{
	pheno_name = "";
	pheno_variance = 0;
	failed = false;
	ctype = coloc_type::COLOC_NONE;
	n_from_cmd = 0;
	n_case_from_cmd = 0;

	vector<string>().swap(snp_name);
	vector<string>().swap(allele1);
	vector<string>().swap(allele2);
	vector<double>().swap(freq);
	vector<double>().swap(beta);
	vector<double>().swap(se);
	vector<double>().swap(pval);
	vector<double>().swap(n);
	vector<double>().swap(n_case);
	vector<double>().swap(Vp_v);
}

/*
 * Initialises the "exposure" phenotype class.
 * This phenotype is called exposure for ease of differentiation between the
 * other phenotype ("outcome").
 * This is phenotype 1 in the code.
 * @param cond_analysis *p_analysis Pointed to conditional analysis class
 * @param string filename Name of the phenotype file which should follow
 * the following format:
 * SNP	Allele1	Allele2	Freq	Beta	SE	Pvalue	N
 * The data MUST follow this order.
 * The header ONLY has to begin with the word "SNP" - the other column
 * names are flexible.
 * @ret void
 */
phenotype *init_pheno(string filename, string pheno_name, double n, double n_case)
{
	phenotype *pheno = new phenotype(pheno_name, n, n_case);
	pheno->read_phenofile(filename);
	return pheno;
}

/*
 * Reads the phenotypical file (regardless if "exposure" or "outcome") and formats the data accordingly.
 * The data will undergo linkage with the genotypic data from the .bim file.
 * This is done by matching SNPs from the phenotype to the genotype based on alleles.
 * This means the data all MUST use the same genomic build. If the SNP in the phenotype file is not found, then it will
 * be output in a file for the user to use as they want.
 * @param cond_analysis *p_analysis Class of current conditional analysis
 * @param phenotype *pheno The phenotype that is being currently read
 * @param string filename Path to the phenotype file
 * @ret void
 */
void phenotype::read_phenofile(string filename)
{
	string line;
	map<string, int>::iterator iter;
	int count = 0;
	double h = 0.0, Vp = 0.0;
	bool delim_found = false, n_warning = false;
	char sep = '\t';

	// Buffers
	string snp_name_buf, allele1_buf, allele2_buf;
	double freq_buf = -1.0, beta_buf, se_buf = 0.0, pval_buf, n_buf, nc_buf = 0.0, s;

	ifstream pfile(filename.c_str());
	if (!pfile) {
		spdlog::critical("Phenotype file cannot be opened for reading: {}.", filename);
		failed = true;
		return;
	}
	spdlog::info("Reading data from phenotype file: {}.", filename);
	//phenotype_clear(pheno);

	// Iterate through file and save data
	while (getline(pfile, line)) {
		istringstream ss(line);
		string substr;
		int tab = 0;

		// Replace buffers
		snp_name_buf = allele1_buf = allele2_buf = "";
		se_buf = 0.0; 
		freq_buf = -1;

		// Check delimiter
		if (ss.str().find('\t', 0) != string::npos) {
			sep = '\t';
			delim_found = true;
		}
		else if (ss.str().find(',', 0) != string::npos) {
			sep = ',';
			delim_found = true;
		}
		else if (ss.str().find(' ', 0) != string::npos) {
			sep = ' ';
			delim_found = true;
		}

		if (!delim_found) {
			spdlog::critical("Cannot determine delimiter for phenotype file \"{}\". Use either tab, comma or space-delimited.", filename);
			failed = true;
			return;
		}

		while (getline(ss, substr, sep)) {
			if (string2upper(substr) == "SNP") // Skip header
				break;

			switch (tab) {
			case 0: // First column contains SNP name
				if (substr == "." || substr == "NA")
					break;
				snp_name_buf = substr;
				break;
			case 1: // Second column contains allele 1
			case 2: // Third column contains allele 2
				if (substr == "." || substr == "NA")
					break;
				if (tab == 1)
					allele1_buf = string2upper(substr);
				else
					allele2_buf = string2upper(substr);
				break;
			case 3: // Fourth column contains allele frequency of A1
				checkEntry(substr, &freq_buf);
				break;
			case 4: // Fifth column contains beta
				checkEntry(substr, &beta_buf);
				break;
			case 5: // Sixth column contains SE
				checkEntry(substr, &se_buf);
				break;
			case 6: // Seventh column contains P value
				checkEntry(substr, &pval_buf);
				break;
			case 7: // Eighth column contains N total
				checkEntry(substr, &n_buf);
				if (n_buf && n_from_cmd && !n_warning) {
					n_warning = true;
					spdlog::warn("N is given at command line (n = {}) but the file also contains these values. Using values from command line instead.", n_from_cmd);
					n_buf = n_from_cmd;
				}
				else if (n_buf && n_from_cmd) {
					n_buf = n_from_cmd;
				}
				break;
			case 8: // Nineth column contains N of cases (optional)
				checkEntry(substr, &nc_buf);
				if (nc_buf && n_case_from_cmd && !n_warning) {
					n_warning = true;
					spdlog::warn("N or n_case is given at command line (n_case = {}) but the file also contains these values. Using values from command line instead.", n_case_from_cmd);
					nc_buf = n_case_from_cmd;
				}
				else if (nc_buf && n_case_from_cmd) {
					nc_buf = n_case_from_cmd;
				}
				ctype = coloc_type::COLOC_CC;
				break;
			}
			tab++;
		}
		if (se_buf == 0.0 || se_buf == 1.0
			|| freq_buf == -1.0 || beta_buf == 1.0)
			continue;

		if (!n_buf && n_from_cmd)
			n_buf = n_from_cmd;
		if (!nc_buf && n_case_from_cmd)
			nc_buf = n_case_from_cmd;

		// Add SNPs
		snp_name.push_back(snp_name_buf);
		allele1.push_back(allele1_buf);
		allele2.push_back(allele2_buf);
		freq.push_back(freq_buf);
		beta.push_back(beta_buf);
		se.push_back(se_buf);
		pval.push_back(pval_buf);
		n.push_back(n_buf);
		s = n_buf != 0 ? nc_buf / n_buf : nc_buf;
		if (s < 0 || s >= 1) {
			spdlog::warn("SNP {} in phenotyle file {} has case proportion outside of range [0, 1) - capping.", snp_name_buf, s);
			s = s < 0 ? 0 : s >= 1 ? s = 1 - 1e-8 : s;
		}
		n_case.push_back(s); // n_case will contain proportion of cases to total sample size

		// Calculate variance of phenotype
		h = 2.0 * freq[count] * (1.0 - freq[count]);
		Vp = h * n[count] * se[count] * se[count] + h * beta[count] * beta[count] * n[count] / (n[count] - 1.0);
		if (Vp < 0.0) {
			spdlog::critical("Error in reading phenotype file {}: variance is less than zero (Vp = {:.2f}).", filename, Vp);
			failed = true;
		}
		Vp_v.push_back(Vp);

		count++;
	}

	if (ctype == coloc_type::COLOC_NONE) {
		ctype = coloc_type::COLOC_QUANT; // n contains sample size
	}

	pfile.close();
	pheno_variance = v_calc_median(Vp_v);

	spdlog::info("Read a total of: {} lines in phenotype file {}.", snp_name.size(), filename);
	spdlog::info("Phenotypic variance estimated from summary statistcs of all SNPs: {:.2f}", pheno_variance);
}

/*
 * Force calculates the Phenotypic variance
 */
double phenotype::calc_variance(vector<size_t> idx)
{
	vector<double> Vp_temp;
	for (auto &i : idx) {
		Vp_temp.push_back(Vp_v[i]);
	}
	pheno_variance = v_calc_median(Vp_temp);
	spdlog::info("New phenotypic variances estimated from SNPs included in analysis is: {:.2f}.", pheno_variance);
	return pheno_variance;
}

/*
 * Reference data constructor
 */
reference::reference(string out, unsigned short chr)
{
	a_out = out;
	a_chr = chr;
}

/*
 * Reference data default constructor
 */
reference::reference()
{
	a_out = "";
	a_chr = -1;
}

void reference::reference_clear()
{
	vector<double>().swap(mu);
	bim_clear();
	fam_clear();
}

/*
 * This function will read genotypic data from the
 * provided .bim file from Plink. This function will
 * store all data from that .bim file for use in
 * the conditional analysis.
 * @param string bimfile File name and path to .bim file
 * @ret int 0 if failed, 1 if successful
 */
int reference::read_bimfile(string bimfile)
{
	string line;
	size_t i = 0;
	vector<string>::iterator iter;

	// Temporary containers
	unsigned short bim_chr_buf;
	int bim_bp_buf;
	double bim_genet_dst_buf;
	string bim_snp_name_buf = "", bim_allele1_buf = "", bim_allele2_buf = "";

	// Rubbish SNPs
	vector<string> bad_snp;
	vector<string> bad_allele1;
	vector<string> bad_allele2;
	vector<string> bad_refA;

	// Prepare for bim reading
	ifstream bim(bimfile.c_str());
	if (bim.fail()) {
		spdlog::critical("Bim file cannot be opened for reading: {}.", bimfile);
		failed = true;
		return 0;
	}
	spdlog::info("Reading data from bim file: {}.", bimfile);
	bim_clear();
	
	while (bim.good()) {
		bim >> bim_chr_buf >> bim_snp_name_buf >> bim_genet_dst_buf >> bim_bp_buf >> bim_allele1_buf >> bim_allele2_buf;

		if (bim.eof())
			break;

		if (bim_chr_buf != a_chr) {
			i++;
			continue;
		}

		bim_chr.push_back(bim_chr_buf);
		bim_snp_name.push_back(bim_snp_name_buf);
		//bim_genet_dst.push_back(bim_genet_dst_buf);
		bim_bp.push_back(bim_bp_buf);
		transform(bim_allele1_buf.begin(), bim_allele1_buf.end(), bim_allele1_buf.begin(), ::toupper); // @TODO Is it quicker to just apply this to the entire vector after?
		bim_allele1.push_back(bim_allele1_buf);
		transform(bim_allele2_buf.begin(), bim_allele2_buf.end(), bim_allele2_buf.begin(), ::toupper);
		bim_allele2.push_back(bim_allele2_buf);
		bim_og_pos.push_back(i++);
	}
	bim.close();

	num_snps = bim_chr.size();
	ref_A = bim_allele1;
	other_A = bim_allele2;
	
	spdlog::info("Number of SNPs read from .bim file: {}.", num_snps);
	return 1;
}

/*
 * Matches .bim data to the matched data from the initial coloc analysis
 * Doing this early cuts down on the memory and processing footprint of the program.
 * @param vector<string> names SNP names from first dataset
 * @param vector<string> names2 SNP names from second dataset
 * @ret void
 */
void reference::match_bim(vector<string> &names, vector<string> &names2)
{
	// Temporary containers
	vector<string> bim_snp_name_t,
		bim_allele1_t,
		bim_allele2_t;
	vector<unsigned char> bim_chr_t;
	vector<int> bim_bp_t;
	//vector<double> bim_genet_dst_t;
	vector<int> positions;
	vector<string> snp_names = names;
	
	copy(names2.begin(), names2.end(), back_inserter(snp_names));
	v_remove_dupes(snp_names);

	positions.resize(snp_names.size());
	fill(positions.begin(), positions.end(), -1);
#pragma omp parallel for
	for (int i = 0; i < snp_names.size(); ++i) {
		vector<string>::iterator it;
		if ((it = find(bim_snp_name.begin(), bim_snp_name.end(), snp_names[i])) == bim_snp_name.end())
			continue;
		positions[i] = bim_og_pos[it - bim_snp_name.begin()];
	}
	bim_og_pos = v_remove_nans(positions);

	for (auto &p : bim_og_pos) {
		bim_snp_name_t.push_back(bim_snp_name[p]);
		bim_allele1_t.push_back(bim_allele1[p]);
		bim_allele2_t.push_back(bim_allele2[p]);
		bim_chr_t.push_back(bim_chr[p]);
		bim_bp_t.push_back(bim_bp[p]);
		//bim_genet_dst_t.push_back(bim_genet_dst[p]);
	}
	bim_snp_name.swap(bim_snp_name_t);
	bim_allele1.swap(bim_allele1_t);
	bim_allele2.swap(bim_allele2_t);
	bim_chr.swap(bim_chr_t);
	bim_bp.swap(bim_bp_t);
	//bim_genet_dst.swap(bim_genet_dst_t);

	num_snps_matched = bim_chr.size();
	ref_A = bim_allele1;
	other_A = bim_allele2;

	spdlog::info("Number of SNPs matched from .bim file to the phenotype data: {}.", num_snps_matched);
}

/*
 * Clear any .bim file-related data stored
 */
void reference::bim_clear() {
	bim_chr.clear();
	bim_snp_name.clear();
	bim_genet_dst.clear();
	bim_bp.clear();
	bim_allele1.clear();
	bim_allele2.clear();
}

/*
 * Helper function called from read_bimfile which will sanitise
 * the SNP list from the .bim file. This includes renaming any
 * duplicated SNPs.
 * @ret void
 */
void reference::sanitise_list()
{
	size_t i;
	to_include.clear();
	to_include.resize(num_snps_matched);
	snp_map.clear();

	for (i = 0; i < num_snps_matched; ++i) {
		stringstream ss;
		to_include[i] = bim_og_pos[i];

		if (snp_map.find(bim_snp_name[i]) != snp_map.end()) {
			spdlog::info("Duplicated SNP ID (pos = {}), name: {}.", i, bim_snp_name[i]);

			ss << bim_snp_name[i] << "_" << i + 1;
			bim_snp_name[i] = ss.str();

			spdlog::info("This SNP has been changed to {}.", bim_snp_name[i]);
		}
		else
			ss << bim_snp_name[i];

		snp_map.insert(pair<string, size_t>(ss.str(), i));
	}
	//stable_sort(to_include.begin(), to_include.end());
}

/*
 * This function will read individual data from the
 * provided .fam file from Plink. This function will
 * store all data from that .fam file for use in
 * the conditional analysis.
 * @param string famfile File name and path to .fam file
 * @ret int 0 if failed, 1 if successful
 */
int reference::read_famfile(string famfile)
{
	string line, fid_buf, iid_buf, fa_buf, ma_buf;
	unsigned short sex_buf, pheno_buf;
	int count = 0;

	// Prepare for bim reading
	ifstream fam(famfile.c_str());
	if (!fam) {
		spdlog::critical("Fam file cannot be opened for reading: {}.", famfile);
		failed = true;
		return 0;
	}
	spdlog::info("Reading data from fam file: {}.", famfile);
	fam_clear();

	while (fam) {
		fam >> fid_buf >> iid_buf >> fa_buf >> ma_buf >> sex_buf >> pheno_buf;
		if (fam.eof()) 
			break;
		fam_fid.push_back(fid_buf);
		fam_iid.push_back(iid_buf);
		fam_fa_id.push_back(fa_buf);
		fam_mo_id.push_back(ma_buf);
		fam_sex.push_back(sex_buf);
		fam_pheno.push_back(pheno_buf);
	}

	fam.close();

	individuals = fam_fid.size();
	spdlog::info("Successfully read {} from .fam file.", individuals);

	pair_fam();
	return 1;
}

/*
 * Pairs FIDs and IIDs from .fam file and builds unique map
 * Used for later analysis.
 * @ret void
 */
void reference::pair_fam()
{
	size_t i = 0, size = 0;

	fam_ids_inc.clear();
	fam_ids_inc.resize(individuals);
	fam_map.clear();

	for (i = 0; i < individuals; i++) {
		fam_ids_inc[i] = i;
		fam_map.insert(pair<string, size_t>(fam_fid[i] + ":" + fam_iid[i], i));
		if (size == fam_map.size())
			spdlog::warn("Duplicate individual in .fam file found: {}, {}.", fam_fid[i], fam_iid[i]); /// Include IDs here probably
		size = fam_map.size();
	}
}

/*
 * Clear any .fam file-related data stored
 */
void reference::fam_clear() {
	fam_fid.clear();
	fam_iid.clear();
	fam_fa_id.clear();
	fam_mo_id.clear();
	fam_sex.clear();
	fam_pheno.clear();
}

/*
 * Filters user datasets by MAF according to a user-defined
 * threshold. 
 * @param double maf Minor allele frequency threshold
 * @ret int 1 if successful, 0 if failed
 */
int reference::filter_snp_maf(double maf)
{
	map<string, size_t> id_map(snp_map);
	map<string, size_t>::iterator it, end = id_map.end();
	vector<size_t> bad_idx;
	vector<double> bad_maf;
	size_t prev_size = to_include.size();
	double f = 0.0;

	spdlog::info("Filtering SNPs with MAF > {}.", maf);

	to_include.clear();
	snp_map.clear();

	for (it = id_map.begin(); it != end; ++it) {
		f = mu[it->second] * 0.5;
		if (f <= maf || (1.0 - f) <= maf) {
			bad_idx.push_back(it->second);
			bad_maf.push_back(f);
			continue;
		}
		snp_map.insert(*it);
		to_include.push_back(it->second);
	}
	
	if (to_include.size() == 0) {
		spdlog::critical("After MAF filtering, no SNPs are retained for the analysis.");
		return 0;
	}

	//stable_sort(to_include.begin(), to_include.end());
	spdlog::info("After filtering based on MAF, there are {} SNPs remaining in the analysis (amount removed: {}).", to_include.size(), prev_size - to_include.size());
	return 1;
}

/*
 * Calculates allele frequencies based on .fam data
 * @ret void
 */
void reference::calculate_allele_freq()
{
	const size_t n = to_include.size(),
		m = fam_ids_inc.size();

	spdlog::info("Calculating allele frequencies from .fam data.");
	mu.clear();
	mu.resize(n);

#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		double fcount = 0.0;
		for (int j = 0; j < m; j++) {
			double f = 0.0;
			if (!bed_snp_1[i][fam_ids_inc[j]] || bed_snp_2[i][fam_ids_inc[j]]) {
				double snp1 = bed_snp_1[i][fam_ids_inc[j]] ? 1.0 : 0.0,
					snp2 = bed_snp_2[i][fam_ids_inc[j]] ? 1.0 : 0.0;
				f = snp1 + snp2;
				if (bim_allele2[i] == ref_A[i])
					f = 2.0 - f;
				mu[i] += f;
				fcount += 1.0;
			}
		}

		if (fcount > 0.0)
			mu[i] /= fcount;
	}
}

/*
 * Creates read_individuals vector from .fam data.
 * @param vector<int> &read_individuals Vector reference to input data
 * @ret void
 */
void reference::get_read_individuals(vector<int> &read_individuals)
{
	read_individuals.clear();
	read_individuals.resize(individuals);
	for (size_t i = 0; i < individuals; i++) {
		if (fam_map.find(fam_fid[i] + ":" + fam_iid[i]) != fam_map.end())
			read_individuals[i] = 1;
		else
			read_individuals[i] = 0;
	}
}

/* 
 * Sub-function that work asynchronously to read the .bed buffer for faster
 * reading.
 * Only available for non-MS compilers due to VS not having up-to-date OpenMP.
 * @param char *buf Buffer
 * @param size_t snp_idx Index of the current SNP being read
 * @param vector<int> read_individuals Which samples have been read and included.
 * @ret void
 */
void reference::parse_bed_data(char *buf, size_t snp_idx, vector<int> read_individuals)
{
	size_t j, k, ind_idx;
	int buf_count = 0, bytes = (int)ceil(individuals / 4.0);
	bitset<8> b;

	for (j = 0, ind_idx = 0; j < individuals && buf_count < bytes; ) {
		b = buf[buf_count];
		k = 0;
		while (k < 7 && j < individuals) { // 11 for AA; 00 for BB
			if (read_individuals[j] == 0)
				k += 2;
			else {
				snp_2[snp_idx][ind_idx] = !b[k++];
				snp_1[snp_idx][ind_idx] = !b[k++];
				ind_idx++;
			}
			j++;
		}
		buf_count++;
	}
}

/*
 * Function that produces the .bed file into buffers and hands these off
 * to the sub-function to process asynchronously.
 * Only available for non-MS compilers due to VS not having up-to-date OpenMP.
 * @param string bedfile Path to bedfile
 * @ret int 0 if failed, 1 if successful
 */

int reference::read_bedfile(string bedfile)
{
	size_t i;
	const size_t bim_size = to_include.size(),
		fam_size = fam_ids_inc.size();
	vector<int> read_individuals;

	snp_1.resize(bim_size);
	snp_2.resize(bim_size);
	for (i = 0; i < bim_size; i++) {
		snp_1[i].resize(fam_size);
		snp_2[i].resize(fam_size);
	}

	// First three bytes are used for the file format and are not read
	char ch[3];
	fstream BIT(bedfile.c_str(), ios::in | ios::binary);	
	if (!BIT) {
		throw("Bed file cannot be opened for reading: {}.", bedfile);
		return 0;
	}
	BIT.read(ch, 3); 

	get_read_individuals(read_individuals);

	// Producer/consumer async to quickly read and process the bed file
	char* buf;
	int bytes = (int)ceil(individuals / 4.0);
#ifndef _MSC_VER
#pragma omp parallel shared(read_individuals)
#pragma omp for ordered
	for (int ii = 0; ii < num_snps; ii++) {
#pragma omp ordered
		{
			buf = new char[bytes];
			BIT.read(buf, bytes);

			// SNP in the matched SNP list?
			vector<size_t>::iterator it;
			if ((it = find(to_include.begin(), to_include.end(), ii)) != to_include.end())
			{
				size_t snp_idx = it - to_include.begin();
#pragma omp task
				parse_bed_data(buf, snp_idx, read_individuals);
			}
		}
	}
#pragma omp taskwait
#else
	for (int ii = 0; ii < num_snps; ii++) {
		{
			buf = new char[bytes];
			BIT.read(buf, bytes);

			// SNP in the matched SNP list?
			vector<size_t>::iterator it;
			if ((it = find(to_include.begin(), to_include.end(), ii)) != to_include.end())
			{
				size_t snp_idx = it - to_include.begin();
				parse_bed_data(buf, snp_idx, read_individuals);
			}
		}
	}
#endif

	BIT.clear();
	BIT.close();

	bed_snp_1.swap(snp_1);
	bed_snp_2.swap(snp_2);

	spdlog::info("Genotype data for {} individuals and {} SNPs read.", fam_size, bim_size);
	return 1;
}

/*
 * Updates the inclusion list of SNPs based on an index-based vector.
 * @param const vector<size_t> idx Index of SNPs
 * @param const vector<string> snps SNP names
 * @ret void
 */
void reference::update_inclusion(const vector<size_t> idx, const vector<string> snps)
{
	size_t i, n = idx.size();
	map<string, size_t> snp_map_buffer(snp_map), snp_map_keep(snp_map);
	map<string, size_t>::iterator iter;

	for (i = 0; i < n; i++) {
		snp_map_buffer.erase(snps[i]);
	}

	for (iter = snp_map_buffer.begin(); iter != snp_map_buffer.end(); iter++) {
		snp_map_keep.erase(iter->first);
	}

	to_include.clear();
	for (iter = snp_map_keep.begin(); iter != snp_map_keep.end(); iter++) {
		to_include.push_back(iter->second);
	}
	//stable_sort(to_include.begin(), to_include.end());
}

/*
 * Initialise matched data class from two phenotypes
 */
mdata::mdata(phenotype *ph1, phenotype *ph2)
{
	size_t n = ph1->snp_name.size();
	vector<string>::iterator it;

	// Match SNPs
	for (it = ph1->snp_name.begin(); it != ph1->snp_name.end(); it++) {
		size_t d1 = distance(ph1->snp_name.begin(), it);
		if (ph1->beta[d1] == 0.0 || isnan(ph1->beta[d1]) || isinf(ph1->beta[d1]) || !isfinite(ph1->beta[d1])
			|| isnan(ph1->se[d1]) || isinf(ph1->se[d1]) || !isfinite(ph1->se[d1]))
		{
			continue;
		}
		vector<string>::iterator it2;

		if ((it2 = find(ph2->snp_name.begin(), ph2->snp_name.end(), *it)) != ph2->snp_name.end()) {
			size_t d2 = distance(ph2->snp_name.begin(), it2);
			if (ph2->beta[d2] == 0.0 || isnan(ph2->beta[d2]) || isinf(ph2->beta[d2]) || !isfinite(ph2->beta[d2])
				|| isnan(ph2->se[d2]) || isinf(ph2->se[d2]) || !isfinite(ph2->se[d2]))
			{
				// Clean your own data!
				continue;
			}

			snp_map.insert(pair<size_t, size_t>(d1, d2));
			ph1->matched_idx.push_back(distance(ph1->snp_name.begin(), it));
			ph2->matched_idx.push_back(distance(ph2->snp_name.begin(), it2));
		}
	}

	// Extract data we want
	map<size_t, size_t>::iterator itmap = snp_map.begin();
	size_t m = snp_map.size();
	while (itmap != snp_map.end()) {
		snps1.push_back(ph1->snp_name[itmap->first]);
		betas1.push_back(ph1->beta[itmap->first]);
		ses1.push_back(ph1->se[itmap->first]);
		pvals1.push_back(ph1->pval[itmap->first]);
		mafs1.push_back(ph1->freq[itmap->first]);
		ns1.push_back(ph1->n[itmap->first]);
		s1.push_back(ph1->n_case[itmap->first]);

		snps2.push_back(ph2->snp_name[itmap->second]);
		betas2.push_back(ph2->beta[itmap->second]);
		ses2.push_back(ph2->se[itmap->second]);
		pvals2.push_back(ph2->pval[itmap->second]);
		mafs2.push_back(ph2->freq[itmap->second]);
		ns2.push_back(ph2->n[itmap->second]);
		s2.push_back(ph2->n_case[itmap->second]);

		itmap++;
	}
	type1 = ph1->get_coloc_type();
	type2 = ph2->get_coloc_type();
}

/*
 * Matched data default constructor
 */
mdata::mdata()
{
	snps1.clear();
	snps2.clear();
	betas1.clear();
	betas2.clear();
	ses1.clear();
	ses2.clear();
	pvals1.clear();
	pvals2.clear();
	mafs1.clear();
	mafs2.clear();
	ns1.clear();
	ns2.clear();
	type1 = type2 = coloc_type::COLOC_NONE;
}
