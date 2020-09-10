#include "data.h"

/*
 * Sub-class phenotype constructor
 */
phenotype::phenotype(string name)
{
	pheno_name = name;
	pheno_variance = 0;
}

/*
 * Sub-class phenotype default constructor
 */
phenotype::phenotype()
{
	pheno_name = "";
	pheno_variance = 0;
}

void phenotype::phenotype_clear()
{
	pheno_name = "";
	pheno_variance = 0;

	vector<string>().swap(snp_name);
	vector<string>().swap(allele1);
	vector<string>().swap(allele2);
	vector<double>().swap(freq);
	vector<double>().swap(beta);
	vector<double>().swap(se);
	vector<double>().swap(pval);
	vector<double>().swap(n);
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
phenotype *init_pheno(string filename, string pheno_name)
{
	phenotype *pheno = new phenotype(pheno_name);
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

	// Buffers
	string snp_name_buf, allele1_buf, allele2_buf;
	double freq_buf = -1.0, beta_buf, se_buf = 0.0, pval_buf, n_buf, nc_buf;

	ifstream pfile(filename.c_str());
	if (!pfile) {
		throw ("IO Error: phenotype file cannot be opened for reading: " + filename + ".");
	}
	cout << "Reading data from phenotype file: " + filename + "." << endl;
	//phenotype_clear(pheno);

	// Iterate through bim file and save data
	while (getline(pfile, line)) {
		istringstream ss(line);
		string substr;
		int tab = 0;

		// Replace buffers
		snp_name_buf = allele1_buf = allele2_buf = "";
		se_buf = 0.0; 
		freq_buf = -1;

		while (getline(ss, substr, '\t')) {
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
			case 7: // Eighth column contains N
				checkEntry(substr, &n_buf);
				break;
			}
			tab++;
		}
		if (se_buf == 0.0 || se_buf == 1.0
			|| freq_buf == -1.0 || beta_buf == 1.0)
			continue;

		// Add SNPs
		snp_name.push_back(snp_name_buf);
		allele1.push_back(allele1_buf);
		allele2.push_back(allele2_buf);
		freq.push_back(freq_buf);
		beta.push_back(beta_buf);
		se.push_back(se_buf);
		pval.push_back(pval_buf);
		n.push_back(n_buf);
		n_case.push_back(nc_buf);

		// Calculate variance of phenotype
		h = 2.0 * freq[count] * (1.0 - freq[count]);
		Vp = h * n[count] * se[count] * se[count] + h * beta[count] * beta[count] * n[count] / (n[count] - 1.0);
		if (Vp < 0.0) {
			ShowError("Error in reading phenotype file \"" + filename + "\": variance is less than zero.");
		}
		Vp_v.push_back(Vp);

		count++;
	}
	pfile.close();
	pheno_variance = v_calc_median(Vp_v);

	cout << "Read a total of: " << snp_name.size() << " lines in phenotype file \"" << filename << "\"." << endl;
	cout << "Phenotypic variance estimated from summary statistcs of all SNPs: " << pheno_variance << "." << endl;
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
	cout << "New phenotypic variances estimated from SNPs included in analysis is: " << pheno_variance << "." << endl;
	return pheno_variance;
}

/*
 * Sub-class reference data constructor
 */
reference::reference(string out, unsigned short chr, bool verbose)
{
	a_out = out;
	a_verbose = verbose;
	a_chr = chr;
}

/*
 * Sub-class phenotype default constructor
 */
reference::reference()
{
	a_out = "";
	a_verbose = false;
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
 * @ret void
 */
void reference::read_bimfile(string bimfile)
{
	string line;
	int count = 0;
	//size_t i;
	vector<string>::iterator iter;

	// Temporary containers
	int bim_chr_buf, bim_bp_buf;
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
		throw ("IO Error: bim file cannot be opened for reading: " + bimfile + ".");
	}
	cout << "Reading data from bim file: " + bimfile + "." << endl;
	bim_clear();

	
	while (bim.good()) {
		bim >> bim_chr_buf;
		if (bim.eof())
			break;
		bim_chr.push_back(bim_chr_buf);

		bim >> bim_snp_name_buf;
		bim_snp_name.push_back(bim_snp_name_buf);

		bim >> bim_genet_dst_buf;
		bim_genet_dst.push_back(bim_genet_dst_buf);

		bim >> bim_bp_buf;
		bim_bp.push_back(bim_bp_buf);
		
		bim >> bim_allele1_buf;
		transform(bim_allele1_buf.begin(), bim_allele1_buf.end(), bim_allele1_buf.begin(), ::toupper);
		bim_allele1.push_back(bim_allele1_buf);

		bim >> bim_allele2_buf;
		transform(bim_allele2_buf.begin(), bim_allele2_buf.end(), bim_allele2_buf.begin(), ::toupper);
		bim_allele2.push_back(bim_allele2_buf);

	}
	bim.close();
	
	/*
	// Iterate through bim file and save data
	while (getline(bim, line)) {
		istringstream ss(line);
		string substr;
		bool control_break = false, control_continue = false;;
		while (getline(ss, substr, '\t')) {
			switch (count) {
			case 0: // First column contains chromosome
				try {
					bim_chr_buf = stoi(substr);

					if (bim_chr_buf != a_chr && a_chr != -1) {
						control_continue = true;
						break;
					}
				}
				catch (...) {
					ShowWarning("Chromosome is not integer - skipping rest of file.", a_verbose);
					control_break = true;
				}
				break;
			case 1: // Second column contains SNP name
				bim_snp_name_buf = substr;
				break;
			case 2: // Third column contains genetic position in morgans
				bim_genet_dst_buf = stod(substr);
				break;
			case 3: // Fourth column contains genetic position on base-pairs
				bim_bp_buf = stoi(substr);
				break;
			case 4: // Fifth column contains allele 1
				bim_allele1_buf = string2upper(substr);
				break;
			case 5: // Sixth column contains allele 2
				bim_allele2_buf = string2upper(substr);
				break;
			}
			count++;

			if (control_break || control_continue)
				break;
		}

		count = 0;
		if (control_break)
			break;
		else if (control_continue)
			continue;

		// Check if this SNP matches with the phenotype SNP
		iter = find(pheno->snp_name.begin(), pheno->snp_name.end(), bim_snp_name_buf);
		if (iter == pheno->snp_name.end())
			continue;
		i = iter - pheno->snp_name.begin();
		//if (find(pheno->matched_idx.begin(), pheno->matched_idx.end(), i) == pheno->matched_idx.end()) // Matched SNPs between the two datasets ONLY!
		//	continue;

		if (bim_allele1_buf != pheno->allele1[i]
			&& bim_allele2_buf != pheno->allele1[i]
			)
		{
			bad_snp.push_back(bim_snp_name_buf);
			bad_allele1.push_back(bim_allele1_buf);
			bad_allele2.push_back(bim_allele2_buf);
			bad_refA.push_back(pheno->allele1[i]);
			continue;
		}

		// Good SNP, add it to the vectors
		bim_chr.push_back(bim_chr_buf);
		bim_snp_name.push_back(bim_snp_name_buf);
		bim_genet_dst.push_back(bim_genet_dst_buf);
		bim_bp.push_back(bim_bp_buf);
		bim_allele1.push_back(bim_allele1_buf);
		bim_allele2.push_back(bim_allele2_buf);
	}

	bim.close();

	// Report on reading
	if (!bad_snp.empty()) {
		string badname = a_out + ".badsnps";
		ofstream badfile(badname.c_str());

		badfile << "SNP\tA1\tA2\tRefA" << endl;
		for (i = 0; i < bad_snp.size(); i++) {
			badfile << bad_snp[i] << "\t" << bad_allele1[i] << "\t" << bad_allele2[i] << "\t" << bad_refA[i] << endl;
		}
		badfile.close();
		cout << "Warning: a number of SNPs from the phenotype file could not be matched to the refrence GWAS data. These SNPs are saved in \"" + badname + "\"." << endl;
	}
	*/
	num_snps = bim_chr.size();
	ref_A = bim_allele1;
	other_A = bim_allele2;
	
	cout << "Number of SNPs read from .bim file: " << num_snps << "." << endl;
}

/*
 * Matches .bim data to the matched data from the initial coloc analysis
 * Doing this early cuts down on the memory and processing footprint of the program.
 */
void reference::match_bim(vector<string> &names, vector<string> &names2)
{
	// Temporary containers
	vector<string> bim_snp_name_t,
		bim_allele1_t,
		bim_allele2_t;
	vector<int> bim_chr_t,
		bim_bp_t;
	vector<double> bim_genet_dst_t;
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
		positions[i] = (int)(it - bim_snp_name.begin());
	}
	bim_og_pos = v_remove_nans(positions);

	for (auto &p : bim_og_pos) {
		bim_snp_name_t.push_back(bim_snp_name[p]);
		bim_allele1_t.push_back(bim_allele1[p]);
		bim_allele2_t.push_back(bim_allele2[p]);
		bim_chr_t.push_back(bim_chr[p]);
		bim_bp_t.push_back(bim_bp[p]);
		bim_genet_dst_t.push_back(bim_genet_dst[p]);
	}
	bim_snp_name.swap(bim_snp_name_t);
	bim_allele1.swap(bim_allele1_t);
	bim_allele2.swap(bim_allele2_t);
	bim_chr.swap(bim_chr_t);
	bim_bp.swap(bim_bp_t);
	bim_genet_dst.swap(bim_genet_dst_t);

	num_snps_matched = bim_chr.size();
	ref_A = bim_allele1;
	other_A = bim_allele2;

	cout << "Number of SNPs matched from .bim file to the phenotype data: " << num_snps_matched << "." << endl;
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
			cout << "Duplicated SNP ID (pos = " << i << "), name: " << bim_snp_name[i] << endl;

			ss << bim_snp_name[i] << "_" << i + 1;
			bim_snp_name[i] = ss.str();

			ShowWarning("This SNP has been changed to \"" + bim_snp_name[i] + "\".", a_verbose);
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
 * @ret void
 */
void reference::read_famfile(string famfile)
{
	string line, tstr;
	int count = 0;

	// Prepare for bim reading
	ifstream fam(famfile.c_str());
	if (!fam) {
		throw ("IO Error: fam file cannot be opened for reading: " + famfile + ".");
	}
	cout << "Reading data from fam file: " + famfile + "." << endl;
	fam_clear();

	while (fam) {
		fam >> tstr;
		if (fam.eof()) break;
		fam_fid.push_back(tstr);
		fam >> tstr;
		fam_iid.push_back(tstr);
		fam >> tstr;
		fam_fa_id.push_back(tstr);
		fam >> tstr;
		fam_mo_id.push_back(tstr);
		fam >> tstr;
		fam_sex.push_back(atoi(tstr.c_str()));
		fam >> tstr;
		fam_pheno.push_back(atoi(tstr.c_str()));

	}

	/* Slow
	// Iterate through fam file and save data
	while (getline(fam, line)) {
		istringstream ss(line);
		string substr;
		while (getline(ss, substr, ' ')) {
			switch (count) {
			case 0: // First column contains family ID
				fam_fid.push_back(substr);
				break;
			case 1: // Second column contains within-family ID
				fam_iid.push_back(substr);
				break;
			case 2: // Third column contains ID of father
				fam_fa_id.push_back(substr);
				break;
			case 3: // Fourth column contains ID of mother
				fam_mo_id.push_back(substr);
				break;
			case 4: // Fifth column contains sex code
				fam_sex.push_back((unsigned short)atoi(substr.c_str()));
				break;
			case 5: // Sixth column contains phenotype value
				int temp = atoi(substr.c_str());
				if (temp != 1 && temp != 2)
					temp = 0;
				fam_pheno.push_back((unsigned short)temp);
				break;
			}
			count++;
		}
		count = 0;
	}
	*/
	fam.close();

	individuals = fam_fid.size();
	cout << "Successfully read " << individuals << " from .fam file." << endl;

	pair_fam();
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
			ShowError("Duplicate individual in .fam file found: \"" + fam_fid[i] + "\t" + fam_iid[i] + "\"."); /// Include IDs here probably
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

void reference::filter_snp_maf(double maf)
{
	map<string, size_t> id_map(snp_map);
	map<string, size_t>::iterator it, end = id_map.end();
	vector<size_t> bad_idx;
	vector<double> bad_maf;
	size_t prev_size = to_include.size();
	double f = 0.0;

	cout << "Filtering SNPs with MAF > " << maf << "." << endl;

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
	
	if (to_include.size() == 0)
		ShowError("Data Error: After MAF filtering, no SNPs are retained for the analysis.");

	//stable_sort(to_include.begin(), to_include.end());
	cout << "After filtering based on MAF, there are " << to_include.size() << " SNPs remaining in the analysis (amount removed: " << prev_size - to_include.size() << ")." << endl;
}

/*
 * Calculates allele frequencies based on .fam data
 */
void reference::calculate_allele_freq()
{
	const size_t n = to_include.size(),
		m = fam_ids_inc.size();

	cout << "Calculating allele frequencies from .fam data." << endl;
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
 * This function will read the allelic information
 * from the .bed file and use it to calculate allele
 * frequencies for the reference data.
 * @param string bedfile File name and path to .bed file
 * @ret void
 */
void reference::read_bedfile(string bedfile)
{
	size_t i, j, k,
		snp_idx, ind_idx;
	const size_t bim_size = to_include.size(),
		fam_size = fam_ids_inc.size();
	vector<int> read_individuals;
	
	char ch[1];
	bitset<8> b;
	fstream BIT(bedfile.c_str(), ios::in | ios::binary);
	cout << "Reading .bed file..." << endl;

	get_read_individuals(read_individuals);

	if (!bim_size || !fam_size) {
		ShowError("No data from either .bim or .fam files to continue with.");
		return;
	}

	if (!BIT) {
		ShowError("Cannot open " + bedfile + " to read.");
		return;
	}

	vector<vector<bool>>snp_1(bim_size, vector<bool>(fam_size, false));
	vector<vector<bool>>snp_2(bim_size, vector<bool>(fam_size, false));

	//bed_snp_1.resize(bim_size);
	//bed_snp_2.resize(bim_size);
	//for (i = 0; i < bim_size; i++) {
	//	bed_snp_1[i].resize(fam_size);
	//	bed_snp_2[i].resize(fam_size);
	//}

	for (i = 0; i < 3; i++) {
		BIT.read(ch, 1); // First three bytes are used for the file format and are not read
	}

	for (i = 0, snp_idx = 0; i < num_snps; i++) {
		// 00: homozygote AA
		// 11: homozygote BB
		// 10: missing
		// 01: heterozygous
		
		// SNP not found in the matched SNP list
		vector<size_t>::iterator it;
		if ((it = find(to_include.begin(), to_include.end(), i)) == to_include.end()) {
			for (j = 0; j < individuals; j += 4)
				BIT.read(ch, 1);
			continue;
		}

		snp_idx = it - to_include.begin();
		// SNP found - calculator allele frequency
		for (j = 0, ind_idx = 0; j < individuals;) {
			BIT.read(ch, 1);
			if (!BIT)
				ShowError("Cannot read the .bed file?");

			b = ch[0];
			k = 0;
			while (k < 7 && j < individuals) { // 11 for AA; 00 for BB
				if (read_individuals[j] == 0)
					k += 2;
				else {
					snp_2[snp_idx][ind_idx] = (!b[k++]);
					snp_1[snp_idx][ind_idx] = (!b[k++]);
					ind_idx++;
				}
				j++;
			}
		}

		//if (snp_idx == bim_size)
		//	break;
		//snp_idx++;
	}

	BIT.clear();
	BIT.close();

	bed_snp_1.swap(snp_1);
	bed_snp_2.swap(snp_2);

	cout << "Genotype data for " << fam_size << " individuals and " << bim_size << " SNPs read." << endl;
}

/*
 * Updates the inclusion list of SNPs based on an index-based vector.
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
				// Clean your own damn data!
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

		snps2.push_back(ph2->snp_name[itmap->second]);
		betas2.push_back(ph2->beta[itmap->second]);
		ses2.push_back(ph2->se[itmap->second]);
		pvals2.push_back(ph2->pval[itmap->second]);
		mafs2.push_back(ph2->freq[itmap->second]);
		ns2.push_back(ph2->n[itmap->second]);

		itmap++;
	}
}

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
}
