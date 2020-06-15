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

void phenotype::phenotype_clear(phenotype *pheno)
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
phenotype *init_exposure(string filename, string pheno_name)
{
	phenotype *pheno = new phenotype("pheno_name");
	pheno->read_phenofile(pheno, filename);

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
void phenotype::read_phenofile(phenotype *pheno, string filename)
{
	string line;
	map<string, int>::iterator iter;
	int count = 0, i;
	double h = 0.0, Vp = 0.0;

	// Buffers
	string snp_name_buf, allele1_buf, allele2_buf;
	double freq_buf = -1.0, beta_buf, se_buf, pval_buf, n_buf;

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
		int tab = 0, j;

		// Replace buffers
		snp_name_buf = allele1_buf = allele2_buf = "";
		freq_buf = -1;

		while (getline(ss, substr, '\t')) {
			if (substr == "SNP") // Skip header
				break;

			switch (tab) {
			case 0: // First column contains SNP name
				if (substr == "." || substr == "NA")
					break;
				snp_name_buf = substr;
				break;
			case 1: // Second column contains allele 1
			case 2: // Third column contains allele 2
				transform(substr.begin(), substr.end(), substr.begin(), ::toupper);
				if (substr == "." || substr == "NA")
					break;
				if (tab == 1)
					allele1_buf = substr;
				else
					allele2_buf = substr;
				break;
			case 3: // Fourth column contains allele frequency of A1
				freq_buf = atof(substr.c_str());
				break;
			case 4: // Fifth column contains beta
				beta_buf = atof(substr.c_str());
				break;
			case 5: // Sixth column contains SE
				se_buf = atof(substr.c_str());
				break;
			case 6: // Seventh column contains P value
				pval_buf = atof(substr.c_str());
				break;
			case 7: // Eighth column contains N
				n_buf = atof(substr.c_str());
				break;
			}
			tab++;
		}
		if (freq_buf == -1.0)
			continue;

		// Add SNPs
		pheno->snp_name.push_back(snp_name_buf);
		pheno->allele1.push_back(allele1_buf);
		pheno->allele2.push_back(allele2_buf);
		pheno->freq.push_back(freq_buf);
		pheno->beta.push_back(beta_buf);
		pheno->se.push_back(se_buf);
		pheno->pval.push_back(pval_buf);
		pheno->n.push_back(n_buf);

		// Calculate variance of phenotype
		h = 2.0 * pheno->freq[count] * (1.0 - pheno->freq[count]);
		Vp = h * pheno->n[count] * pow(pheno->se[count], 2) + h * pow(pheno->beta[count], 2) * pheno->n[count] / (pheno->n[count] - 1.0);
		if (Vp < 0.0) {
			throw ("Error in reading phenotype file \"" + filename + "\": variance is less than zero.");
		}
		pheno->Vp_v.push_back(Vp);
		pheno->Vp_v.push_back(Vp);

		count++;
	}
	pfile.close();
	pheno_variance = v_calc_median(pheno->Vp_v);
	// If want to adjust P-values, add GC here

	cout << "Read a total of: " << pheno->snp_name.size() << " lines in phenotype file \"" << filename << "\"." << endl;
	cout << "Phenotypic variance estimated from summary statistcs of all SNPs: " << pheno_variance << "." << endl;
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
	vector<double>().swap(mu); /// What am I?
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
void reference::read_bimfile(string bimfile, phenotype *pheno)
{
	string line;
	int count = 0, i;
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
	if (!bim) {
		throw ("IO Error: bim file cannot be opened for reading: " + bimfile + ".");
	}
	cout << "Reading data from bim file: " + bimfile + "." << endl;
	bim_clear();

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
				bim_allele1_buf = substr;
				break;
			case 5: // Sixth column contains allele 2
				bim_allele2_buf = substr;
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
		if (bim_allele1_buf != pheno->allele1[i]
			&& bim_allele2_buf != pheno->allele2[i])
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
	// First bad SNPs so we can free that memory
	if (!bad_snp.empty()) {
		string badname = a_out + ".badsnps";
		ofstream badfile(badname.c_str());

		badfile << "SNP\tA1\tA2\tRefA" << endl;
		for (i = 0; i < bad_snp.size(); i++) {
			badfile << bad_snp[i] << "\t" << bad_allele1[i] << "\t" << bad_allele2[i] << "\t" << bad_refA[i] << endl;
		}
		badfile.close();
		cout << "Warning: a number of SNPs from the phenotype file could not be matched to the refrence GWAS data. These SNPs are saved in \"" + badname + "\"." << endl;

		vector<string>().swap(bad_snp);
		vector<string>().swap(bad_allele1);
		vector<string>().swap(bad_allele2);
		vector<string>().swap(bad_refA);
	}

	num_snps = bim_chr.size();
	ref_A = bim_allele1;
	other_A = bim_allele2;

	cout << "Number of SNPs read from .bim file: " << num_snps << "." << endl;

	sanitise_list();
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
	int i;

	to_include.clear(); // These to_includes here were commented out - can they remain like this??????
	to_include.resize(num_snps);
	snp_map.clear();

	for (i = 0; i < num_snps; i++) {
		stringstream ss;
		to_include[i] = i;

		if (snp_map.find(bim_snp_name[i]) != snp_map.end()) {
			ShowWarning("Duplicated SNP ID \"" + bim_snp_name[i] + "\".", a_verbose);

			ss << bim_snp_name[i] << "_" << i + 1;
			bim_snp_name[i] = ss.str();

			ShowWarning("This SNP has been changed to \"" + bim_snp_name[i] + "\".", a_verbose);
		}
		else
			ss << bim_snp_name[i];

		snp_map.insert(pair<string, int>(ss.str(), i));
	}
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
	string line;
	int count = 0;

	// Prepare for bim reading
	ifstream fam(famfile.c_str());
	if (!fam) {
		throw ("IO Error: fam file cannot be opened for reading: " + famfile + ".");
	}
	cout << "Reading data from fam file: " + famfile + "." << endl;
	fam_clear();

	// Iterate through bim file and save data
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
	int i = 0, size = 0;

	fam_ids_inc.clear();
	fam_ids_inc.resize(individuals);
	fam_map.clear();

	for (i = 0; i < individuals; i++) {
		fam_ids_inc[i] = i;
		fam_map.insert(pair<string, int>(fam_fid[i] + ":" + fam_iid[i], i));
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

/*
 * Calculates allele frequencies based on .fam data
 */
void reference::calculate_allele_freq()
{
	int i = 0, j = 0;
	const size_t fam_ids_size = fam_ids_inc.size(),
		bim_ids_size = to_include.size();

	cout << "Calculating allele frequencies from .fam data." << endl;
	mu.clear();
	mu.resize(num_snps);

#pragma omp parallel for
	for (i = 0; i < bim_ids_size; i++) {
		double fcount = 0.0, f_buf = 0.0;

		for (j = 0; j < fam_ids_size; j++) {
			if (!bed_snp_1[to_include[i]][fam_ids_inc[j]] || bed_snp_2[to_include[i]][fam_ids_inc[j]]) {
				f_buf = bed_snp_1[to_include[i]][fam_ids_inc[j]] + bed_snp_2[to_include[i]][fam_ids_inc[j]];
				if (bim_allele2[to_include[i]] == ref_A[to_include[i]])
					f_buf = 2.0 - f_buf;
				mu[to_include[i]] += f_buf;
				fcount += 1;
			}
		}

		if (fcount > 0.0)
			mu[to_include[i]] /= fcount;
	}
}

void reference::get_read_individuals(vector<int> &read_individuals)
{
	read_individuals.clear();
	read_individuals.resize(individuals);
	for (int i = 0; i < individuals; i++) {
		if (fam_map.find(fam_fid[i] + ":" + fam_iid[i]) != fam_map.end())
			read_individuals[i] = 1;
		else
			read_individuals[i] = 0;
	}
}

void reference::get_read_snps(vector<int> &read_snps)
{
	read_snps.clear();
	read_snps.resize(num_snps);
	for (int i = 0; i < num_snps; i++) {
		if (snp_map.find(bim_snp_name[i]) != snp_map.end())
			read_snps[i] = 1;
		else
			read_snps[i] = 0;
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
	int i, j, k,
		snp_idx, ind_idx;
	const size_t bim_size = to_include.size(),
		fam_size = fam_ids_inc.size();
	vector<int> read_individuals, read_snps;
	
	char ch[1];
	bitset<8> b;
	fstream BIT(bedfile.c_str(), ios::in | ios::binary);

	get_read_individuals(read_individuals);
	get_read_snps(read_snps);

	if (!bim_size || !fam_size) {
		ShowError("No data from either .bim or .fam files to continue with.");
		return;
	}

	if (!BIT) {
		ShowError("Cannot open " + bedfile + " to read.");
		return;
	}

	bed_snp_1.resize(bim_size);
	bed_snp_2.resize(bim_size);
	for (i = 0; i < bim_size; i++) {
		bed_snp_1[i].resize(fam_size);
		bed_snp_2[i].resize(fam_size);
	}

	for (i = 0; i < 3; i++) {
		BIT.read(ch, 1); // First three bytes are used for the file format and are not read
	}

	for (i = 0, snp_idx = 0; i < num_snps; i++) {
		// 00: homozygote AA
		// 11: homozygote BB
		// 01: heterozygote
		// 10: missing
		
		// SNP not found
		if (read_snps[i] == 0) {
			for (j = 0; j < individuals; j += 4)
				BIT.read(ch, 1);
			continue;
		}

		// SNP found - calculator allele frequency
		for (j = 0, ind_idx = 0; j < individuals;) {
			BIT.read(ch, 1);
			if (!BIT)
				throw ("Cannot read the .bed file?");

			b = ch[0];
			k = 0;
			while (k < 7 && j < individuals) { // 11 for AA; 00 for BB
				if (read_individuals[j] == 0)
					k += 2;
				else {
					bed_snp_2[snp_idx][ind_idx] = (!b[k++]);
					bed_snp_1[snp_idx][ind_idx] = (!b[k++]);
					ind_idx++;
				}
				j++;
			}
		}

		if (snp_idx == bim_size)
			break;
		snp_idx++;
	}

	BIT.clear();
	BIT.close();

	cout << "Genotype data for " << fam_size << " individuals and " << bim_size << " SNPs read." << endl;
}
