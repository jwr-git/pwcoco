/*
 * 
 * PWCOCO
 *
 *
 * With major thanks to the GCTA-COJO and COLOC R package teams!
 */

#include "options.h"

int main(int argc, char* argv[])
{
	cout << " ****************************************************************************************************" << endl;
	cout << " *                                             _..._       .-'''-.           _..._       .-'''-.    *" << endl;
	cout << " *                                          .-'_..._''.   '   _    \\      .-'_..._''.   '   _    \\  *" << endl;
	cout << " * _________   _...._                     .' .'      '.\\/   /` '.   \\   .' .'      '.\\/   /` '.   \\ *" << endl;
	cout << " * \\        |.'      '-.         _     _ / .'          .   |     \\  '  / .'          .   |     \\  ' *" << endl;
	cout << " *  \\        .'```'.    '. /\\    \\\\   //. '            |   '      |  '. '            |   '      |  '*" << endl;
	cout << " *   \\      |       \\     \\`\\\\  //\\\\ // | |            \\    \\     / / | |            \\    \\     / / *" << endl;
	cout << " *    |     |        |    |  \\`//  \\'/  | |             `.   ` ..' /  | |             `.   ` ..' /  *" << endl;
	cout << " *    |      \\      /    .    \\|   |/   . '                '-...-'`   . '                '-...-'`   *" << endl;
	cout << " *    |     |\\`'-.-'   .'      '         \\ '.          .               \\ '.          .              *" << endl;
	cout << " *    |     | '-....-'`                   '. `._____.-'/                '. `._____.-'/              *" << endl;
	cout << " *   .'     '.                              `-.______ /                   `-.______ /               *" << endl;
	cout << " * '-----------'                                     `                             `                *" << endl;
	cout << " ****************************************************************************************************" << endl;

	chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	cout << "Options: " << endl;

	try {
		option(argc, argv);
	}
	catch (const string &err_msg) {
		cerr << "\n" << err_msg << endl;
	}
	catch (const char* err_msg) {
		cerr << "\n" << err_msg << endl;
	}

	cout << "Analysis finished.\nComputational time: " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1000000.0 << " (secs)" << endl; // Change me to make me better please

	return 0;
}

/*
 * Function which handles cmd line options flags
 * @param int option_num Number of arguments
 * @param char *option_str Arguments given
 * @ret
 */
void option(int option_num, char* option_str[])
{
	unsigned short chr = 0;
	int i = 0, top_snp = -1;
	double p_cutoff = 5e-8, collinear = 0.9, maf = 0.0, ld_window = 1.0e7,
		freq_threshold = 0.2;
	string bfile = "", bim_file = "", fam_file = "", bed_file = "",
		phen1_file = "", phen2_file = "",
		out = "", snplist = "",
		opt;
	bool verbose = true;

	/*
	 * Read through option flags given by user.
	 * For a list of these...
	 */
	for (i = 1; i < option_num; i++) {
		opt = option_str[i];
		/* Required flags */
		if (opt == "--bfile") {
			bfile = option_str[++i];
			cout << "Using " << bfile << " as .bfile." << endl;

			if (bfile.substr(bfile.length() - 5) == ".bfile")
				bfile = bfile.substr(0, bfile.length() - 5);

			bim_file = bfile + ".bim";
			fam_file = bfile + ".fam";
			bed_file = bfile + ".bed";

		}
		else if (opt == "--phen1_file") {
			phen1_file = option_str[++i];
			cout << "Using " << phen1_file << " as phenotype 1 file." << endl;
		}
		else if (opt == "--phen2_file") {
			phen2_file = option_str[++i];
			cout << "Phenotype 2 file is: " << phen2_file << endl;
		}

		/* Optional */
		else if (opt == "--snp_list") {
			snplist = option_str[++i];
			cout << "Using " << snplist << " as list of SNPs upon which to condition." << endl;
		}
		else if (opt == "--p_cutoff") {
			p_cutoff = stod(option_str[++i]);
		}
		else if (opt == "--out") {
			out = option_str[++i];
		}
		else if (opt == "--verbose") {
			verbose = true;
		}
		else if (opt == "--chr") {
			chr = (unsigned short)stoi(option_str[++i]);
		}
		else if (opt == "--top_snp") {
			top_snp = stoi(option_str[++i]);
			if (top_snp < 1 || top_snp > 10000) {
				top_snp = (top_snp <= 1 ? 1 : top_snp >= 10000 ? 10000 : top_snp);
			}
		}
		else if (opt == "--ld_window") {
			ld_window = stoi(option_str[++i]);
			if (ld_window > 100000) {
				cout << "LD window was set too high, capping at 100,000." << endl;
				ld_window = 10000;
			}
			ld_window *= 1000;
		}
		else if (opt == "--collinear") {
			collinear = stod(option_str[++i]);
			if (collinear < 0.01 || collinear > 0.99) {
				cout << "Collinearity check should be in the range (0.01, 0.99), capping value given within this range." << endl;
				collinear = (collinear <= 0.01 ? 0.01 : collinear >= 0.99 ? 0.99 : collinear);
			}
		}
		else if (opt == "--maf") {
			maf = stod(option_str[++i]);
			if (maf < 0.0 || maf > 0.5) {
				cout << "MAF flag should be within the range of (0.0, 0.5). Clipping to be within this range." << endl;
				maf = (maf < 0.0 ? 0.0 : maf > 0.5 ? 0.5 : maf);
			}
		}
		else if (opt == "--freq_threshold") {
			freq_threshold = stod(option_str[++i]);
			if (freq_threshold < 0.0 || freq_threshold > 1.0) {
				cout << "Allele frequency threshold is not within the range of (0.0, 1.0). Clipping to be within this range." << endl;
				freq_threshold = (freq_threshold < 0.0 ? 0.0 : freq_threshold > 1.0 ? 1.0 : freq_threshold);
			}
		}
	}
	cout << "Using P cutoff of " << p_cutoff << "." << endl;

	// .bim file MUST be supplied
	if (bim_file.compare("") == 0) {
		ShowError("No .bim file found; a .bim file MUST be supplied!");
	}
	else if (!file_exists(bim_file)) {
		ShowError("IO Error: .bim file \"" + bim_file + "\" cannot be opened.");
	}

	// .fam file MUST be supplied
	if (fam_file.compare("") == 0) {
		ShowError("No .fam file found; a .fam file MUST be supplied!");
	}
	else if (!file_exists(fam_file)) {
		ShowError("IO Error: .fam file \"" + fam_file + "\" cannot be opened.");
	}

	// Check phenotype 1 file before continuing
	if (phen1_file.compare("") == 0) {
		ShowError("Phenotype 1 file has not been supplied, please supply this before continuing.");
	}
	else if (!file_exists(phen1_file)) {
		ShowError("IO Error: Cannot open the phenotype 1 file: " + phen1_file + " for reading.");
	}

	// Check snplist before continuing
	if (snplist.compare("") != 0 && !file_exists(snplist)) {
		ShowError("IO Error: Cannot open the snplist file: " + snplist + " for reading.");
	}

	// Check phenotype 2 file if given
	if (phen2_file.compare("") != 0 && !file_exists(phen2_file)) {
		ShowError("IO Error: Phenotype 2 file given but cannot open: " + phen2_file + " for reading.");
	}

	// Default output if not given
	if (out.compare("") == 0) {
		out = "results";
	}

	if (chr < 1 || chr > 23) {
		ShowError("Chromosome is out of bounds: " + chr);
	}
	else {
		cout << "Limiting analysis to chromosome " << chr << "." << endl;
	}

	// First pass through - match data and perform coloc
	phenotype *exposure = init_pheno(phen1_file, "exposure");
	phenotype *outcome = init_pheno(phen2_file, "outcome");

	mdata *matched = new mdata(exposure, outcome);
	coloc_analysis *initial_coloc = new coloc_analysis(matched, 1e-4, 1e-4, 1e-5);
	initial_coloc->init_coloc();

	if (initial_coloc->pp_abf[H4] >= 0.80) { // TODO user-specified flag
		cout << "Stopping algorthim as H4 for initial colocalisation analysis is already at or above 80%." << endl;
		//return;
	}

	// Initialise the reference
	reference *ref = new reference(out, chr, verbose);

	ref->read_bimfile(bim_file); // No need to do this for the outcome dataset as well as they have already been matched
	ref->match_bim(exposure->get_snp_names(), outcome->get_snp_names());
	ref->sanitise_list();
	ref->read_famfile(fam_file);
	ref->read_bedfile(bed_file);
	ref->calculate_allele_freq();
	if (maf > 0.0)
		ref->filter_snp_maf(maf);

	// Find each independent SNPs for both exposure and outcome data
	cond_analysis *exp_analysis = new cond_analysis(p_cutoff, collinear, ld_window, out, verbose, top_snp, freq_threshold, "exposure");
	exp_analysis->init_conditional(exposure, ref);
	exp_analysis->find_independent_snps(ref);

	cond_analysis *out_analysis = new cond_analysis(p_cutoff, collinear, ld_window, out, verbose, top_snp, freq_threshold, "outcome");
	out_analysis->init_conditional(outcome, ref);
	out_analysis->find_independent_snps(ref);

	cout << "There are " << exp_analysis->get_num_ind() << " selected SNPs in the exposure dataset and " << out_analysis->get_num_ind() << " in the outcome dataset." << endl;
	cout << "Performing " << (exp_analysis->get_num_ind() == 0 ? 1 : exp_analysis->get_num_ind()) * (out_analysis->get_num_ind() == 0 ? 1 : out_analysis->get_num_ind()) << " conditional and colocalisation analyses." << endl;
	
	// If both of the conditionals failed - don't continue
	if (exp_analysis->get_num_ind() == 0 && out_analysis->get_num_ind() == 0) 
	{
		cout << "Both conditional analyses failed to run and so no colocalisation will be run." << endl;
		return;
	}

	// Perform PWCoCo!
	// TODO There's a clearer way to do this for sure
	string exp_snp_name = "", out_snp_name = "";
	if (exp_analysis->get_num_ind() == 0) 
	{
		exp_snp_name = "[Unconditioned]";
		for (int j = 0; j < out_analysis->get_num_ind(); j++)
		{
			out_analysis->pw_conditional(out_analysis->get_num_ind() > 1 ? (int)j : -1, ref); // Be careful not to remove the only independent SNP
			out_snp_name = out_analysis->get_ind_snp_name(j);

			mdata *matched_conditional = new mdata(exp_analysis, out_analysis);
			coloc_analysis *conditional_coloc = new coloc_analysis(matched_conditional, 1e-4, 1e-4, 1e-5);
			initial_coloc->init_coloc(exp_snp_name, out_snp_name);

			delete(matched_conditional);
			delete(conditional_coloc);
		}
	}
	else if (out_analysis->get_num_ind() == 0)
	{
		out_snp_name = "[Unconditioned]";
		for (int i = 0; i < exp_analysis->get_num_ind(); i++)
		{
			exp_analysis->pw_conditional(exp_analysis->get_num_ind() > 1 ? (int)i : -1, ref); // Be careful not to remove the only independent SNP
			exp_snp_name = exp_analysis->get_ind_snp_name(i);

			mdata *matched_conditional = new mdata(exp_analysis, out_analysis);
			coloc_analysis *conditional_coloc = new coloc_analysis(matched_conditional, 1e-4, 1e-4, 1e-5);
			initial_coloc->init_coloc(exp_snp_name, out_snp_name);

			delete(matched_conditional);
			delete(conditional_coloc);
		}
	}
	else {
		for (int i = 0; i < exp_analysis->get_num_ind(); i++)
		{
			exp_analysis->pw_conditional(exp_analysis->get_num_ind() > 1 ? (int)i : -1, ref); // Be careful not to remove the only independent SNP
			exp_snp_name = exp_analysis->get_ind_snp_name(i);
			for (int j = 0; j < out_analysis->get_num_ind(); j++)
			{
				out_analysis->pw_conditional(out_analysis->get_num_ind() > 1 ? (int)j : -1, ref);
				out_snp_name = out_analysis->get_ind_snp_name(j);

				mdata *matched_conditional = new mdata(exp_analysis, out_analysis);
				coloc_analysis *conditional_coloc = new coloc_analysis(matched_conditional, 1e-4, 1e-4, 1e-5);
				initial_coloc->init_coloc(exp_snp_name, out_snp_name);

				delete(matched_conditional);
				delete(conditional_coloc);
			}
		}
	}

	return;
}
