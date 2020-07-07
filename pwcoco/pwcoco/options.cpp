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

	long int time_used = 0, start = time(NULL);
	time_t curr = time(0);
	char time_str[26];

	ctime_s(time_str, sizeof time_str, &curr);
	printf("Analysis started: %s\n", time_str);
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

	curr = time(0);
	ctime_s(time_str, sizeof time_str, &curr);
	printf("\nAnalysis finished: %s", time_str);
	time_used = time(NULL) - start;
	cout << "Computational time: " << time_used / 3600 << ":" << (time_used % 3600) / 60 << ":" << time_used % 60 << endl; // Change me to make me better please

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
	int i = 0, ld_window = 1e7, top_snp = -1;
	double p_cutoff = 5e-8, collinear = 0.9, maf = 0.0;
	string bfile = "", bim_file = "", fam_file = "", bed_file = "",
		phen1_file = "", phen2_file = "",
		out = "", snplist = "";
	bool verbose = true;

	/*
	 * Read through option flags given by user.
	 * For a list of these...
	 */
	for (i = 1; i < option_num; i++) {
		/* Required flags */
		if (_strcmpi(option_str[i], "--bfile") == 0) {
			bfile = option_str[++i];
			cout << "Using " << bfile << " as .bfile." << endl;

			if (bfile.substr(bfile.length() - 5) == ".bfile")
				bfile = bfile.substr(0, bfile.length() - 5);

			bim_file = bfile + ".bim";
			fam_file = bfile + ".fam";
			bed_file = bfile + ".bed";

		}
		else if (_strcmpi(option_str[i], "--phen1_file") == 0) {
			phen1_file = option_str[++i];
			cout << "Using " << phen1_file << " as phenotype 1 file." << endl;
		}
		else if (_strcmpi(option_str[i], "--phen2_file") == 0) {
			phen2_file = option_str[++i];
			cout << "Phenotype 2 file is: " << phen2_file << endl;
		}

		/* Optional */
		else if (_strcmpi(option_str[i], "--snp_list") == 0) {
			snplist = option_str[++i];
			cout << "Using " << snplist << " as list of SNPs upon which to condition." << endl;
		}
		else if (_strcmpi(option_str[i], "--p_cutoff") == 0) {
			p_cutoff = stod(option_str[++i]);
		}
		else if (_strcmpi(option_str[i], "--out") == 0) {
			out = option_str[++i];
		}
		else if (_strcmpi(option_str[i], "--verbose") == 0) {
			verbose = (bool)stoi(option_str[++i]);
		}
		else if (_strcmpi(option_str[i], "--chr") == 0) {
			chr = (unsigned short)stoi(option_str[++i]);
		}
		else if (_strcmpi(option_str[i], "--top_snp") == 0) {
			top_snp = stoi(option_str[++i]);
			if (top_snp < 1 || top_snp > 10000) {
				top_snp = (top_snp <= 1 ? 1 : top_snp >= 10000 ? 10000 : top_snp);
			}
		}
		else if (_strcmpi(option_str[i], "--ld_window") == 0) {
			ld_window = stoi(option_str[++i]);
			if (ld_window > 100000) {
				cout << "LD window was set too high, capping at 100,000." << endl;
				ld_window = 10000;
			}
			ld_window *= 1000;
		}
		else if (_strcmpi(option_str[i], "--collinear") == 0) {
			collinear = stod(option_str[++i]);
			if (collinear < 0.01 || collinear > 0.99) {
				cout << "Collinearity check should be in the range (0.01, 0.99), capping value given within this range." << endl;
				collinear = (collinear <= 0.01 ? 0.01 : collinear >= 0.99 ? 0.99 : collinear);
			}
		}
		else if (_strcmpi(option_str[i], "--maf") == 0) {
			maf = stod(option_str[++i]);
			if (maf < 0.0 || maf > 0.5) {
				cout << "MAF flag should be within the range of (0.0, 0.5). Clipping to be within this range." << endl;
				maf = (maf < 0.0 ? 0.0 : maf > 0.5 ? 0.5 : maf);
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
	phenotype *exposure = init_exposure(phen1_file, "exposure");
	phenotype *outcome = init_exposure(phen2_file, "outcome");

	mdata *matched = new mdata(exposure, outcome);
	coloc_analysis *initial_coloc = new coloc_analysis(matched, 1e-4, 1e-4, 1e-5);
	initial_coloc->init_coloc();

	if (initial_coloc->pp_abf[H4] >= 0.80) { // TODO user-specified flag
		cout << "Stopping algorthim as H4 for initial colocalisation analysis is already at or above 80%." << endl;
		//return;
	}

	// Initialise the reference
	reference *ref = new reference(out, chr, verbose);

	ref->read_bimfile(bim_file, exposure); // No need to do this for the outcome dataset as well as they have already been matched
	ref->read_famfile(fam_file);
	ref->read_bedfile(bed_file);
	ref->calculate_allele_freq();
	if (maf > 0.0)
		ref->filter_snp_maf(maf);

	// Find each independent SNPs for both exposure and outcome data
	cond_analysis *exp_analysis = new cond_analysis(p_cutoff, collinear, ld_window, out, verbose, top_snp);
	cond_analysis *out_analysis = new cond_analysis(p_cutoff, collinear, ld_window, out, verbose, top_snp);

	exp_analysis->init_conditional(exposure, ref);
	//out_analysis->init_conditional(outcome, ref);

	//cout << "There are " << exp_analysis->sw_snps << " selected SNPs in the exposure dataset and " << out_analysis->sw_snps << " in the outcome dataset." << endl;
	//cout << "Performing " << exp_analysis->sw_snps * out_analysis->sw_snps << " conditional  and colocalisation analyses." << endl;
	
	exp_analysis->massoc(ref, snplist);

	return;
}
