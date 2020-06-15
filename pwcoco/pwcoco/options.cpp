/*
 * 
 * PWCOCO
 *
 *
 * Stolen from Jian Yang and Giambartolomei
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
	int i = 0, ld_window = 10000, top_snp = -1;
	double p_cutoff = 1e-3, collinear = 0.9;
	string bfile = "", bim_file = "", fam_file = "", bed_file = "",
		phen1_file = "", phen2_file = "",
		out = "";
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

		/* Optional */
		else if (_strcmpi(option_str[i], "--phen2_file") == 0) {
			phen2_file = option_str[++i];
			cout << "Phenotype 2 file is: " << phen2_file << endl;
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
				cout << "LD window was set too high, capping at 100,000\n" << endl;
				ld_window = 10000;
			}
			ld_window *= 1000;
		}
		else if (_strcmpi(option_str[i], "--collinear") == 0) {
			collinear = stod(option_str[++i]);
			if (collinear < 0.01 || collinear > 0.99) {
				cout << "Collinearity check should be in the range (0.01, 0.99), capping value given within this range.\n" << endl;
				collinear = (collinear <= 0.01 ? 0.01 : collinear >= 0.99 ? 0.99 : collinear);
			}
		}
	}
	cout << "Using P cutoff of " << p_cutoff << "." << endl;

	// .bim file MUST be supplied
	if (bim_file.compare("") == 0) {
		throw("No .bim file found; a .bim file MUST be supplied!");
	}
	else if (!file_exists(bim_file)) {
		throw("IO Error: .bim file \"" + bim_file + "\" cannot be opened.");
	}

	// .fam file MUST be supplied
	if (fam_file.compare("") == 0) {
		throw("No .fam file found; a .fam file MUST be supplied!");
	}
	else if (!file_exists(fam_file)) {
		throw("IO Error: .fam file \"" + fam_file + "\" cannot be opened.");
	}

	// Check phenotype 1 file before continuing
	if (phen1_file.compare("") == 0) {
		throw("Phenotype 1 file has not been supplied, please supply this before continuing.");
	}
	else if (!file_exists(phen1_file)) {
		throw("IO Error: Cannot open the phenotype 1 file: " + phen1_file + " for reading.");
	}

	// Check phenotype 2 file if given
	if (phen2_file.compare("") != 0 && !file_exists(phen2_file)) {
		throw("IO Error: Phenotype 2 file given but cannot open: " + phen2_file + " for reading.");
	}

	// Default output if not given
	if (out.compare("") == 0) {
		out = "results";
	}

	if (chr < 1 || chr > 23) {
		throw("Chromosome is out of bounds: " + chr);
	}
	else {
		cout << "Limiting analysis to chromosome " << chr << "." << endl;
	}

	// Initialise!
	cond_analysis *p_analysis = new cond_analysis(p_cutoff, collinear, ld_window, out, verbose, top_snp);
	reference *ref = new reference(out, chr, verbose);
	phenotype *exposure = init_exposure(phen1_file, "exposure");

	ref->read_bimfile(bim_file, exposure);
	ref->read_famfile(fam_file);
	ref->read_bedfile(bed_file);

	// Do the stuff!
	vector<int> selected, remain;
	eigenVector bC, bC_se, pC;
 	p_analysis->init_conditional(exposure, ref);
	p_analysis->stepwise_select(selected, remain, bC, bC_se, pC, ref);
	p_analysis->massoc(p_analysis);

	this_thread::sleep_for(5000ms);
	return;
}
