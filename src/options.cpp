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
	spdlog::info(" ****************************************************************************************************");
	spdlog::info(" *                                             _..._       .-'''-.           _..._       .-'''-.    *");
	spdlog::info(" *                                          .-'_..._''.   '   _    \\      .-'_..._''.   '   _    \\  *");
	spdlog::info(" * _________   _...._                     .' .'      '.\\/   /` '.   \\   .' .'      '.\\/   /` '.   \\ *");
	spdlog::info(" * \\        |.'      '-.         _     _ / .'          .   |     \\  '  / .'          .   |     \\  ' *");
	spdlog::info(" *  \\        .'```'.    '. /\\    \\\\   //. '            |   '      |  '. '            |   '      |  '*");
	spdlog::info(" *   \\      |       \\     \\`\\\\  //\\\\ // | |            \\    \\     / / | |            \\    \\     / / *");
	spdlog::info(" *    |     |        |    |  \\`//  \\'/  | |             `.   ` ..' /  | |             `.   ` ..' /  *");
	spdlog::info(" *    |      \\      /    .    \\|   |/   . '                '-...-'`   . '                '-...-'`   *");
	spdlog::info(" *    |     |\\`'-.-'   .'      '         \\ '.          .               \\ '.          .              *");
	spdlog::info(" *    |     | '-....-'`                   '. `._____.-'/                '. `._____.-'/              *");
	spdlog::info(" *   .'     '.                              `-.______ /                   `-.______ /               *");
	spdlog::info(" * '-----------'                                     `                             `                *");
	spdlog::info(" ****************************************************************************************************");

	spdlog::info("!! Remember to periodically 'git pull' in your PWCoCo directory to obtain the latest updates and bug fixes !!");

	/*
	 * Option reading
	 */
	unsigned short chr = 0;
	int i = 0;
	double p_cutoff = 5e-8, collinear = 0.9, maf = 0.1, ld_window = 1.0e7,
		freq_threshold = 0.2, init_h4 = 80, top_snp = 1e10,
		p1 = 1e-4, p2 = 1e-4, p3 = 1e-5,
		n1 = 0.0, n2 = 0.0, n1_case = 0.0, n2_case = 0.0;
	string bfile = "", bim_file = "", fam_file = "", bed_file = "",
		phen1_file = "", phen2_file = "",
		out = "pwcoco_out", log = "pwcoco_log", snplist = "",
		opt;
	bool out_cond = false, cond_ssize = false,
		data_folder = false; // Whether the data is in folders or files

	for (i = 1; i < argc; i++) {
		opt = argv[i];

		if (opt == "--bfile") {
			bfile = argv[++i];

			if (bfile.substr(bfile.length() - 5) == ".bfile")
				bfile = bfile.substr(0, bfile.length() - 5);

			bim_file = bfile + ".bim";
			fam_file = bfile + ".fam";
			bed_file = bfile + ".bed";

			spdlog::info("Using {} as reference files.", bfile);
		}
		else if (opt == "--phen1_file" || opt == "--sum_stats1") {
			phen1_file = argv[++i];

			spdlog::info("Using {} summary statistic 1 files or folders.", phen1_file);
		}
		else if (opt == "--phen2_file" || opt == "--sum_stats2") {
			phen2_file = argv[++i];

			spdlog::info("Using {} summary statistic 2 files or folders.", phen2_file);
		}

		else if (opt == "--n1") {
			n1 = stod(argv[++i]);
			if (n1 <= 0) {
				n1 = 0;
			}

			spdlog::info("--n1 {}.", n1);
		}
		else if (opt == "--n2") {
			n2 = stod(argv[++i]);
			if (n2 <= 0) {
				n2 = 0;
			}

			spdlog::info("--n2 {}.", n2);
		}
		else if (opt == "--n1_case") {
			n1_case = stod(argv[++i]);
			if (n1_case <= 0) {
				n1_case = 0;
			}

			spdlog::info("--n1_case {}.", n1_case);
		}
		else if (opt == "--n2_case") {
			n2_case = stod(argv[++i]);
			if (n2_case <= 0) {
				n2_case = 0;
			}

			spdlog::info("--n2_case {}.", n2_case);
		}
		else if (opt == "--log") {
			log = argv[++i];

			spdlog::info("--log {}.", log);
		}
		else if (opt == "--p_cutoff") {
			p_cutoff = stod(argv[++i]);

			spdlog::info("--p_cutoff {}.", p_cutoff);
		}
		else if (opt == "--out") {
			out = argv[++i];

			spdlog::info("--out {}.", out);
		}
		else if (opt == "--chr") {
			chr = (unsigned short)stoi(argv[++i]);

			if (chr < 1 || chr > 23) {
				chr = (chr < 1 ? 1 : chr > 23 ? 23 : chr);
			}

			spdlog::info("--chr {}.", chr);
		}
		else if (opt == "--top_snp") {
			top_snp = stod(argv[++i]);

			if (top_snp < 1 || top_snp > 10000) {
				top_snp = (top_snp <= 1 ? 1 : top_snp >= 10000 ? 10000 : top_snp);
			}

			spdlog::info("--top_snp {}.", top_snp);
		}
		else if (opt == "--ld_window") {
			ld_window = stoi(argv[++i]);

			if (ld_window > 10000) {
				ld_window = 10000;
			}
			ld_window *= 1000;

			spdlog::info("--ld_window {}.", ld_window);
		}
		else if (opt == "--collinear") {
			collinear = stod(argv[++i]);

			if (collinear < 0.01 || collinear > 0.99) {
				collinear = (collinear <= 0.01 ? 0.01 : collinear >= 0.99 ? 0.99 : collinear);
			}

			spdlog::info("--collinear {}.", collinear);
		}
		else if (opt == "--maf") {
			maf = stod(argv[++i]);

			if (maf < 0.0 || maf > 0.5) {
				maf = (maf < 0.0 ? 0.0 : maf > 0.5 ? 0.5 : maf);
			}

			spdlog::info("--maf {}.", maf);
		}
		else if (opt == "--freq_threshold") {
			freq_threshold = stod(argv[++i]);

			if (freq_threshold < 0.0 || freq_threshold > 1.0) {
				freq_threshold = (freq_threshold < 0.0 ? 0.0 : freq_threshold > 1.0 ? 1.0 : freq_threshold);
			}

			spdlog::info("--freq_threshold {}.", freq_threshold);
		}
		else if (opt == "--init_h4") {
			init_h4 = stod(argv[++i]);

			if (init_h4 < 0.0 || init_h4 > 100.0) {
				init_h4 = (init_h4 < 0.0 ? 0.0 : init_h4 > 100.0 ? 100.0 : init_h4);
			}

			spdlog::info("--init_h4 {}.", init_h4);
		}
		else if (opt == "--out_cond") {
			out_cond = true;

			spdlog::info("--out_cond.");
		}
		else if (opt == "--coloc_pp") {
			p1 = stod(argv[++i]);
			p2 = stod(argv[++i]);
			p3 = stod(argv[++i]);

			p1 = (p1 > 1.0 ? 1.0 : p1 < 1e-50 ? 1e-50 : p1);
			p2 = (p2 > 1.0 ? 1.0 : p2 < 1e-50 ? 1e-50 : p2);
			p3 = (p3 > 1.0 ? 1.0 : p3 < 1e-50 ? 1e-50 : p3);

			spdlog::info("--coloc_pp p1 {} p2 {} p3 {}.", p1, p2, p3);
		}
		else if (opt == "--cond_ssize") {
			cond_ssize = true;

			spdlog::info("--cond_ssize.");
		}
	}

	// First set up the logger
	vector<spdlog::sink_ptr> sinks;
#ifdef _WIN32
	sinks.push_back(make_shared<spdlog::sinks::wincolor_stdout_sink_mt>());
#else
	sinks.push_back(make_shared<spdlog::sinks::ansicolor_stdout_sink_mt>());
#endif
	try {
		sinks.push_back(make_shared<spdlog::sinks::basic_file_sink_mt>(log + ".txt"));
		auto logger = make_shared<spdlog::logger>("pwcoco_log", begin(sinks), end(sinks)); 
		spdlog::set_default_logger(logger);
	}
	catch (const spdlog::spdlog_ex &ex) {
		cout << "Log setup failed: " << ex.what() << " - attempting to create default log \"pwcoco_log.txt\"." << endl;
		sinks.push_back(make_shared<spdlog::sinks::basic_file_sink_mt>("pwcoco_log.txt"));
		auto logger = make_shared<spdlog::logger>("pwcoco_log", begin(sinks), end(sinks));
		spdlog::set_default_logger(logger);
	}
	spdlog::flush_on(spdlog::level::info);

	// Some analytics why not?
	chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	// .bim file MUST be supplied
	if (bim_file.compare("") == 0) {
		spdlog::critical("No .bim file found; a .bim file MUST be supplied!");
		return 0;
	}
	//else if (!file_exists(bim_file)) {
	//	spdlog::critical(".bim file {} cannot be opened.", bim_file);
	//	return;
	//}

	// .fam file MUST be supplied
	if (fam_file.compare("") == 0) {
		spdlog::critical("No .fam file found; a .fam file MUST be supplied!");
		return 0;
	}
	//else if (!file_exists(fam_file)) {
	//	spdlog::critical(".fam file {} cannot be opened.", fam_file);
	//	return 0;
	//}

	// Check if summary stats 1 is file or directory
	if (phen1_file.compare("") == 0) {
		spdlog::critical("Summary stats 1 file has not been supplied, please supply this before continuing.");
		return 0;
	}
	else {
		const fs::path path(phen1_file);
		error_code ec;

		if (fs::is_directory(path, ec)) {
			data_folder = true;
			spdlog::info("Summary stats 1 is treated as a folder.");
		}
		else if (fs::is_regular_file(path, ec)) {
			data_folder = false;
			spdlog::info("Summary stats 1 is treated as a file.");
		}
	}

	// Check to see if summary stats 2 matches 1
	if (phen2_file.compare("") == 0) {
		spdlog::critical("Summary stats 2 file has not been supplied, please supply this before continuing.");
		return 0;
	}
	else {
		const fs::path path(phen2_file);
		error_code ec;

		if (data_folder && fs::is_directory(path, ec)) {
			spdlog::info("Summary stats 2 is treated as a folder.");
		}
		else if (data_folder && !fs::is_directory(path, ec)) {
			spdlog::info("Summary stats 2 expected to be a folder but is not. Please fix this before continuing.");
			return 0;
		}
		else if (!data_folder && fs::is_regular_file(path, ec)) {
			spdlog::info("Summary stats 2 is treated as a file.");
		}
		else if (!data_folder && !fs::is_regular_file(path, ec)) {
			spdlog::info("Summary stats 2 expected to be a file but is not. Please fix this before continuing.");
			return 0;
		}
	}

	// Set up for some common variables
	reference *ref = new reference(out, chr); // Reference dataset
	init_h4 /= 100; // coloc returns h4 as a decimal

	// Depending on whether the summary statistics are given as a folder
	// or as separate files, we either:
	// 1. Load the reference panel into memory and preserve it across analyses (folders)
	// 2. Load the necessary SNPs from the reference panel (files)
	// These steps should increase efficiency and speed for the two different cases
	if (data_folder) 
	{
		// Case 1
		// Folders were given
		// therefore we expect the same file name in two different folders
		using recursive_directory_iterator = fs::recursive_directory_iterator;
		for (const auto &dir_entry : recursive_directory_iterator(phen1_file)) 
		{
			string filename{ dir_entry.path().filename().u8string() },
				path_to_file1{ dir_entry.path().u8string() },
				path_to_file2 = phen2_file + "\\" + filename;

			phenotype *exposure = init_pheno(path_to_file1, "exposure", n1, n1_case);
			phenotype *outcome = init_pheno(path_to_file2, "outcome", n2, n2_case);
			if (exposure->has_failed() || outcome->has_failed()) {
				spdlog::error("Reading of either summary statistic files has failed; have these been moved or altered?");
				spdlog::error("File 1: {}", path_to_file1);
				spdlog::error("File 2: {}", path_to_file2);
				continue;
			}

			if (initial_coloc(exposure, outcome, out, p1, p2, p3, init_h4)) {
				continue;
			}

			if (!ref->is_ready()) {
				// Bim-related first
				if (ref->read_bimfile(bim_file) == 0) {
					return 0;
				}
				// In case 1, we load all of the SNPs in the reference data as they may be required
				ref->whole_bim();
				ref->sanitise_list();

				// Fam-related
				if (ref->read_famfile(fam_file) == 0) {
					return 0;
				}

				// Finally bed-related
				if (ref->read_bedfile(bed_file) == 0) {
					return 0;
				}
				//ref->calculate_allele_freq();
				if (maf > 0.0) {
					if (ref->filter_snp_maf(maf) == 0)
						return 0;
				}
			}

			// Do the related conditional and colocalisation analyses
			if (pwcoco_sub(exposure, outcome, ref, p_cutoff, collinear, ld_window, out, top_snp, freq_threshold, cond_ssize, out_cond, p1, p2, p3)) {
				continue;
			}
		}
	}
	else {
		// Case 2
		// Files were given
		phenotype *exposure = init_pheno(phen1_file, "exposure", n1, n1_case);
		phenotype *outcome = init_pheno(phen2_file, "outcome", n2, n2_case);
		if (exposure->has_failed() || outcome->has_failed()) {
			spdlog::critical("Reading of either summary statistic files has failed; have these been moved or altered?");
			return 1;
		}

		if (initial_coloc(exposure, outcome, out, p1, p2, p3, init_h4)) {
			return 0;
		}

		// Bim-related first
		if (ref->read_bimfile(bim_file) == 0) {
			return 0;
		}
		// In case 2, we only need those SNPs which have already been matched between the exposure and the outcome
		ref->match_bim(exposure->get_snp_names(), outcome->get_snp_names());
		ref->sanitise_list();

		// Fam-related
		if (ref->read_famfile(fam_file) == 0) {
			return 0;
		}

		// Finally bed-related
		if (ref->read_bedfile(bed_file) == 0) {
			return 0;
		}
		//ref->calculate_allele_freq();
		if (maf > 0.0) {
			if (ref->filter_snp_maf(maf) == 0)
				return 0;
		}

		// Do the related conditional and colocalisation analyses
		if (pwcoco_sub(exposure, outcome, ref, p_cutoff, collinear, ld_window, out, top_snp, freq_threshold, cond_ssize, out_cond, p1, p2, p3)) {
			return 0;
		}
	}

#ifdef PYTHON_INC
	Py_Finalize();
#endif

	chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	spdlog::info("Analysis finished. Computational time: {} secs", (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1000000.0);
	return 0;
}

/*
 * Common function to do the initial pass of colocalisation
 * @param phenotype exposure read exposure data
 * @param phenotype outcome read outcome data
 * @param double p1 value for coloc p1
 * @param double p2 value for coloc p2
 * @param double p3 value for coloc p3
 * @param double init_h4 H4 threshold
 * @ret int
 */
int initial_coloc(phenotype *exposure, phenotype *outcome, string out, double p1, double p2, double p3, double init_h4)
{
	// First pass through - match data and perform coloc
	mdata *matched = new mdata(exposure, outcome);
	coloc_analysis *initial_coloc = new coloc_analysis(matched, out, p1, p2, p3);
	initial_coloc->init_coloc();

	if (initial_coloc->pp_abf[H4] > init_h4) {
		spdlog::info("Stopping algorthim as H4 for initial colocalisation analysis is already at or above threshold ({}%).", init_h4 * 100);
		return 1;
	}

	delete(initial_coloc);
	return 0;
}

/*
 * Common function to run the subsequent conditional and colocalisation analyses
 */
int pwcoco_sub(phenotype *exposure, phenotype *outcome, reference *ref, double p_cutoff, double collinear, double ld_window, string out, double top_snp,
	double freq_threshold, double cond_ssize, bool out_cond, double p1, double p2, double p3)
{
	// Holder for the conditional matrices
	conditional_dat *exp_cdat = new conditional_dat();
	conditional_dat *out_cdat = new conditional_dat();

	// Find each independent SNPs for both exposure and outcome data
	cond_analysis *exp_analysis = new cond_analysis(p_cutoff, collinear, ld_window, out, top_snp, freq_threshold, "exposure", cond_ssize);
	exp_analysis->init_conditional(exposure, ref);
	exp_analysis->find_independent_snps(exp_cdat, ref);

	cond_analysis *out_analysis = new cond_analysis(p_cutoff, collinear, ld_window, out, top_snp, freq_threshold, "outcome", cond_ssize);
	out_analysis->init_conditional(outcome, ref);
	out_analysis->find_independent_snps(out_cdat, ref);

	// If both of the conditionals failed - don't continue
	if (exp_analysis->get_num_ind() == 0 && out_analysis->get_num_ind() == 0)
	{
		spdlog::warn("Both conditional analyses failed to run or find any conditionally independednt signals, and so no colocalisation will be run.");
		return 1;
	}

	spdlog::info("There are {} selected SNPs in the exposure dataset and {} in the outcome dataset.", exp_analysis->get_num_ind(), out_analysis->get_num_ind());
	spdlog::info("Performing {} conditional and colocalisation analyses.", (exp_analysis->get_num_ind() == 0 ? 1 : exp_analysis->get_num_ind()) * (out_analysis->get_num_ind() == 0 ? 1 : out_analysis->get_num_ind()));

#ifdef PYTHON_INC
	Py_Initialize();

	PyObject *sys = PyImport_ImportModule("sys");
	PyObject *path = PyObject_GetAttrString(sys, "path");
#endif

	// Perform PWCoCo!
	string exp_snp_name = "", out_snp_name = "";
	if (exp_analysis->get_num_ind() == 0)
		exp_snp_name = "unconditioned";
	if (out_analysis->get_num_ind() == 0)
		out_snp_name = "unconditioned";

//#pragma omp parallel for
	for (int i = 0; i < exp_analysis->get_num_ind() + 1; i++)
	{
		if (exp_snp_name != "unconditioned" && i < exp_analysis->get_num_ind()) {
			// New conditional data struct for parallelisation
			conditional_dat *par_exp_cdat = new conditional_dat(*exp_cdat);

			exp_analysis->pw_conditional(exp_analysis->get_num_ind() > 1 ? i : -1, out_cond, par_exp_cdat, ref); // Be careful not to remove the only independent SNP
			exp_snp_name = exp_analysis->get_ind_snp_name(i);

			delete(par_exp_cdat);
		}
		else if ((exp_snp_name != "unconditioned" && i >= exp_analysis->get_num_ind())
			|| (exp_snp_name == "unconditioned" && i > 0))
		{
			break;
		}

		for (int j = 0; j < out_analysis->get_num_ind() + 1; j++)
		{
			if (out_snp_name != "unconditioned" && j < out_analysis->get_num_ind()) {
				conditional_dat *par_out_cdat = new conditional_dat(*out_cdat);

				out_analysis->pw_conditional(out_analysis->get_num_ind() > 1 ? j : -1, out_cond, par_out_cdat, ref);
				out_snp_name = out_analysis->get_ind_snp_name(j);

				delete(par_out_cdat);
			}
			else if ((out_snp_name != "unconditioned" && j >= out_analysis->get_num_ind())
				|| (out_snp_name == "unconditioned" && j > 0))
			{
				break;
			}

			mdata *matched_conditional;
			if (exp_snp_name == "unconditioned")
				matched_conditional = new mdata(out_analysis, exposure);
			else if (out_snp_name == "unconditioned")
				matched_conditional = new mdata(exp_analysis, outcome);
			else
				matched_conditional = new mdata(exp_analysis, out_analysis);

			coloc_analysis *conditional_coloc = new coloc_analysis(matched_conditional, out, p1, p2, p3);
			conditional_coloc->init_coloc(exp_snp_name, out_snp_name);

			delete(matched_conditional);
			delete(conditional_coloc);
		}
	}

	return 0;
}
