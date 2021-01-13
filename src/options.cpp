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

	try {
		option(argc, argv);
	}
	catch (const string &err_msg) {
		cerr << "\n" << err_msg << endl;
	}
	catch (const char* err_msg) {
		cerr << "\n" << err_msg << endl;
	}
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
	int i = 0;
	double p_cutoff = 5e-8, collinear = 0.9, maf = 0.1, ld_window = 1.0e7,
		freq_threshold = 0.2, init_h4 = 80, top_snp = 1e10;
	string bfile = "", bim_file = "", fam_file = "", bed_file = "",
		phen1_file = "", phen2_file = "",
		out = "pwcoco_out", log = "pwcoco_log", snplist = "",
		opt;
	bool out_cond = false;

	/*
	 * Read through option flags given by user.
	 * For a list of these...
	 */
	for (i = 1; i < option_num; i++) {
		opt = option_str[i];
		/* Required flags */
		if (opt == "--bfile") {
			bfile = option_str[++i];

			if (bfile.substr(bfile.length() - 5) == ".bfile")
				bfile = bfile.substr(0, bfile.length() - 5);

			bim_file = bfile + ".bim";
			fam_file = bfile + ".fam";
			bed_file = bfile + ".bed";

		}
		else if (opt == "--phen1_file") {
			phen1_file = option_str[++i];
		}
		else if (opt == "--phen2_file") {
			phen2_file = option_str[++i];
		}

		/* Optional */
		else if (opt == "--log") {
			log = option_str[++i];
		}
		else if (opt == "--p_cutoff") {
			p_cutoff = stod(option_str[++i]);
		}
		else if (opt == "--out") {
			out = option_str[++i];
		}
		else if (opt == "--chr") {
			chr = (unsigned short)stoi(option_str[++i]);

			if (chr < 1 || chr > 23) {
				chr = (chr < 1 ? 1 : chr > 23 ? 23 : chr);
			}
		}
		else if (opt == "--top_snp") {
			top_snp = stod(option_str[++i]);

			if (top_snp < 1 || top_snp > 10000) {
				top_snp = (top_snp <= 1 ? 1 : top_snp >= 10000 ? 10000 : top_snp);
			}
		}
		else if (opt == "--ld_window") {
			ld_window = stoi(option_str[++i]);

			if (ld_window > 10000) {
				ld_window = 10000;
			}
			ld_window *= 1000;
		}
		else if (opt == "--collinear") {
			collinear = stod(option_str[++i]);

			if (collinear < 0.01 || collinear > 0.99) {
				collinear = (collinear <= 0.01 ? 0.01 : collinear >= 0.99 ? 0.99 : collinear);
			}
		}
		else if (opt == "--maf") {
			maf = stod(option_str[++i]);

			if (maf < 0.0 || maf > 0.5) {
				maf = (maf < 0.0 ? 0.0 : maf > 0.5 ? 0.5 : maf);
			}
		}
		else if (opt == "--freq_threshold") {
			freq_threshold = stod(option_str[++i]);

			if (freq_threshold < 0.0 || freq_threshold > 1.0) {
				freq_threshold = (freq_threshold < 0.0 ? 0.0 : freq_threshold > 1.0 ? 1.0 : freq_threshold);
			}
		}
		else if (opt == "--init_h4") {
			init_h4 = stod(option_str[++i]);

			if (init_h4 < 0.0 || init_h4 > 100.0) {
				init_h4 = (init_h4 < 0.0 ? 0.0 : init_h4 > 100.0 ? 100.0 : init_h4);
			}
		}
		else if (opt == "--out_cond") {
			out_cond = stoi(option_str[++i]) ? true : false;
		}
	}

	// First set up the logger
	vector<spdlog::sink_ptr> sinks;
#ifdef _WIN32
	sinks.push_back(make_shared<spdlog::sinks::wincolor_stdout_sink_st>());
#else
	sinks.push_back(make_shared<spdlog::sinks::ansicolor_stdout_sink_st>());
#endif
	try {
		sinks.push_back(make_shared<spdlog::sinks::basic_file_sink_st>(log + ".txt"));
		auto logger = make_shared<spdlog::logger>("pwcoco_log", begin(sinks), end(sinks)); 
		spdlog::set_default_logger(logger);
	}
	catch (const spdlog::spdlog_ex &ex) {
		cout << "Log setup failed: " << ex.what() << " - attempting to create default log \"pwcoco_log.txt\"." << endl;
		sinks.push_back(make_shared<spdlog::sinks::basic_file_sink_st>("pwcoco_log.txt"));
		auto logger = make_shared<spdlog::logger>("pwcoco_log", begin(sinks), end(sinks));
		spdlog::set_default_logger(logger);
	}
	spdlog::flush_every(std::chrono::seconds(10));

	// Some analytics why not?
	chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	spdlog::info("Using {} as reference files (.bfile, etc.)", bfile);
	spdlog::info("Using {} and {} as phenotype files.", phen1_file, phen2_file);
	spdlog::info("Using P cutoff of {}", p_cutoff);

	// .bim file MUST be supplied
	if (bim_file.compare("") == 0) {
		spdlog::critical("No .bim file found; a .bim file MUST be supplied!");
		return;
	}
	else if (!file_exists(bim_file)) {
		spdlog::critical(".bim file {} cannot be opened.", bim_file);
		return;
	}

	// .fam file MUST be supplied
	if (fam_file.compare("") == 0) {
		spdlog::critical("No .fam file found; a .fam file MUST be supplied!");
		return;
	}
	else if (!file_exists(fam_file)) {
		spdlog::critical(".fam file {} cannot be opened.", fam_file);
		return;
	}

	// Check phenotype 1 file before continuing
	if (phen1_file.compare("") == 0) {
		spdlog::critical("Phenotype 1 file has not been supplied, please supply this before continuing.");
		return;
	}
	else if (!file_exists(phen1_file)) {
		spdlog::critical("Cannot open the phenotype 1 file: {} for reading.", phen1_file);
		return;
	}

	// Check phenotype 2 file if given
	if (phen2_file.compare("") != 0 && !file_exists(phen2_file)) {
		spdlog::critical("Phenotype 2 file given but cannot open: {} for reading.", phen2_file);
		return;
	}

	if (chr) {
		spdlog::info("Limiting analysis to chromosome {}.", chr);
	}
	// coloc returns h4 as a decimal
	init_h4 /= 100;

	// First pass through - match data and perform coloc
	phenotype *exposure = init_pheno(phen1_file, "exposure");
	phenotype *outcome = init_pheno(phen2_file, "outcome");
	if (exposure->has_failed() || outcome->has_failed()) {
		spdlog::critical("Reading of either phenotype file failed; have these been moved or altered?");
		return;
	}

	mdata *matched = new mdata(exposure, outcome);
	coloc_analysis *initial_coloc = new coloc_analysis(matched, out, 1e-4, 1e-4, 1e-5);
	initial_coloc->init_coloc();

	if (initial_coloc->pp_abf[H4] > init_h4) {
		spdlog::info("Stopping algorthim as H4 for initial colocalisation analysis is already at or above threshold ({}%).", init_h4*100);
		return;
	}
	delete(initial_coloc);

	// Initialise the reference
	reference *ref = new reference(out, chr);

	// Bim-related first
	if (ref->read_bimfile(bim_file) == 0) { // No need to do this for the outcome dataset as well as they have already been matched
		return;
	}
	ref->match_bim(exposure->get_snp_names(), outcome->get_snp_names());
	ref->sanitise_list();

	// Fam-related
	if (ref->read_famfile(fam_file) == 0) {
		return;
	}

	// Finally bed-related
	if (ref->read_bedfile(bed_file) == 0) {
		return;
	}
	ref->calculate_allele_freq();
	if (maf > 0.0) {
		if (ref->filter_snp_maf(maf) == 0)
			return;
	}

	// Find each independent SNPs for both exposure and outcome data
	cond_analysis *exp_analysis = new cond_analysis(p_cutoff, collinear, ld_window, out, top_snp, freq_threshold, "exposure");
	exp_analysis->init_conditional(exposure, ref);
	exp_analysis->find_independent_snps(ref);

	cond_analysis *out_analysis = new cond_analysis(p_cutoff, collinear, ld_window, out, top_snp, freq_threshold, "outcome");
	out_analysis->init_conditional(outcome, ref);
	out_analysis->find_independent_snps(ref);

	spdlog::info("There are {} selected SNPs in the exposure dataset and {} in the outcome dataset.", exp_analysis->get_num_ind(), out_analysis->get_num_ind());
	spdlog::info("Performing {} conditional and colocalisation analyses.", (exp_analysis->get_num_ind() == 0 ? 1 : exp_analysis->get_num_ind()) * (out_analysis->get_num_ind() == 0 ? 1 : out_analysis->get_num_ind()));
	
	// If both of the conditionals failed - don't continue
	if (exp_analysis->get_num_ind() == 0 && out_analysis->get_num_ind() == 0) 
	{
		spdlog::warn("Both conditional analyses failed to run and so no colocalisation will be run.");
		return;
	}

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
			exp_analysis->pw_conditional(exp_analysis->get_num_ind() > 1 ? (int)i : -1, out_cond, ref); // Be careful not to remove the only independent SNP
			exp_snp_name = exp_analysis->get_ind_snp_name(i);
		}
		else if ((exp_snp_name != "unconditioned" && i >= exp_analysis->get_num_ind())
			|| (exp_snp_name == "unconditioned" && i > 0))
		{
			break;
		}

		for (int j = 0; j < out_analysis->get_num_ind(); j++)
		{
			if (out_snp_name != "unconditioned" && i < out_analysis->get_num_ind()) {
				out_analysis->pw_conditional(out_analysis->get_num_ind() > 1 ? (int)j : -1, out_cond, ref);
				out_snp_name = out_analysis->get_ind_snp_name(j);
			}
			else if ((out_snp_name != "unconditioned" && i >= out_analysis->get_num_ind())
				|| (out_snp_name == "unconditioned" && i > 0))
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

			coloc_analysis *conditional_coloc = new coloc_analysis(matched_conditional, out, 1e-4, 1e-4, 1e-5);
			conditional_coloc->init_coloc(exp_snp_name, out_snp_name);

			delete(matched_conditional);
			delete(conditional_coloc);
		}
	}

#ifdef PYTHON_INC
	Py_Finalize();
#endif

	chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	spdlog::info("Analysis finished. Computational time: {} secs", (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1000000.0);
	return;
}
