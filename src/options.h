#pragma once

#include <algorithm>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <time.h>
#include <thread>
#include <filesystem>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

#include "data.h"
#include "conditional.h"
#include "coloc.h"
#include "helper_funcs.h"

using namespace std;
namespace fs = std::filesystem;

int initial_coloc(phenotype *exposure, phenotype *outcome, string out, double p1, double p2, double p3, double init_h4);
int pwcoco_sub(phenotype *exposure, phenotype *outcome, reference *ref, double p_cutoff1, double p_cutoff2, double collinear, double ld_window, string out, double top_snp,
	double freq_threshold, double cond_ssize, bool out_cond, double p1, double p2, double p3);
