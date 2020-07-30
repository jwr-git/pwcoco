#pragma once

#include <algorithm>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <time.h>
#include <thread>

#include "data.h"
#include "conditional.h"
#include "coloc.h"
#include "helper_funcs.h"

using namespace std;

void option(int option_num, char* option_str[]);
