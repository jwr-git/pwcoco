#pragma once

#include <algorithm>
#include <complex>
#include <cstring>
#include <iostream>
#include <fcntl.h>
#include <fstream>
#include <map>
#include <numeric>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <sys/stat.h>
#include "dcdflib.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

#define FLOATERR std::numeric_limits<double>::epsilon();

bool file_exists(const std::string &name);
void checkEntry(std::string txt, double *val);

bool isFloatEqual(double lhs, double rhs);
std::string string2upper(const std::string &str);
double pchisq(double x, double df);

double v_calc_median(const std::vector<double> &x);
std::vector<std::size_t> v_sort_indices(const std::vector<std::string> &v);
std::vector<std::size_t> v_remove_nans(std::vector<int> &v);
std::vector<std::size_t> v_remove_nans(std::vector<size_t> &v);
std::vector<std::string> v_merge_nodupes(std::vector<std::string> v1, std::vector<std::string> v2);
void v_remove_dupes(std::vector<std::string> &v);
void v_remove_dupes(std::vector<size_t> &v);
void eigenVector2Vector(Eigen::VectorXd &x, std::vector<double> &y);
double logsum(const std::vector<double> &x);
double logdiff(double x, double y);

double lm(const std::vector<double> &x, const std::vector<double> &y);
double lm_fixed(const std::vector<double> &x, const std::vector<double> &y);

bool file_is_empty(std::ifstream &pFile);
bool isNumber(std::string s);

//template<typename KeyType, typename LeftValue, typename RightValue>
//std::map<KeyType, std::pair<LeftValue, RightValue> > IntersectMaps(const std::map<KeyType, LeftValue> &left, const std::map<KeyType, RightValue> &right);

//std::map<std::string, int> vm_intersect(const std::map<std::string, int> &left, const std::vector<std::string> &right);
