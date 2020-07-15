#include "helper_funcs.h"

/*
 * Determines wheter a file exists.
 * @param const string name File path and name
 * @ret bool True if exists, false if not
 */
bool file_exists(const std::string &name)
{
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

/*
 * Helper function that shows warnings from the tool.
 * These warnings will not end the program.
 * @param const string msg Warning message
 * @param bool verbose Whether to show the warning or not.
 * @ret void
 */
void ShowWarning(const std::string msg, bool verbose)
{
	if (verbose) {
		std::cout << "\033[0;31m" << msg << "\033[0m" << std::endl;
	}
}

/*
 * Helper function that shows errors from the tool.
 * These warnings will end the program.
 * @param const string msg Error message
 * @ret void
 */
void ShowError(const std::string msg)
{
	throw(msg);
}

double absi(const double &x)
{
	std::complex<double> cld(x);
	double ldAbs = abs(cld);
	return ldAbs;
}

bool isFloatEqual(double lhs, double rhs)
{
	return absi(lhs - rhs) < FLOATERR;
}

double pchisq(double x, double df)
{
	double p, q;
	int st = 0; // Error variable
	int w = 1; // Function variable
	double bnd = 1; // Boundary function

	if (x < 0)
		return -9;

	// NCP is set to 0
	cdfchi(&w, &p, &q, &x, &df, &st, &bnd);

	if (st != 0)
		return -9;

	return q;
}

/*
 * Vector function that will find the median of the given vector.
 * @param vector<double> &x Vector whose median is required
 * @ret double Median value
 */
double v_calc_median(const std::vector<double> &x)
{
	std::vector<double> b(x);
	size_t size = b.size();
	if (size == 1)
		return b[0];
	std::stable_sort(b.begin(), b.end());
	if (size % 2 == 0)
		return (b[size / 2 - 1] + b[size / 2]) / 2;
	else
		return b[(size - 1) / 2];
}

std::vector<std::size_t> v_sort_indices(const std::vector<std::string> &v)
{
	// Initialise vector of original indices
	std::vector<std::size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// Sort indices based on comparing values in v
	std::sort(idx.begin(), idx.end(),
		[&v](std::size_t i1, std::size_t i2) { return v[i1] < v[i2]; });

	return idx;
}

void eigenVector2Vector(Eigen::VectorXd &x, std::vector<double> &y)
{
	y.resize(x.size());
	for (int i = 0; i < x.size(); i++)
		y[i] = x[i];
}

/*
 * Implementation of logsum from mvc library in R
 * TODO This differs slightly from the R result - need to check this out
 */
double logsum(const std::vector<double> &x)
{
	std::vector<double> temp = x;
	double sum = 0.0;
	auto m_it = *std::max_element(temp.begin(), temp.end());
	for (auto& it : temp) {
		sum += exp(it - m_it);
	}
	return m_it + log(sum);
}

double logdiff(double x, double y)
{
	double m = std::max(x, y);
	return m + log(exp(x - m) - exp(y - m));
}

/*
template<typename KeyType, typename LeftValue, typename RightValue>
std::map<KeyType, pair<LeftValue, RightValue> > IntersectMaps(const std::map<KeyType, LeftValue> &left, const std::map<KeyType, RightValue> &right)
{
	std::map<KeyType, pair<LeftValue, RightValue> > result;
	typename std::map<KeyType, LeftValue>::const_iterator il = left.begin();
	typename std::map<KeyType, RightValue>::const_iterator ir = right.begin();
	while (il != left.end() && ir != right.end())
	{
		if (il->first < ir->first)
			++il;
		else if (ir->first < il->first)
			++ir;
		else
		{
			result.insert(make_pair(il->first, make_pair(il->second, ir->second)));
			++il;
			++ir;
		}
	}
	return result;
}

std::map<std::string, int> vm_intersect(const std::map<std::string, int> &left, const std::vector<std::string> &right)
{
	std::map<std::string, int> result;
	typename std::map<std::string, int>::const_iterator il = left.begin();
	typename std::vector<std::string>::const_iterator ir = right.begin();

	while (il != left.end() && ir != right.end())
	{
		if (il->first < ir->c_str())
			++il;
		else if (ir->c_str() < il->first)
			++ir;
		else
		{
			result.insert(make_pair(il->first, il->second));
			++il;
			++ir;
		}
	}
	return result;
}

std::vector<int> vm_intersect(const std::map<std::string, int> &m, const std::vector<std::string> &v)
{
	std::vector<int> result;
	std::vector<std::string>::iterator iter = v.begin();

	// Check if string vector is in map and store the position
	while (iter != v.end())
	{
		if (m.count(iter.c_str()))
			result.push_back(m[iter.c_str()]);
		iter++;
	}
	return result;
}
*/

std::vector<double> lm(const std::vector<double> &x, const std::vector<double> &y)
{
	if (x.size() != y.size()) {
		ShowError("Cannot compute regression for colocalisation due to differing sizes of vectors.");
	}
	std::size_t n = x.size();
	double sumX = std::accumulate(x.begin(), x.end(), 0.0);
	double sumY = std::accumulate(y.begin(), y.end(), 0.0);
	double meanX = sumX / n;
	double meanY = sumY / n;
	std::vector<double> a = x;
	std::vector<double> b = y;
	transform(a.begin(), a.end(), a.begin(), [=](double r) { return r - meanX; });
	transform(b.begin(), b.end(), b.begin(), [=](double r) { return r - meanY; });

	double sdXX = 0.0, sdXY = 0.0;
	for (int i = 0; i < n; i++) {
		sdXX += a[i] * a[i];
		sdXY += a[i] * b[i];
	}

	double m = sdXY / sdXX;
	double c = meanY - (m * meanX);
	
	std::vector<double> res = x;
	transform(res.begin(), res.end(), res.begin(), [=](double r) { return m * r + c; });
	return res;

}

std::vector<double> lm_fixed(const std::vector<double> &x, const std::vector<double> &y)
{
	if (x.size() != y.size()) {
		throw("Cannot compute regression for colocalisation due to differing sizes of vectors.");
	}
	size_t i, n = x.size();
	double sumXY = 0.0;
	double sumXX = 0.0;
	for (size_t i = 0; i < n; i++) {
		sumXY += x[i] * y[i];
		sumXX += x[i] * x[i];
	}
	double m = sumXY / sumXX;

	std::vector<double> res = x;
	transform(res.begin(), res.end(), res.begin(), [=](double r) { return m * r; });
	return res;
}

std::string string2upper(const std::string &str)
{
	std::string t = str;
	transform(t.begin(), t.end(), t.begin(), ::toupper);
	return t;
}
