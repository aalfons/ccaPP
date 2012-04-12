/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#include <R.h>
#include "utils.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


// C++ isnan() is not portable and gives error on Windows systems
// use R macro ISNAN() instead

// --------------
// median and MAD
// --------------

// median using std::vector
// order of observations is messed up
double median(vector<double>& x) {
	int n = x.size();
	// find median
	int half = (n + 1) / 2;	// divide by integer to get truncated result
	half--;					// indices start with 0
	double med;
	if((n % 2) == 1) {
		// odd number of observations, take the middle one
		nth_element(x.begin(), x.begin()+half, x.end());
		med = x[half];
	} else {
		// even number of observations, take the mean of the two middle ones
		nth_element(x.begin(), x.begin()+half, x.end());
		double tmp = x[half];
		nth_element(x.begin(), x.begin()+half+1, x.end());
		med = 0.5 * (tmp + x[half+1]);
	}
	return med;
}

// median using arma::vec
double median(const vec& x) {
	uword n = x.n_elem;
	// make sure it returns NA in the presence of NA's
	for(uword i = 0; i < n; i++) {
		if(ISNAN(x(i))) return(NA_REAL);
	}
	// copy data to std::vector to not mess up the order of observations
	vector<double> xx(n);
	for(uword i = 0; i < n; i++) {
		xx[i] = x(i);
	}
	return median(xx);
}

// R interface to median()
SEXP R_fastMedian(SEXP R_x) {
	NumericVector Rcpp_x(R_x);						// convert data to Rcpp type
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// reuse memory
	return wrap(median(x));	// call arma version and wrap result
}

// MAD
// median is returned through the corresponding parameter
double mad(const vec& x, const double& constant, double& center) {
	uword n = x.n_elem;
	// make sure it returns NA in the presence of NA's
	for(uword i = 0; i < n; i++) {
		if(ISNAN(x(i))) return(NA_REAL);
	}
	// copy data to std::vector to not mess up the order of observations
	vector<double> xx(n);
	for(uword i = 0; i < n; i++) {
		xx[i] = x(i);
	}
	// find median
	center = median(xx);
	// compute MAD
	for(uword i = 0; i < n; i++) {
		xx[i] = abs(xx[i] - center);
	}
	return constant * median(xx);
}
double mad(const vec& x, double& center) {
	return mad(x, 1.4826, center);
}

// R interfaces to mad()
SEXP R_fastMAD(SEXP R_x, SEXP R_constant) {
	NumericVector Rcpp_x(R_x);						// convert data to Rcpp type
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// reuse memory
	double constant = as<double>(R_constant), center;
	double MAD = mad(x, constant, center);			// call arma version
	return List::create(
			Named("center") = wrap(center),
			Named("MAD") = wrap(MAD)
			);
}


// ---------------
// order and ranks
// ---------------

// define new data type for sorting
// .first component is index (type double since ties are broken by averaging)
// .second component is value
typedef pair<uword, double> sortData;

// unlike R, C++ returns false in comparisons with NA/NaN
bool sortDataLess(const sortData& left, const sortData& right) {
	return left.second < right.second;
}

// compute order of observations
// stable sorting is not necessary since ties are broken by averaging
uvec order(const vec& x) {
	// initialize data structure for sorting
	const uword n = x.n_elem;
	vector<sortData> foo(n);
	for(uword i = 0; i < n; i++) {
		foo[i] = sortData(i, x(i));
	}
	// call STL's sort()
	sort(foo.begin(), foo.end(), sortDataLess);
	// construct and return vector of indices
	uvec indices(n);
	for(uword i = 0; i < n; i++) {
		sortData bar = foo[i];
		indices(i) = bar.first;
	}
	return indices;
}

// compute ranks of observations in a vector
vec rank(const vec& x) {
	const uword n = x.n_elem;
	uword i, j, k;
	// compute order of observations
	uvec ord = order(x);
	// compute ranks (break ties by taking average)
	vec ranks(n);
	for(i = 0; i < n; i = j+1) {
		j = i;
		// check if there is a series of equal values
		while((j < n-1) && (x(ord(j)) == x(ord(j+1)))) {
			j++;
		}
		// break ties by average rank, otherwise this gives just the rank
		for(k = i; k <= j; k++) {
			ranks(ord(k)) = (i + j + 2) / 2.0;
		}
	}
	// return ranks
	return ranks;
}

// R interface to rank() (for testing)
SEXP R_rank(SEXP R_x) {
	NumericVector Rcpp_x(R_x);						// convert data to Rcpp type
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// convert data to arma type
	vec ranks = rank(x);							// call arma version
	return wrap(ranks.memptr(), ranks.memptr() + ranks.n_elem);
}
