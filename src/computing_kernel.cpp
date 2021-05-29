#include "Rcpp.h"
#include "kernelforleon/computing_kernel.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix CPP_computing_kernel(const std::vector<std::string>& x,
                                         const unsigned long long& max_kmer_length,
                                         const double exponential) {
    return correlation_kernel_3(x, max_kmer_length, exponential);
}
