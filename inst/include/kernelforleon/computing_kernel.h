#pragma once

#include <cmath>
#include "Rcpp.h"
#include "BLOSUM62_2.h"

inline double kernel_1(const char &codon_1, const char &codon_2, const double &exponential = 0.1) {
    return pow(internal::read_BLOSUM62_2(codon_1, codon_2), exponential);
}

inline double kernels_2_3(const std::string &sequence_1,
                          const std::string &sequence_2,
                          const unsigned long long &max_kmer_length,
                          const double &exponential = 0.1) {
    typedef std::vector<std::string>::size_type LenSq;

    const LenSq sequence_1_size = sequence_1.size();
    const LenSq sequence_2_size = sequence_2.size();

    // Constructs a matrix of kernel^2_1 scores, but rotated by 45 degrees.
    // The reason is that k-mers for k >= 2 can be built based on 1-mers, but the shift between sequences must be preserved for all k-mer elements.
    // Such matrix ensures that each vector has all possible k-mers for given shift.
    // E.g. one of the vectors in the middle contains the following pairs: [(1,1), (2,2), (3,3), (4,4), ...],
    // while the next contains: [(1,2), (2,3), (3,4), ...].
    std::vector<std::vector<double>> kernel_2_1(sequence_1_size + sequence_2_size - 1);

    // Index keeps track of which vector is being filled
    LenSq index = 0;
    // Filling the first half of the matrix
    for (; index < sequence_1_size; ++index) {
        const LenSq sequence_1_start_index = sequence_1_size - index - 1;
        const LenSq sequence_2_start_index = 0;
        kernel_2_1[index] = std::vector<double>(std::min(index + 1, sequence_2_size));

        for (LenSq j = 0; j < kernel_2_1[index].size(); ++j) {
            kernel_2_1[index][j] = kernel_1(
                    sequence_1[sequence_1_start_index + j],
                    sequence_2[sequence_2_start_index + j],
                    exponential);
        }
    }

    // Filling the second half
    for (; index < kernel_2_1.size(); ++index) {
        const LenSq sequence_1_start_index = 0;
        // We start from (0, 1) here, because (0, 0) was already computed above
        const LenSq sequence_2_start_index = index - sequence_1_size + 1;
        kernel_2_1[index] = std::vector<double>(std::min(sequence_1_size + sequence_2_size - index - 1, sequence_1_size));

        for (LenSq j = 0; j < kernel_2_1[index].size(); ++j) {
            kernel_2_1[index][j] = kernel_1(
                    sequence_1[sequence_1_start_index + j],
                    sequence_2[sequence_2_start_index + j],
                    exponential);
        }
    }

    double sum = 0.0;
    index = 0;
    for (; index < kernel_2_1.size(); ++index) {
        for (LenSq i = 0; i < kernel_2_1[index].size(); ++i) {
            double kmer_product = 1.0;
            for (LenSq j = 0; j < max_kmer_length; ++j) {
                if (i + j >= kernel_2_1[index].size()) {
                    break;
                }
                kmer_product *= kernel_2_1[index][i + j];
                sum += kmer_product;
            }
        }
    }

    return sum;
}

inline double correlation_kernel_3(const std::string &sequence_1,
                                   const std::string &sequence_2,
                                   const double &self_similarity_1,
                                   const double &self_similarity_2,
                                   const unsigned long long &max_kmer_length,
                                   const double &exponential = 0.1) {
    return kernels_2_3(sequence_1, sequence_2, max_kmer_length, exponential) /
        sqrt(self_similarity_1 * self_similarity_2);
}

inline Rcpp::NumericMatrix correlation_kernel_3(const std::vector<std::string> &sq,
                                                             const unsigned long long &max_kmer_length,
                                                             const double &exponential = 0.1) {
    typedef std::vector<std::string>::size_type LenSq;

    // Initializes kernel_3 scores for similarity of sequences to themselves so that they aren't computed for each standardization
    std::vector<double> self_similarity(sq.size(), 0);
    for (LenSq i = 0; i < sq.size(); ++i) {
        self_similarity[i] = kernels_2_3(sq[i], sq[i], max_kmer_length, exponential);
    }

    // Initializes returned correlation matrix with 1s on diagonal, because correlation of a sequence to itself is equal to 1
    Rcpp::NumericMatrix ret = Rcpp::NumericMatrix::diag(sq.size(), 1);
    for (LenSq i = 0; i < sq.size(); ++i) {
        for (LenSq j = i + 1; j < sq.size(); ++j) {
            const double correlation_score = correlation_kernel_3(
                    sq[i], sq[j], self_similarity[i], self_similarity[j], max_kmer_length, exponential);
            ret(i, j) = correlation_score;
            ret(j, i) = correlation_score;
        }
    }
    return ret;
}
