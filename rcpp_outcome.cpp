/*The C++ code defined a parallel processing routine using the Rcpp and RcppParallel 
libraries to efficiently search for specific ICD codes within large datasets. It started by 
importing the necessary libraries, Rcpp.h for seamless R and C++ integration and RcppParallel.h 
for parallel execution capabilities. The function returned the result matrix, which contained 
the indices and dates where the specified ICD code was found. This approach leveraged parallel 
computing to significantly speed up the processing of large datasets, making the search for ICD 
codes much more efficient.
by Luo wenjin
date 2024-01-12
*/

#include <Rcpp.h>
using namespace Rcpp;
#include <RcppParallel.h>
using namespace RcppParallel;

// Worker class to find rows containing the code and extract dates
struct FindDates : public Worker {
  // Input matrices and search code
  const CharacterMatrix mat; // Matrix of codes
  const CharacterMatrix dates; // Matrix of corresponding dates
  const std::string code; // Code to search for
  // Output matrix
  CharacterMatrix result; // Matrix to store the results
  
  // Constructor
  FindDates(const CharacterMatrix& mat, const CharacterMatrix& dates, const std::string& code, CharacterMatrix& result)
  : mat(mat), dates(dates), code(code), result(result) {}
  
  // Overload operator() for parallel execution
  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t i = begin; i < end; ++i) {
      for (int j = 0; j < mat.ncol(); ++j) {
        if (std::string(mat(i, j)).find(code) != std::string::npos) {
          result(i, 0) = std::to_string(i + 1); // Store row index (1-based)
          result(i, 1) = dates(i, j); // Store corresponding date
          break;
        }
      }
    }
  }
};

// [[Rcpp::export]]
CharacterMatrix findDatesWithCode(CharacterMatrix mat, CharacterMatrix dateMat, std::string code) {
  CharacterMatrix result(mat.nrow(), 2); // Initialize result matrix with 2 columns
  std::fill(result.begin(), result.end(), NA_STRING); // Fill with NA values
  
  FindDates worker(mat, dateMat, code, result); // Create worker instance
  parallelFor(0, mat.nrow(), worker); // Execute parallel processing
  
  return result; // Return the result matrix
}
