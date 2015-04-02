/* @(#)hfHelper.hpp
 */

#ifndef _HFHELPER_HPP
#define _HFHELPER_HPP 1

#include <algorithm>
#include <unordered_map>
#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include "hapSamp.hpp"
#include "utils.hpp"

#define EPSILON 0.001 // precision of input floats

namespace HapfuseHelper {
enum class WeightingStyle { AVERAGE, LINEAR, STEP };
enum class fileType { WTCCC, BCF };

struct init {
  int verbosity = 0;
  std::string genderFile;
  bool is_x = false;
  std::vector<std::string> outputFiles;
  std::string mode = "b";
  WeightingStyle ws = WeightingStyle::STEP;
  std::string wtcccHapFilesFile;
  std::string wtcccSampFilesFile;
  std::vector<std::string> cmdLineInputFiles;
  std::unordered_map<std::string, bool> in_format_tags{
      {"GT", false}, {"GP", false}, {"APP", false}};
  std::unordered_map<std::string, bool> out_format_tags{
      {"GT", false}, {"GP", false}, {"APP", false}};
  bool unmatchedSitesOK = false;
  std::string alignMapFile;
  std::string assumeChrom;
};

void load_files_from_file(const std::string &fileFile,
                          std::vector<std::string> &inFiles);

// converts a probability to phred scale
double prob2Phred(double prob);

// converts a phred scaled probability back to a probability
double phred2Prob(double phred);

void wtccc_hap_order(std::vector<std::string> &wtccc_hap_files,
                     std::vector<std::string> &wtccc_samp_files);

void bcf_order(std::vector<std::string> &bcf_files);

template <class T>
void
order_by_first_pos(std::vector<std::pair<unsigned, std::size_t>> &first_pos,
                   std::vector<T> &toOrder) {

  sort(first_pos.begin(), first_pos.end());

  // sort hap files
  std::vector<T> tmp;
  tmp.reserve(toOrder.size());
  for (auto o : first_pos)
    tmp.push_back(std::move(toOrder[o.second]));
  std::swap(tmp, toOrder);
}

std::string::size_type tokenize_partial(std::string &str,
                                        std::size_t n_max_tokens,
                                        std::vector<std::string> &tokens);
int tokenize_from(const std::string &str, std::string::size_type p_last,
                  std::vector<std::string> &tokens);
}

#endif /* _HFHELPER_HPP */
