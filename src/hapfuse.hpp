/* @(#)hapfuse.hpp
 */

#ifndef _HAPFUSE_HPP
#define _HAPFUSE_HPP 1

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <dirent.h>
#include <stdint.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <zlib.h>
#include <list>
#include <set>
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <cfloat>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <future>
#include <memory>
#include <sys/stat.h>
#include <getopt.h>

#include "site.hpp"
#include "hapSamp.hpp"
#include "utils.hpp"

#define EPSILON 0.001 // precision of input floats

namespace HapfuseHelper {
struct init {
  bool is_x = false;
  std::string outputFile = "";
  std::string mode = "v";
  bool useLinearWeighting = false;
  std::string wtcccHapFiles = "";
  std::string wtcccSampFiles = "";
  std::vector<std::string> inFiles;
};

void load_files_from_file(const std::string &fileFile,
                          std::vector<std::string> &inFiles);

// converts a probability to phred scale
double prob2Phred(double prob);

// converts a phred scaled probability back to a probability
double phred2Prob(double phred);

enum class fileType { WTCCC, BCF };
}

class hapfuse {
private:
  list<Site> site;
  std::vector<std::string> file; // only used for bcfs
  std::vector<std::string> m_wtcccHapFiles;
  std::vector<std::string> m_wtcccSampFiles;
  HapfuseHelper::fileType m_inputFileType = HapfuseHelper::fileType::BCF;
  set<std::string> male;
  std::vector<std::string> m_names;
  const HapfuseHelper::init m_init;

  inline size_t numSamps() { return m_names.size(); }

  void write_vcf_head();
  std::vector<Site> load_chunk(size_t chunkIdx, bool first);
  std::vector<Site> load_chunk_WTCCC(const std::string &hapFile,
                                     const std::string &sampFile, bool first);
  std::vector<Site> load_chunk_bcf(const std::string &inFile, bool first);

  std::tuple<float, float> extractGP(float *gp, int &gtA, int &gtB);
  void merge_chunk(std::vector<Site> chunk);

  htsFile *m_fusedVCF = NULL;
  bcf_hdr_t *m_hdr_out = NULL;
  void write_site(const Site &osite) const;
  void write_sites(const list<Site> &outSites) {
    for (auto s : outSites)
      write_site(s);
  }

public:
  static void document(void);
  bool is_x() { return m_init.is_x; }
  bool gender(const char *F);
  bool load_dir(const char *D);
  void work();

  // destroy output header and output file descriptor
  hapfuse(HapfuseHelper::init init);
  ~hapfuse() {
    bcf_hdr_destroy(m_hdr_out);
    hts_close(m_fusedVCF);
  }
};

#endif /* _HAPFUSE_HPP */
