/* @(#)hapfuse.hpp
 */

#ifndef _HAPFUSE_HPP
#define _HAPFUSE_HPP 1

#if __cplusplus <= 199711L
#error This library needs at least a C++11 compliant compiler
#endif

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
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <cfloat>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <future>
#include <memory>
#include <sys/stat.h>
#include <getopt.h>

#include "hfHelper.hpp"
#include "hapSamp.hpp"
#include "site.hpp"
#include "utils.hpp"
#include "writer.hpp"

class hapfuse {
private:
  const HapfuseHelper::init m_init;
  size_t m_numInputChunks;
  const bool m_out_GT = false;
  const bool m_out_GP = false;
  const bool m_out_APP = false;
  const bool m_in_GT = false;
  const bool m_in_GP = false;
  const bool m_in_APP = false;

  // reader
  std::vector<std::string> m_bcfFiles; // only used for bcfs
  std::vector<std::string> m_wtcccHapFiles;
  std::vector<std::string> m_wtcccSampFiles;

  std::vector<Site> load_chunk_WTCCC(const std::string &hapFile,
                                     const std::string &sampFile, bool first);
  std::vector<Site> load_chunk_bcf(const std::string &inFile, bool first);
  std::tuple<float, float> extractGP(float *gp, int gtA, int gtB);

  // writer
  hf::Writer m_writer;
  hf::WriterHelper::init m_writerInit;

  list<Site> site;
  HapfuseHelper::fileType m_inputFileType = HapfuseHelper::fileType::BCF;
  HapfuseHelper::fileType m_outputFileType = HapfuseHelper::fileType::BCF;
  std::vector<std::string> m_names;

  inline size_t numSamps() { return m_names.size(); }

  void find_overlap(std::vector<Site> chunk, std::list<Site>::iterator &first,
                    std::list<Site>::iterator &mid,
                    std::list<Site>::iterator &last, size_t &chunkMidIdx);
  void merge_chunk(std::vector<Site> chunk);

public:
  bool load_dir(const char *D);
  void work();
  std::vector<Site> load_chunk(size_t chunkIdx, bool first);

  // destroy output header and output file descriptor
  hapfuse(HapfuseHelper::init init);
};

#endif /* _HAPFUSE_HPP */
