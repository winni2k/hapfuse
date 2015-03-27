/* wkretzsch@gmail.com writer.hPP
 */

#ifndef _WRITER_HPP
#define _WRITER_HPP 1

#include <exception>
#include <list>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <htslib/hts.h>
#include <htslib/vcf.h>

#include "hfHelper.hpp"
#include "utils.hpp"
#include "site.hpp"

namespace hf {

namespace WriterHelper {

struct init {
  std::vector<std::string> outputFiles;
  std::string mode = "b";
  bool GT = false;
  bool GP = false;
  bool APP = false;
  std::string chrom;
  std::vector<std::string> sampNames;
  bool is_x = false;
  std::string genderFile;
};
}

class Writer {

private:
  WriterHelper::init m_init;
  HapfuseHelper::fileType m_outputFileType = HapfuseHelper::fileType::BCF;
  std::string m_outputBCFFile;
  std::string m_outputWTCCCHapsFile;
  std::string m_outputWTCCCSampleFile;
  set<std::string> m_maleIDs;

  htsFile *m_fusedVCF = nullptr;
  bcf_hdr_t *m_hdr_out = nullptr;
  ofile m_fusedWTCCCHaps;
  ofile m_fusedWTCCCSample;
  void write_bcf_site(const Site &osite, const std::vector<unsigned> &gts,
                      const std::vector<float> &lineGPs,
                      const std::vector<float> &lineAPPs) const;
  void write_wtccc_site(const Site &osite, const std::vector<unsigned> &gts);

  void write_head();
  void write_vcf_head();
  void write_wtccc_sample();

  void loadGenderFile();

public:
  Writer(){};
  ~Writer() {
    if (m_hdr_out != nullptr)
      bcf_hdr_destroy(m_hdr_out);
    if (m_fusedVCF != nullptr)
      hts_close(m_fusedVCF);
  }

  void init(hf::WriterHelper::init init);
  void write_site(const Site &osite);
  void write_sites(const std::list<Site> &outSites) {
    for (auto s : outSites)
      write_site(s);
  }
};
}
#endif /* _WRITER_HPP */
