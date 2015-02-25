/* @(#)hapSamp.hpp
 */

#ifndef _HAPSAMP_HPP
#define _HAPSAMP_HPP 1

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "site.hpp"
#include "utils.hpp"
#include "hapfuse.hpp"

class HapSamp {

private:
  const std::string m_hapFile;
  ifile m_hapFD;
  std::vector<std::string> m_samps;

  std::string m_chrom;
  const bool m_multiSpaceAllowed = false;

  Site m_firstSite;
  unsigned m_linesRead = 0;
  std::string m_buffer; // unprocessed data

  void FillFirstSite();

public:
  HapSamp(std::string hapFile, const std::string &sampFile,
          bool multiSpaceAllowed = false);

  unsigned GetFirstPos();
  Site GetSite();
  void CheckSamps(const std::vector<std::string> &trueSamps);
  std::vector<std::string> GetSamps() { return m_samps; }
  std::string GetChrom();
};

namespace HapfuseHelper {

void load_files_from_file(const std::string &fileFile,
                          std::vector<std::string> &inFiles);

double prob2Phred(double prob);
double phred2Prob(double phred);
void wtccc_hap_order(std::vector<std::string> &wtccc_hap_files,
                     std::vector<std::string> &wtccc_samp_files);
}
#endif /* _HAPSAMP_HPP */
