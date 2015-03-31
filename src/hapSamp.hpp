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
#include "hfHelper.hpp"

class HapSamp {

private:
  const std::string m_hapFile;
  ifile m_hapFD;
  std::vector<std::string> m_samps;

  const bool m_multiSpaceAllowed = false;

  // chromosome
  std::string m_chrom;
  bool m_checkChrom = true;

  Site m_firstSite;
  unsigned m_linesRead = 0;
  std::string m_buffer; // unprocessed data

  void FillFirstSite();

public:
  HapSamp(std::string hapFile, const std::string &sampFile,
          bool multiSpaceAllowed = false, bool checkChrom = true,
          std::string chrom = "");

  unsigned GetFirstPos();
  Site GetSite();
  void CheckSamps(const std::vector<std::string> &trueSamps);
  std::vector<std::string> GetSamps() { return m_samps; }
  std::string GetChrom();
};

#endif /* _HAPSAMP_HPP */
