/* @(#)hapSamp.cpp
 */

#include "hapSamp.hpp"

using namespace std;
using namespace sutils;

HapSamp::HapSamp(string hapFile, const string &sampFile)
    : m_hapFile(std::move(hapFile)), m_hapFD(m_hapFile) {

  if (!m_hapFD.isGood())
    throw std::runtime_error("Could not open file [" + m_hapFile + "]");

  // might as well read in samples at this time
  ifile sampFD(sampFile);
  if (!sampFD.isGood())
    throw std::runtime_error("Could not open file [" + sampFile + "]");
  string line;

  // discard header
  getline(sampFD, line);
  getline(sampFD, line);
  vector<string> tokens;
  while (getline(sampFD, line)) {
    tokens.clear();
    tokenize(line, tokens, 1);
    m_samps.push_back(std::move(tokens[0]));
  }
  sampFD.close();
}

unsigned HapSamp::GetFirstPos() {
  FillFirstSite();
  return m_firstSite.pos;
}

void HapSamp::FillFirstSite() {

  if (!m_firstSite.empty())
    return;

  // otherwise this is the first line, let's load it and find the first position
  assert(m_linesRead == 0);
  if (!m_hapFD.isGood())
    throw std::runtime_error("Could not open file [" + m_hapFile + "]");

  assert(m_buffer.empty());
  getline(m_hapFD, m_buffer);

  if (m_buffer.empty())
    throw std::runtime_error("First line file [" + m_hapFile + "] is empty");

  vector<string> tokens;
  tokenize(m_buffer, tokens, 5);

  vector<string> alls;
  alls.push_back(std::move(tokens.at(3)));
  alls.push_back(std::move(tokens.at(4)));
  try {
    m_firstSite.init(std::move(tokens.at(0)), stoul(tokens.at(2)),
                     std::move(alls));
  }
  catch (std::invalid_argument &e) {
    throw std::runtime_error("At line " + to_string(m_linesRead) +
                             "\nCould not read position of line: " + m_buffer);
  }
}

// check if site is empty() for eof
Site HapSamp::GetSite() {

  if (!m_hapFD.isGood())
    throw std::runtime_error("Could not open file [" + m_hapFile + "]");

  vector<string> tokens;
  if (!m_buffer.empty()) {
    tokenize(m_buffer, tokens);
    m_buffer.clear();
  } else {
    string line;
    if(!getline(m_hapFD, line))
        return Site();
    tokenize(line, tokens);
  }

  if (tokens.size() != 5 + 2 * m_samps.size())
    throw std::runtime_error("Input line " + to_string(m_linesRead + 1) +
                             " has wrong number of columns: " +
                             to_string(tokens.size()));

  if (!m_firstSite.empty() && tokens[0] != m_firstSite.chr)
    throw std::runtime_error("Encountered line with chromosome [" + tokens[0] +
                             "], which differs from first line [" +
                             m_firstSite.chr + "]");

  Site outSite;
  outSite.hap.reserve(m_samps.size() * 2);
  vector<string> alls;
  alls.push_back(std::move(tokens[3]));
  alls.push_back(std::move(tokens[4]));
  try {
    outSite.init(std::move(tokens[0]), stoul(tokens[2]), std::move(alls));
  }
  catch (std::invalid_argument &e) {
    throw std::runtime_error("At line " + to_string(m_linesRead) +
                             "\nCould not read position: " + tokens[2]);
  }

  for (size_t colNum = 5; colNum < tokens.size(); ++colNum) {
    if (tokens[colNum] == "0")
      outSite.hap.push_back(0.0);
    else if (tokens[colNum] == "1")
      outSite.hap.push_back(1.0);
    else
      throw std::runtime_error("Encountered unexpected allele (" +
                               tokens[colNum] + ") at line " +
                               to_string(m_linesRead + 1));
  }

  ++m_linesRead;
  return outSite;
}

void HapSamp::CheckSamps(const vector<string> &samps) {

  if (samps.size() != m_samps.size())
    throw std::runtime_error("Chunk has wrong number of samples");

  for (size_t sIdx = 0; sIdx < m_samps.size(); ++sIdx)
    if (samps.at(sIdx) != m_samps.at(sIdx))
      throw std::runtime_error("sample mismatch found. " + samps[sIdx] +
                               " != " + m_samps[sIdx]);
}

string HapSamp::GetChrom() {
  FillFirstSite();
  return m_firstSite.chr;
}
