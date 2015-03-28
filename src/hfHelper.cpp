
#include "hfHelper.hpp"

namespace HapfuseHelper {

void load_files_from_file(const std::string &fileFile,
                          std::vector<std::string> &inFiles) {

  ifile fdFiles(fileFile);
  std::string line;
  if (!fdFiles.isGood())
    throw std::runtime_error("Could not open file: " + fileFile);

  while (getline(fdFiles, line))
    inFiles.push_back(line);
}

// converts a probability to phred scale
double prob2Phred(double prob) {
  assert(prob >= 0 - EPSILON);

  if (prob < 0)
    prob = 0;

  if (prob == 0)
    return 999;
  else if (prob >= 1)
    return 2 * DBL_MIN; // 0 comes out negative zero on my system, not cool!
  else
    return -10 * log10(prob);
}

// converts a phred scaled probability back to a probability
double phred2Prob(double phred) {
  assert(phred >= 0);

  if (phred >= DBL_MAX_10_EXP * 10)
    return DBL_MIN;
  else
    return pow(10, -phred / 10);
}

void wtccc_hap_order(vector<string> &wtccc_hap_files,
                     vector<string> &wtccc_samp_files) {
  assert(wtccc_hap_files.size() == wtccc_samp_files.size());

  // fill vector with order and first position
  vector<pair<unsigned, size_t>> wtccc_first_pos;
  wtccc_first_pos.reserve(wtccc_hap_files.size());
  for (size_t fileNum = 0; fileNum < wtccc_hap_files.size(); ++fileNum) {
    HapSamp wtccc(wtccc_hap_files[fileNum], wtccc_samp_files[fileNum]);
    wtccc_first_pos.push_back(make_pair(wtccc.GetFirstPos(), fileNum));
  }

  order_by_first_pos(wtccc_first_pos, wtccc_hap_files);
  order_by_first_pos(wtccc_first_pos, wtccc_samp_files);
}

void bcf_order(vector<string> &bcf_files) {

  assert(!bcf_files.empty());

  std::unique_ptr<bcf1_t, void (*)(bcf1_t *)> rec(
      bcf_init1(), [](bcf1_t *b) { bcf_destroy1(b); });

  // fill vector with order and first position
  vector<pair<unsigned, size_t>> bcf_first_pos;
  bcf_first_pos.reserve(bcf_files.size());

  size_t fileNum = 0;
  for (auto file : bcf_files) {
    std::unique_ptr<htsFile, void (*)(htsFile *)> fp(
        hts_open(file.c_str(), "r"), [](htsFile *f) { hts_close(f); });
    if (!fp)
      throw myException("Could not open file: " + file);

    std::unique_ptr<bcf_hdr_t, void (*)(bcf_hdr_t *)> hdr(
        bcf_hdr_read(fp.get()), [](bcf_hdr_t *h) { bcf_hdr_destroy(h); });

    if (bcf_read1(fp.get(), hdr.get(), rec.get()) < 0)
      throw runtime_error("BCF file [" + file + "] seems to have no entries");

    // get the first five columns
    bcf_unpack(rec.get(), BCF_UN_STR);

    // site will store the site's information
    bcf_first_pos.push_back(make_pair(rec->pos + 1, fileNum));
    ++fileNum;
  }

  order_by_first_pos(bcf_first_pos, bcf_files);
}

// this returns a vector filled with the first n_max_tokens -1 tokens of str
// the last token contains the rest of str
string::size_type tokenize_partial(string &str, size_t n_max_tokens,
                                   vector<string> &tokens) {
  tokens.clear();
  tokens.reserve(n_max_tokens);
  string::size_type p_last = str.find_first_not_of(" \t", 0);
  string::size_type p_curr = str.find_first_of(" \t", p_last);
  if (n_max_tokens > 1)
    while ((string::npos != p_curr || string::npos != p_last) &&
           tokens.size() + 1 != n_max_tokens) {
      tokens.push_back(str.substr(p_last, p_curr - p_last));
      p_last = str.find_first_not_of(" \t", p_curr);
      p_curr = str.find_first_of(" \t", p_last);
    }

  // push in string as last token
  tokens.push_back(std::move(str));
  return p_last;
}

int tokenize_from(const string &str, string::size_type p_last,
                  vector<string> &tokens) {
  string new_str = str.substr(p_last, string::npos);
  return sutils::tokenize(new_str, tokens);
}
}
