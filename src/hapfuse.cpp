/**

License
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the 'Software'), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

@file hapfuse.cpp
@brief joint chunked haplotypes into chromosome wide haplotypes
@author Yi Wang
@date 04/01/2011

Later changes were made by Warren Kretzschmar <wkretzsch@gmail.com>

*/

#include "hapfuse.hpp"

using namespace std;

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

  sort(wtccc_first_pos.begin(), wtccc_first_pos.end());

  // sort hap files
  vector<string> tmp_wtccc_files;
  tmp_wtccc_files.reserve(wtccc_hap_files.size());
  for (auto o : wtccc_first_pos)
    tmp_wtccc_files.push_back(std::move(wtccc_hap_files[o.second]));
  swap(tmp_wtccc_files, wtccc_hap_files);

  // do same for samples
  tmp_wtccc_files.clear();
  for (auto o : wtccc_first_pos)
    tmp_wtccc_files.push_back(std::move(wtccc_samp_files[o.second]));
  swap(tmp_wtccc_files, wtccc_samp_files);
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

hapfuse::hapfuse(HapfuseHelper::init init)
    : m_init(std::move(init)), m_bcfFiles(m_init.cmdLineInputFiles) {

  // Input files will be WTCCC style
  if (!m_init.wtcccHapFilesFile.empty()) {
    m_inputFileType = HapfuseHelper::fileType::WTCCC;

    HapfuseHelper::load_files_from_file(m_init.wtcccHapFilesFile,
                                        m_wtcccHapFiles);
    if (m_init.wtcccSampFilesFile.empty())
      throw std::runtime_error(
          "Need to define WTCCC sample files list with -H if using -h");
    else
      HapfuseHelper::load_files_from_file(m_init.wtcccSampFilesFile,
                                          m_wtcccSampFiles);
    if (m_wtcccHapFiles.size() != m_wtcccSampFiles.size())
      throw std::runtime_error(
          "-h and -H need to contain same number of files");

    m_numInputChunks = m_wtcccHapFiles.size();

    // order input hap files by first position
    HapfuseHelper::wtccc_hap_order(m_wtcccHapFiles, m_wtcccSampFiles);
  }

  // else load BCFs from command line
  else {
    m_inputFileType = HapfuseHelper::fileType::BCF;
    if (m_bcfFiles.empty())
      throw std::runtime_error("No input files given");

    m_numInputChunks = m_bcfFiles.size();

    // order input BCFs by first position
    HapfuseHelper::bcf_order(m_bcfFiles);
  }

  // if multiple input files are given then error out
  if (!m_bcfFiles.empty() && !m_wtcccHapFiles.empty())
    throw std::runtime_error("BCF and WTCCC style input files specified. "
                             "Please only specify one input file type.");

  if (m_numInputChunks == 0)
    throw std::runtime_error("No input files specified");

  struct stat buffer;
  for (auto oneFile : m_bcfFiles)
    if (stat(oneFile.c_str(), &buffer) != 0)
      throw runtime_error("Input file does not exist [" + oneFile + "]");

  for (auto oneFile : m_wtcccHapFiles)
    if (stat(oneFile.c_str(), &buffer) != 0)
      throw runtime_error("Input file does not exist [" + oneFile + "]");

  for (auto oneFile : m_wtcccSampFiles)
    if (stat(oneFile.c_str(), &buffer) != 0)
      throw runtime_error("Input file does not exist [" + oneFile + "]");

  // open output file for writing
  assert(!m_fusedVCF);
  if (m_init.outputFile.empty())
    throw runtime_error("Please specify an output file");

  size_t found = m_init.mode.find_first_of("buzv");
  if (m_init.mode.size() > 0 &&
      (m_init.mode.size() != 1 || found == string::npos))
    throw runtime_error("-O needs to be 'b' for compressed bcf, 'u' for "
                        "uncompressed bcf, 'z' for compressed vcf or 'v' for "
                        "uncompressed vcf");
  string cmode = "w" + m_init.mode;
  m_fusedVCF = hts_open(m_init.outputFile.c_str(), cmode.c_str());
}

bool hapfuse::gender(const char *F) {
  ifstream fi(F);

  if (!fi) {
    cerr << "fail to open " << F << endl;
    return false;
  }

  string s, g;

  for (fi >> s >> g; !fi.eof(); fi >> s >> g)
    if (g == "male")
      male.insert(s);

  fi.close();
  return true;
}

// extract the GP fields, convert them to allelic probabilities and store in
// pHap1,2
std::tuple<float, float> hapfuse::extractGP(float *gp, int &gtA, int &gtB) {

  // convert GPs to probabilities
  vector<float> GPs;
  GPs.reserve(3);
  double sum = 0;
  for (int i = 0; i < 3; ++i, ++gp) {
    GPs.push_back(HapfuseHelper::phred2Prob(*gp));
    sum += GPs[i];
  }

  float pHap1 = 0, pHap2 = 0;

  if (fabs(sum - 1) >= EPSILON)
    throw std::runtime_error("GPs don't sum to 1: " + to_string(GPs[0]) + " " +
                             to_string(GPs[1]) + " " + to_string(GPs[2]) +
                             "; sum = " + to_string(sum));
  pHap1 = GPs[2];
  pHap2 = GPs[2];

  // swap alleles if evidence exists in GT field
  if (gtA != gtB) {
    if (gtB == 0)
      pHap1 += GPs[1];
    else
      pHap2 += GPs[1];
  } else {
    pHap1 += GPs[1] / 2;
    pHap2 += GPs[1] / 2;
  }
  return make_tuple(pHap1, pHap2);
}

vector<Site> hapfuse::load_chunk(size_t chunkIdx, bool first) {

  if (m_inputFileType == HapfuseHelper::fileType::BCF) {
    return load_chunk_bcf(m_bcfFiles.at(chunkIdx), first);
  } else if (m_inputFileType == HapfuseHelper::fileType::WTCCC) {
    return load_chunk_WTCCC(m_wtcccHapFiles.at(chunkIdx),
                            m_wtcccSampFiles.at(chunkIdx), first);
  } else
    throw std::logic_error("Unknown input chunk type");
}

vector<Site> hapfuse::load_chunk_WTCCC(const string &hapFile,
                                       const string &sampFile, bool first) {
  // load and check samples
  HapSamp chunk(std::move(hapFile), sampFile);

  // Open output VCF for writing
  if (first) {
    assert(!m_hdr_out);
    m_hdr_out = bcf_hdr_init("w");

    assert(m_names.empty());
    m_names = chunk.GetSamps();

    // populate header with sample names
    for (auto sampName : m_names)
      bcf_hdr_add_sample(m_hdr_out, sampName.c_str());

    bcf_hdr_add_sample(m_hdr_out, NULL);

    // add contig
    bcf_hdr_printf(m_hdr_out, "##contig=<ID=%s>", chunk.GetChrom().c_str());

    write_vcf_head();
  } else
    chunk.CheckSamps(m_names);

  // load haplotypes
  vector<Site> sites;
  while (true) {
    Site line = chunk.GetSite();
    line.weight = 1;
    if (line.empty())
      break;
    sites.push_back(std::move(line));
  }
  return sites;
}

vector<Site> hapfuse::load_chunk_bcf(const string &inFile, bool first) {

  //  auto del = [](htsFile *f) { hts_close(f); };
  //  std::unique_ptr<htsFile, decltype(del)> fp(hts_open(inFile.c_str(), "r"),
  //                                             del);
  std::unique_ptr<htsFile, void (*)(htsFile *)> fp(
      hts_open(inFile.c_str(), "r"), [](htsFile *f) { hts_close(f); });
  if (!fp)
    throw myException("Could not open file: " + inFile);

  std::unique_ptr<bcf_hdr_t, void (*)(bcf_hdr_t *)> hdr(
      bcf_hdr_read(fp.get()), [](bcf_hdr_t *h) { bcf_hdr_destroy(h); });

  std::unique_ptr<bcf1_t, void (*)(bcf1_t *)> rec(
      bcf_init1(), [](bcf1_t *b) { bcf_destroy1(b); });

  // save first header as template for output file
  if (first) {
    assert(!m_hdr_out);
    m_hdr_out = bcf_hdr_dup(hdr.get());

    // save names to internal string vector
    assert(m_names.size() == 0);
    m_names.reserve(bcf_hdr_nsamples(m_hdr_out));
    for (int i = 0; i < bcf_hdr_nsamples(m_hdr_out); ++i)
      m_names.push_back(m_hdr_out->samples[i]);

    write_vcf_head();
  }

  vector<string> name;
  vector<Site> chunk;

  // skip headers

  // parse #CHROM header line
  if (bcf_hdr_nsamples(hdr.get()) != static_cast<int>(m_names.size()))
    throw runtime_error("Unequal number of samples in vcf files: " +
                        to_string(bcf_hdr_nsamples(hdr.get())) + " != " +
                        to_string(m_names.size()));
  for (size_t i = 0; i < m_names.size(); ++i) {
    string name(hdr->samples[i]);
    if (name != m_names[i])
      throw runtime_error("Sample names do not match: " + name + " != " +
                          m_names[i]);
  }
  unsigned in = m_names.size();

  // parsing each line of data
  unsigned cnt_lines = 0;
  while (bcf_read1(fp.get(), hdr.get(), rec.get()) >= 0) {

    // get the first five and sample columns
    bcf_unpack(rec.get(), BCF_UN_STR | BCF_UN_IND);

    // site will store the site's information
    Site site;
    site.hap.resize(in * 2);
    site.weight = 1; // still need to weight correctly
    string a, b;

    ++cnt_lines;

    // tokenize on "\t"
    // fill in site information (s)
    //    vector<string> tokens;
    //    boost::split(tokens, buffer, boost::is_any_of("\t"));
    site.chr = bcf_hdr_id2name(hdr.get(), rec->rid);
    site.pos = (rec->pos + 1);

    // make sure site is biallelic
    assert(rec->n_allele == 2);
    string a1(rec->d.allele[0]);
    string a2(rec->d.allele[1]);

    // make sure site is snp
    site.all.push_back(a1);
    site.all.push_back(a2);

    if (!bcf_get_fmt(hdr.get(), rec.get(), "GT"))
      throw std::runtime_error("expected GT field in VCF");
    if (!bcf_get_fmt(hdr.get(), rec.get(), "GP") &&
        !bcf_get_fmt(hdr.get(), rec.get(), "APP"))
      throw std::runtime_error("expected GP or APP field");

    // read sample specific data
    // get genotype array
    int ngt, *gt_arr = NULL, ngt_arr = 0;
    ngt = bcf_get_genotypes(hdr.get(), rec.get(), &gt_arr, &ngt_arr);

    if (ngt != 2 * bcf_hdr_nsamples(hdr.get()))
      throw std::runtime_error("number of genotypes found = " + to_string(ngt) +
                               "; with values: " + to_string(gt_arr[0]) + " " +
                               to_string(gt_arr[1]));

    // get APP array
    int n_arr = 0, m_arr = 0;
    float *arr = NULL;
    int stride = 0;
    if (bcf_get_fmt(hdr.get(), rec.get(), "APP")) {
      stride = 2;
      n_arr = bcf_get_format_float(hdr.get(), rec.get(), "APP", &arr, &m_arr);
      assert(n_arr / stride == bcf_hdr_nsamples(hdr.get()));
    }
    // get GP array
    else if (bcf_get_fmt(hdr.get(), rec.get(), "GP")) {
      stride = 3;
      n_arr = bcf_get_format_float(hdr.get(), rec.get(), "GP", &arr, &m_arr);
    } else {
      free(arr);
      throw std::runtime_error("could not read record");
    }

    // cycle through sample vals
    for (int sampNum = 0; sampNum < bcf_hdr_nsamples(hdr.get()); ++sampNum) {

      assert(gt_arr[sampNum * 2] != bcf_gt_missing);
      assert(gt_arr[sampNum * 2] != bcf_int32_vector_end);
      assert(gt_arr[sampNum * 2 + 1] != bcf_gt_missing);
      assert(gt_arr[sampNum * 2 + 1] != bcf_int32_vector_end);

      // this assumes the second gt val carries the phase information
      if (!bcf_gt_is_phased(gt_arr[sampNum * 2 + 1]))
        throw std::runtime_error("Error in GT data, genotype is not phased.");

      //    for (uint tokenColIdx = firstSampIdx; tokenColIdx < tokens.size();
      //         ++tokenColIdx) {

      // parse haps

      int gtA = bcf_gt_allele(gt_arr[sampNum * 2]);
      int gtB = bcf_gt_allele(gt_arr[sampNum * 2 + 1]);

      // extract/estimate allelic probabilities

      float pHap1;
      float pHap2;
      // parse APPs
      if (bcf_get_fmt(hdr.get(), rec.get(), "APP")) {
        pHap1 = HapfuseHelper::phred2Prob(arr[sampNum * stride]);
        pHap2 = HapfuseHelper::phred2Prob(arr[sampNum * stride + 1]);
      }
      // parse GPs
      else if (bcf_get_fmt(hdr.get(), rec.get(), "GP")) {
        std::tie(pHap1, pHap2) = extractGP(&arr[sampNum * stride], gtA, gtB);
      } else
        throw std::runtime_error("could not load GP or APP field ");

      assert(pHap1 + pHap2 < 2 + EPSILON);

      // make sure pHap1 and 2 are greater than zero
      if (pHap1 < 0)
        pHap1 = 0;

      if (pHap2 < 0)
        pHap2 = 0;

      site.hap[sampNum * 2] = pHap1;
      site.hap[sampNum * 2 + 1] = pHap2;
    }
    chunk.push_back(std::move(site));
    free(arr);
    free(gt_arr);
  }

  // calculate correct weights
  if (m_init.ws == HapfuseHelper::WeightingStyle::LINEAR) {
    assert(!chunk.empty());
    unsigned chunkStartPos = chunk.front().pos;
    unsigned chunkEndPos = chunk.back().pos;
    for (auto &cSite : chunk) {
      cSite.weight =
          min(cSite.pos - chunkStartPos, chunkEndPos - cSite.pos) + 1;

      // adjust haps according to weights
      for (auto &h : cSite.hap)
        h = h * cSite.weight;
    }
  }

  return chunk;
}

void hapfuse::write_vcf_head() {

  assert(m_hdr_out);
  assert(m_fusedVCF);

  if (m_init.out_format_tags.at("GT") == true)
    if (bcf_hdr_append(
            m_hdr_out,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"))
      throw runtime_error("could not append VCF GT format header");

  if (m_init.out_format_tags.at("APP") == true)
    if (bcf_hdr_append(m_hdr_out,
                       "##FORMAT=<ID=APP,Number=2,Type=Float,Description="
                       "\"Phred-scaled allelic probability, "
                       "P(Allele=1|Haplotype)\">"))
      throw runtime_error("Could not append VCF APP format header");

  if (m_init.out_format_tags.at("GP") == true)
    if (bcf_hdr_append(
            m_hdr_out,
            "##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Phred-scaled "
            "genotype posterior probabilities\">"))
      throw runtime_error("could not append VCF GP format header");

  ostringstream version;
  version << "##source=UoO:SNPTools:hapfuseV" << PACKAGE_VERSION;
  if (bcf_hdr_append(m_hdr_out, version.str().c_str()))
    throw runtime_error("could not append hapfuse header");

  bcf_hdr_write(m_fusedVCF, m_hdr_out);
}

void hapfuse::write_site(const Site &osite) const {

  assert(m_hdr_out);
  assert(m_fusedVCF);
  // fill empty record with data and then print
  std::unique_ptr<bcf1_t, void (*)(bcf1_t *)> rec(
      bcf_init1(), [](bcf1_t *b) { bcf_destroy1(b); });
  string alleles;
  for (auto a : osite.all)
    alleles += a + ",";
  alleles.pop_back();
  bcf_update_alleles_str(m_hdr_out, rec.get(), alleles.c_str());
  rec->pos = osite.pos - 1;
  rec->rid = bcf_hdr_name2id(m_hdr_out, osite.chr.c_str());
  assert(rec->rid >= 0);

  size_t numSamps = osite.hap.size() / 2;
  double k = 1.0 / osite.weight;
  bool is_par = (osite.pos >= 60001 && osite.pos <= 2699520) ||
                (osite.pos >= 154931044 && osite.pos <= 155270560);

  vector<unsigned> gts;
  vector<float> lineGPs;
  vector<float> lineAPPs;
  vector<double> GPs(3);

  for (uint i = 0; i < numSamps; i++) {
    const double p0 = osite.hap[i * 2] * k, p1 = osite.hap[i * 2 + 1] * k;

    if (m_init.out_format_tags.at("APP") == true) {
      lineAPPs.push_back(HapfuseHelper::prob2Phred(p0));
      lineAPPs.push_back(HapfuseHelper::prob2Phred(p1));
    }

    const double prr = (1 - p0) * (1 - p1), pra = (1 - p0) * p1 + p0 * (1 - p1),
                 paa = p0 * p1;

    unsigned a, b;
    if (m_init.is_x && (male.find(m_names[i]) != male.end()) && !is_par) {
      if (prr >= paa)
        a = b = 0;
      else
        a = b = 1;
    } else {
      if (prr >= pra && prr >= paa)
        a = b = 0;
      else if (pra >= prr && pra >= paa) {
        if (p0 > p1) {
          a = 1;
          b = 0;
        } else {
          a = 0;
          b = 1;
        }
      } else
        a = b = 1;
    }

    // defined GT
    if (m_init.out_format_tags.at("GT") == true) {
      gts.push_back(bcf_gt_phased(a));
      gts.push_back(bcf_gt_phased(b));
    }

    if (m_init.out_format_tags.at("GP") == true) {
      GPs[0] = (1 - p0) * (1 - p1);
      GPs[1] = (1 - p0) * p1 + p0 * (1 - p1);
      GPs[2] = p0 * p1;

      for (auto &GP : GPs)
        GP = HapfuseHelper::prob2Phred(GP);

      for (auto g : GPs)
        lineGPs.push_back(g);
    }
  } // end numSamps
  if (m_init.out_format_tags.at("GT") == true)
    bcf_update_genotypes(m_hdr_out, rec.get(), gts.data(), gts.size());
  if (m_init.out_format_tags.at("GP") == true)
    bcf_update_format_float(m_hdr_out, rec.get(), "GP", lineGPs.data(),
                            lineGPs.size());
  if (m_init.out_format_tags.at("APP") == true)
    bcf_update_format_float(m_hdr_out, rec.get(), "APP", lineAPPs.data(),
                            lineAPPs.size());

  bcf_write(m_fusedVCF, m_hdr_out, rec.get());
}

bool hapfuse::load_dir(const char *D) {
  string d = D;

  if (d.size() && d[d.size() - 1] != '/')
    d += '/';

  std::unique_ptr<DIR, void (*)(DIR *)> dir(opendir(D),
                                            [](DIR *d) { closedir(d); });

  if (dir == nullptr) {
    cerr << "fail to open " << D << "\n";
    return false;
  }

  struct dirent *ptr;

  while ((ptr = readdir(dir.get())) != NULL) {
    string s = d + ptr->d_name;

    if (s.find(".vcf.gz") != string::npos)
      m_bcfFiles.push_back(s);
  }

  sort(m_bcfFiles.begin(), m_bcfFiles.end());
  return true;
}

void hapfuse::work() {

  vector<double> sum;
  list<Site> outputSites;
  std::future<void> outputFut;
  list<std::future<vector<Site>>> chunkFutures;
  for (uint i = 0; i < m_numInputChunks; i++) {

    // add chunks to the future vector
    // make sure the first chunk is loaded synchronously so that m_name can be
    // filled and the header can be written without hiccups
    if (i == 0)
      chunkFutures.push_back(
          std::async(launch::async, &hapfuse::load_chunk, this, i, true));

    // load up a few chunks in parallel
    // make sure not to go off the end...
    else if (i == 1)
      for (unsigned j = 0; j != 2 && j + i < m_numInputChunks; ++j)
        chunkFutures.push_back(std::async(launch::async, &hapfuse::load_chunk,
                                          this, i + j, false));
    else if (i + 1 < m_numInputChunks)
      chunkFutures.push_back(
          std::async(launch::async, &hapfuse::load_chunk, this, i + 1, false));

    assert(!chunkFutures.empty());
    cout << "Waiting on chunk to load..." << flush;
    vector<Site> chunk = std::move(chunkFutures.front().get());
    cout << " done" << endl;

    chunkFutures.pop_front();

    // write out any sites that have previously been loaded,
    // but are not in the chunk we just loaded
    if (!site.empty()) {
      if (i > 1) {
        cout << "Waiting for merged sites to be written.." << flush;
        outputFut.get();
        cout << " done" << endl;
      }
      outputSites.clear();

      // find first site with the same position as first chunk position
      auto li = site.begin();
      for (; li != site.end(); ++li)
        if (li->pos >= chunk[0].pos)
          break;
      // splice all the site sites that are before chunk[0] in position
      // into outputSites for printing
      outputSites.splice(outputSites.end(), site, site.begin(), li);
      outputFut =
          std::async(launch::async, &hapfuse::write_sites, this, outputSites);
    }

    // find the correct phase
    sum.assign(numSamps() * 2, 0);

    for (list<Site>::iterator li = site.begin(); li != site.end(); ++li)
      for (uint m = 0; m < chunk.size(); m++)
        if (*li == chunk[m]) {
          double *p = &(li->hap[0]), *q = &(chunk[m].hap[0]);

          for (uint j = 0; j < numSamps(); j++) {
            sum[j * 2] += p[j * 2] * q[j * 2] + p[j * 2 + 1] * q[j * 2 + 1];
            sum[j * 2 + 1] += p[j * 2] * q[j * 2 + 1] + p[j * 2 + 1] * q[j * 2];
          }
        }

    // swap phase if needed
    for (uint j = 0; j < numSamps(); j++)
      if (sum[j * 2] < sum[j * 2 + 1]) {
        for (uint m = 0; m < chunk.size(); m++) {
          double t = chunk[m].hap[j * 2];
          chunk[m].hap[j * 2] = chunk[m].hap[j * 2 + 1];
          chunk[m].hap[j * 2 + 1] = t;
        }
      }

    // add chunk to buffer
    merge_chunk(std::move(chunk));
    cout << "Merged chunk " << i + 1 << " of " << m_numInputChunks << endl;
  }

  if (m_numInputChunks > 1) {
    cout << "Waiting for merged sites to be written.." << flush;
    outputFut.get();
    cout << " done" << endl;
  }
  for (auto li : site)
    write_site(li);
}

/*
  finding the beginning, middle, and end of overlap according to list
*/

void hapfuse::find_overlap(vector<Site> chunk, list<Site>::iterator &first,
                           list<Site>::iterator &mid,
                           list<Site>::iterator &last, size_t &chunkMidIdx) {

  assert(!chunk.empty());
  auto li = site.begin();
  size_t m = 0;
  bool firstFound = false;
  vector<size_t> overlapChunkSiteNums;
  overlapChunkSiteNums.reserve(chunk.size());
  while (li != site.end()) {

    if (m == chunk.size())
      break;

    // chunk site is before current site
    // and did not match a previous site
    if (chunk[m].pos < li->pos) {
      if (!m_init.unmatchedSitesOK)
        throw std::runtime_error("Encountered site in overlap region that "
                                 "does not exist in both chunks: " +
                                 chunk[m].chr + ":" + to_string(chunk[m].pos) +
                                 " " + chunk[m].all[0] + " " + chunk[m].all[1]);
      if (!firstFound) {
        first = li;
        --first;
      }
      ++m;
    }

    // chunk site matches list site and needs to be merged in
    else if (chunk[m] == *li) {
      overlapChunkSiteNums.push_back(m);
      ++m;
      ++li;
    }

    // chunk site might match a later site in site
    else
      ++li;
  }

  last = li;

  assert(!overlapChunkSiteNums.empty());
  chunkMidIdx = (overlapChunkSiteNums.size() - 1) >> 1;
  assert(overlapChunkSiteNums.size() > chunkMidIdx);

  // go search for that halfway point again...
  mid = first;
  while (*mid != chunk[chunkMidIdx])
    ++mid;
}

/* How merging works
   . = sites kept
   x = deleted sites
   | = interval between sites at which to cut

   STEP

   Only the number of overlap sites are counted
   The old/left haplotypes are used for the center site
   if the number of overlap sites is odd.
   The example assumes that non-overlapping sites
   in the overlap interval are allowed

   Diplotype 1: . . . .. .|xx x
   Diplotype 2:      xx xx| .  . . . .

*/

void hapfuse::merge_chunk(vector<Site> chunk) {

  if (site.empty()) {
    // move all the sites over to site
    for (size_t m = 0; m != chunk.size(); ++m)
      site.insert(site.end(), std::move(chunk[m]));
  } else {

    if (m_init.ws == HapfuseHelper::WeightingStyle::STEP) {
      auto start = site.begin();
      auto mid = site.begin();
      auto last = site.begin();
      size_t chunkMidIdx = 0;
      find_overlap(chunk, start, mid, last, chunkMidIdx);

      // now delete all sites after halfway point
      auto deleteMe = mid;
      ++deleteMe;
      if (deleteMe != site.end())
        site.erase(deleteMe, site.end());

      // and now insert all sites from chunk after halfway point
      size_t m = chunkMidIdx + 1;
      for (; m < chunk.size(); ++m)
        site.insert(site.end(), std::move(chunk[m]));
    }

    // alternative way of merging
    if (m_init.ws == HapfuseHelper::WeightingStyle::AVERAGE ||
        m_init.ws == HapfuseHelper::WeightingStyle::LINEAR) {
      auto li = site.begin();
      size_t m = 0;
      bool firstOverlapFound = false;
      while (li != site.end()) {

        // exit if we have merged in all chunk sites
        if (m == chunk.size())
          break;

        // chunk site is before current site
        // and did not match a previous site
        if (chunk[m].pos < li->pos) {
          if (!m_init.unmatchedSitesOK)
            throw std::runtime_error("Encountered site in overlap region that "
                                     "does not exist in both chunks: " +
                                     chunk[m].chr + ":" +
                                     to_string(chunk[m].pos) + " " +
                                     chunk[m].all[0] + " " + chunk[m].all[1]);
          site.insert(li, std::move(chunk[m]));
          ++m;
        }
        // chunk site matches list site and needs to be merged in
        else if (chunk[m] == *li) {
          // just measure extent of overlap and fuse together later
          // if we are just cutting overlap in half

          li->weight += chunk[m].weight;

          double *p = &(li->hap[0]), *q = &(chunk[m].hap[0]);
          for (uint j = 0; j < numSamps() * 2; j++)
            p[j] += q[j];

          ++m;
          ++li;
          firstOverlapFound = true;
        }

        // chunk site might match a later site in site
        else {
          if (firstOverlapFound)
            if (!m_init.unmatchedSitesOK)
              throw std::runtime_error(
                  "Encountered site in overlap region that "
                  "does not exist in both chunks: " +
                  chunk[m].chr + ":" + to_string(chunk[m].pos) + " " +
                  chunk[m].all[0] + " " + chunk[m].all[1]);

          ++li;
        }
      }
      // load up all the remaining sites in chunk that are past the last site in
      // site
      for (; m < chunk.size(); ++m)
        site.insert(site.end(), std::move(chunk[m]));
    }
  }
}
