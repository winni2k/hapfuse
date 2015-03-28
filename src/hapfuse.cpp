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

2015-04-24: Later changes were made by Warren Kretzschmar <wkretzsch@gmail.com>

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

hapfuse::hapfuse(HapfuseHelper::init init)
    : m_init(std::move(init)),
      m_out_GT(m_init.out_format_tags.at("GT") == true),
      m_out_GP(m_init.out_format_tags.at("GP") == true),
      m_out_APP(m_init.out_format_tags.at("APP") == true),
      m_in_GT(m_init.in_format_tags.at("GT") == true),
      m_in_GP(m_init.in_format_tags.at("GP") == true),
      m_in_APP(m_init.in_format_tags.at("APP") == true),
      m_bcfFiles(m_init.cmdLineInputFiles) {

  // make sure requested input and output format tags makes sense
  if (m_out_APP && !m_in_APP)
    throw runtime_error(
        "Need to specify APP as input tag if specifying APP as output tag");
  if (m_out_GP && (!m_in_APP && !m_in_GP))
    throw runtime_error("Need to specify APP or GP as input tag if specifying "
                        "GP as output tag");
  if (m_out_GT && (!m_in_APP && !m_in_GP && !m_in_GT))
    throw runtime_error("Need to specify GT, APP or GP as input tag if "
                        "specifying GT as output tag");

  // Input files will be WTCCC style
  if (!m_init.wtcccHapFilesFile.empty()) {
    m_inputFileType = HapfuseHelper::fileType::WTCCC;

    HapfuseHelper::load_files_from_file(m_init.wtcccHapFilesFile,
                                        m_wtcccHapFiles);
    if (m_init.wtcccSampFilesFile.empty())
      throw std::runtime_error(
          "Need to define WTCCC sample files list with -s if using -h");
    else
      HapfuseHelper::load_files_from_file(m_init.wtcccSampFilesFile,
                                          m_wtcccSampFiles);
    if (m_wtcccHapFiles.size() != m_wtcccSampFiles.size())
      throw std::runtime_error(
          "-h and -s need to contain same number of files");

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

  // start filling writer init struct
  m_writerInit.outputFiles = m_init.outputFiles;
  m_writerInit.mode = m_init.mode;
  m_writerInit.GT = m_out_GT;
  m_writerInit.GP = m_out_GP;
  m_writerInit.APP = m_out_APP;
  m_writerInit.is_x = m_init.is_x;
  m_writerInit.genderFile = m_init.genderFile;
}

// extract the GP fields, convert them to allelic probabilities and store in
// pHap1,2
std::tuple<float, float> hapfuse::extractGP(float *gp, int gtA, int gtB) {

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

  if (!(m_in_GT && m_out_GT) || m_in_APP || m_out_APP || m_in_GP || m_out_GP)
    throw runtime_error("Please specify only GT as input and output "
                        "tags when using WTCCC hap files");

  // load and check samples
  HapSamp chunk(std::move(hapFile), sampFile);

  // Open output VCF for writing
  if (first) {
    assert(m_names.empty());
    m_names = chunk.GetSamps();
    m_writerInit.sampNames = m_names;
    m_writerInit.chrom = chunk.GetChrom();
    m_writer.init(std::move(m_writerInit));
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

    // save names to internal string vector
    assert(m_names.empty());
    m_names.reserve(bcf_hdr_nsamples(hdr.get()));
    for (int i = 0; i < bcf_hdr_nsamples(hdr.get()); ++i)
      m_names.push_back(hdr.get()->samples[i]);

    m_writerInit.chrom = bcf_hdr_id2name(hdr.get(), rec->rid);
    m_writerInit.sampNames = m_names;
    m_writer.init(std::move(m_writerInit));
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
  //  unsigned in = m_names.size();

  // parsing each line of data
  unsigned cnt_lines = 0;
  while (bcf_read1(fp.get(), hdr.get(), rec.get()) >= 0) {

    // get the first five and sample columns
    bcf_unpack(rec.get(), BCF_UN_STR | BCF_UN_IND);

    // site will store the site's information
    Site site;
    site.hap.resize(m_names.size() * 2);
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

    if (!m_in_GT)
      throw runtime_error(
          "Need to specify input GT tag if fusing BCF/VCF files");

    if (m_in_GT && !bcf_get_fmt(hdr.get(), rec.get(), "GT"))
      throw std::runtime_error("expected GT field in VCF");
    if (m_in_GP && !bcf_get_fmt(hdr.get(), rec.get(), "GP"))
      throw std::runtime_error("expected GP field in VCF");
    if (m_in_APP && !bcf_get_fmt(hdr.get(), rec.get(), "APP"))
      throw std::runtime_error("expected APP field in VCF");

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
    if (m_in_GP || m_in_APP) {
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
    } else if (!m_in_GT)
      throw runtime_error("No output tag defined");

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
      if (m_in_GP || m_in_APP) {
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
      } else if (m_in_GT) {
        pHap1 = static_cast<float>(gtA);
        pHap2 = static_cast<float>(gtB);
      } else
        throw runtime_error("no output tag specified");

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
    clog << "Waiting on chunk to load..." << flush;
    vector<Site> chunk = std::move(chunkFutures.front().get());
    clog << " done" << endl;

    chunkFutures.pop_front();

    // write out any sites that have previously been loaded,
    // but are not in the chunk we just loaded
    if (!site.empty()) {
      if (i > 1) {
        clog << "Waiting for merged sites to be written.." << flush;
        outputFut.get();
        clog << " done" << endl;
      }
      clog << "Seeking to beginning of overlap region..." << flush;

      // find first site with the same position as first chunk position
      auto li = site.begin();
      for (; li != site.end(); ++li)
        if (li->pos >= chunk[0].pos)
          break;
      // splice all the site sites that are before chunk[0] in position
      // into outputSites for printing
      list<Site> outputSites;
      outputSites.splice(outputSites.end(), site, site.begin(), li);
      outputFut = std::async(launch::async, &hf::Writer::write_sites,
                             &(this->m_writer), std::move(outputSites));

      clog << " done" << endl;
    }

    // find the correct phase
    clog << "Matching up haplotypes between chunks..." << flush;
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
    clog << " done" << endl;

    // add chunk to buffer
    if (chunk.empty())
      clog << "Skipping chunk " << i + 1 << " of " << m_numInputChunks
           << " because it is empty" << endl;
    else
      clog << "Merging chunk " << i + 1 << " of " << m_numInputChunks
           << "\n\tChunk region: " << chunk.front().chr << ":"
           << chunk.front().pos << "-" << chunk.back().pos << " ... " << flush;
    merge_chunk(std::move(chunk));
    clog << "done" << endl;
  }

  if (m_numInputChunks > 1) {
    clog << "Waiting for merged sites to be written.." << flush;
    outputFut.get();
    clog << " done" << endl;
  }
  for (auto li : site)
    m_writer.write_site(li);
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

// let's return the range of the chunk merged
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
