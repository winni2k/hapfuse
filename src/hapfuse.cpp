/**
@file hapfuse.cpp
@brief joint chunked haplotypes into chromosome wide haplotypes
@author Yi Wang
@date 04/01/2011
*/
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

#include "version.hpp"
#include "utils.hpp"

#define BUFFER_SIZE                                                            \
  1000000 // 65536 chars is not big enough for more than ~50,000 samples
#define EPSILON 0.001 // precision of input floats

using namespace std;

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

struct Site {
  vector<double> hap;
  string chr;
  uint32_t pos;
  uint16_t cov;
  vector<string> all;

  void write(ofile &fusedVCF);
};

class hapfuse {
private:
  list<Site> site;
  vector<string> file;
  set<string> male;
  bool m_is_x;
  vector<string> m_names;

  void write_vcf_head();
  vector<Site> load_chunk(const char *F, bool first);
  std::tuple<float, float> extractGP(float *gp, int &gtA, int &gtB);

  htsFile *m_fusedVCF = NULL;
  bcf_hdr_t *m_hdr_out = NULL;
  void write_site(const Site &osite) const;
  void write_sites(const list<Site> &outSites) {
    for (auto s : outSites)
      write_site(s);
  }

public:
  static void document(void);
  bool is_x() { return m_is_x; }
  bool gender(const char *F);
  bool load_dir(const char *D);
  bool load_files(const vector<string> &inFiles);
  void work();

  // destroy output header and output file descriptor
  hapfuse(const string &outputFile, const string &mode, bool is_x);
  ~hapfuse() {
    bcf_hdr_destroy(m_hdr_out);
    hts_close(m_fusedVCF);
  }
};

hapfuse::hapfuse(const string &outputFile, const string &mode, bool is_x)
    : m_is_x(is_x) {

  // open output file for writing
  assert(!m_fusedVCF);
  if (outputFile.size() == 0)
    throw runtime_error("Please specify an output file");

  size_t found = mode.find_first_of("buzv");
  if (mode.size() > 0 && (mode.size() != 1 || found == string::npos))
    throw runtime_error("-O needs to be 'b' for compressed bcf, 'u' for "
                        "uncompressed bcf, 'z' for compressed vcf or 'v' for "
                        "uncompressed vcf");
  string cmode = "w" + mode;
  m_fusedVCF = hts_open(outputFile.c_str(), cmode.c_str());
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
    GPs.push_back(phred2Prob(*gp));
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

vector<Site> hapfuse::load_chunk(const char *F, bool first) {

  // open F and make sure it opened ok
  string inFile(F);
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
    write_vcf_head();
  }

  string temp;
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

    // s will store the site's information
    Site s;
    s.hap.resize(in * 2);
    s.cov = 1;
    string a, b;

    ++cnt_lines;

    // tokenize on "\t"
    // fill in site information (s)
    //    vector<string> tokens;
    //    boost::split(tokens, buffer, boost::is_any_of("\t"));
    s.chr = bcf_hdr_id2name(hdr.get(), rec->rid);
    s.pos = (rec->pos + 1);

    // make sure site is biallelic
    assert(rec->n_allele == 2);
    string a1(rec->d.allele[0]);
    string a2(rec->d.allele[1]);

    // make sure site is snp
    s.all.push_back(a1);
    s.all.push_back(a2);

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
        pHap1 = phred2Prob(arr[sampNum * stride]);
        pHap2 = phred2Prob(arr[sampNum * stride + 1]);
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

      s.hap[sampNum * 2] = pHap1;
      s.hap[sampNum * 2 + 1] = pHap2;
    }
    chunk.push_back(s);
    free(arr);
    free(gt_arr);
  }

  return chunk;
}

void hapfuse::write_vcf_head() {

  assert(m_hdr_out);
  assert(m_fusedVCF);

  if (bcf_hdr_append(m_hdr_out,
                     "##FORMAT=<ID=APP,Number=2,Type=Float,Description="
                     "\"Phred-scaled allelic probability, "
                     "P(Allele=1|Haplotype)\">"))
    throw runtime_error("Could not append VCF APP format header");

  if (bcf_hdr_append(
          m_hdr_out,
          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"))
    throw runtime_error("could not append VCF GP format header");

  if (bcf_hdr_append(
          m_hdr_out,
          "##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Phred-scaled "
          "genotype posterior probabilities\">"))
    throw runtime_error("could not append VCF GP format header");

  ostringstream version;
  version << "##source=UoO:SNPTools:hapfuseV" << hapfuse_VERSION_MAJOR << "."
          << hapfuse_VERSION_MINOR;
  if (bcf_hdr_append(m_hdr_out, version.str().c_str()))
    throw runtime_error("could not append hapfuse header");

  bcf_hdr_write(m_fusedVCF, m_hdr_out);

  // save names to internal string vector
  assert(m_names.size() == 0);
  m_names.reserve(bcf_hdr_nsamples(m_hdr_out));
  for (int i = 0; i < bcf_hdr_nsamples(m_hdr_out); ++i)
    m_names.push_back(m_hdr_out->samples[i]);
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

  uint in = osite.hap.size() / 2;
  double k = 1.0 / osite.cov;
  bool is_par = (osite.pos >= 60001 && osite.pos <= 2699520) ||
                (osite.pos >= 154931044 && osite.pos <= 155270560);

  vector<unsigned> gts;
  vector<float> lineGPs;
  vector<float> lineAPPs;

  for (uint i = 0; i < in; i++) {
    uint a, b;
    double p0 = osite.hap[i * 2] * k, p1 = osite.hap[i * 2 + 1] * k;
    double prr = (1 - p0) * (1 - p1), pra = (1 - p0) * p1 + p0 * (1 - p1),
           paa = p0 * p1;

    if (m_is_x && (male.find(m_names[i]) != male.end()) && !is_par) {
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

    vector<double> GPs(3);
    GPs[0] = (1 - p0) * (1 - p1);
    GPs[1] = (1 - p0) * p1 + p0 * (1 - p1);
    GPs[2] = p0 * p1;

    for (auto &GP : GPs)
      GP = prob2Phred(GP);

    p0 = prob2Phred(p0);
    p1 = prob2Phred(p1);

    gts.push_back(bcf_gt_phased(a));
    gts.push_back(bcf_gt_phased(b));
    for (auto g : GPs)
      lineGPs.push_back(g);
    lineAPPs.push_back(p0);
    lineAPPs.push_back(p1);
    //    fusedVCF << "\t" << a << "|" << b << ":" << GPs[0] << "," << GPs[1] <<
    // ","
    //           << GPs[2] << ":" << p0 << "," << p1;
  }
  bcf_update_genotypes(m_hdr_out, rec.get(), gts.data(), gts.size());
  bcf_update_format_float(m_hdr_out, rec.get(), "GP", lineGPs.data(),
                          lineGPs.size());
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
      file.push_back(s);
  }

  sort(file.begin(), file.end());
  return true;
}

bool hapfuse::load_files(const vector<string> &inFiles) {

  for (auto inFile : inFiles)
    file.push_back(inFile);

  return true;
}

void hapfuse::work() {

  vector<double> sum;
  list<Site> outputSites;
  std::future<void> outputFut;
  list<std::future<vector<Site> > > chunkFutures;
  for (uint i = 0; i < file.size(); i++) {

    // add chunks to the future vector
    // make sure the first chunk is loaded synchronously so that m_name can be
    // filled and the header can be written without hiccups
    if (i == 0)
      chunkFutures.push_back(std::async(launch::async, &hapfuse::load_chunk,
                                        this, file[i].c_str(), i == 0));
    // load up a few chunks in parallel
    // make sure not to go off the end...
    else if (i == 1)
      for (unsigned j = 0; j != 2 && j + i < file.size(); ++j)
        chunkFutures.push_back(std::async(launch::async, &hapfuse::load_chunk,
                                          this, file[i + j].c_str(), false));
    else if (i + 1 < file.size())
      chunkFutures.push_back(std::async(launch::async, &hapfuse::load_chunk,
                                        this, file[i + 1].c_str(), false));

    assert(!chunkFutures.empty());
    vector<Site> chunk = std::move(chunkFutures.front().get());
    chunkFutures.pop_front();
    unsigned in = m_names.size();

    // write out any sites that have previously been loaded,
    // but are not in the chunk we just loaded
    if (i != 0)
      outputFut.get();
    outputSites.clear();
    while (!site.empty() && site.front().pos < chunk[0].pos) {
      outputSites.splice(outputSites.end(), site, site.begin());
    }
    outputFut =
        std::async(launch::async, &hapfuse::write_sites, this, outputSites);

    // find the correct phase
    sum.assign(in * 2, 0);

    for (list<Site>::iterator li = site.begin(); li != site.end(); ++li)
      for (uint m = 0; m < chunk.size(); m++)
        if (li->pos == chunk[m].pos) {
          double *p = &(li->hap[0]), *q = &(chunk[m].hap[0]);

          for (uint j = 0; j < in; j++) {
            sum[j * 2] += p[j * 2] * q[j * 2] + p[j * 2 + 1] * q[j * 2 + 1];
            sum[j * 2 + 1] += p[j * 2] * q[j * 2 + 1] + p[j * 2 + 1] * q[j * 2];
          }
        }

    // swap phase if needed
    for (uint j = 0; j < in; j++)
      if (sum[j * 2] < sum[j * 2 + 1]) {
        for (uint m = 0; m < chunk.size(); m++) {
          double t = chunk[m].hap[j * 2];
          chunk[m].hap[j * 2] = chunk[m].hap[j * 2 + 1];
          chunk[m].hap[j * 2 + 1] = t;
        }
      }

    // add chunk to buffer
    for (uint m = 0; m < chunk.size(); m++) {
      bool found = false;

      for (list<Site>::iterator li = site.begin(); li != site.end(); li++)
        if (li->pos == chunk[m].pos) {
          found = true;
          li->cov++;
          double *p = &(li->hap[0]), *q = &(chunk[m].hap[0]);

          for (uint j = 0; j < in * 2; j++)
            p[j] += q[j];
        }

      if (!found)
        site.push_back(chunk[m]);
    }

    cout << file[i] << endl;
  }

  outputFut.get();
  for (list<Site>::iterator li = site.begin(); li != site.end(); li++)
    write_site(*li);
}

void hapfuse::document(void) {
  cerr << "\nhapfuse v" << hapfuse_VERSION_MAJOR << "."
       << hapfuse_VERSION_MINOR;
  cerr << "\njoint chunked haplotypes into chromosome wide haplotypes";
  cerr << "\nauthor Warren W Kretzschmar @ Marchini Group @ U of Oxford";
  cerr << "\nbased on code by Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
  cerr << "\n\nUsage:\thapfuse [options] <-o out.vcf> <VCF/BCF files to "
          "process in "
          "order>";

  cerr << "\n\n\t-o <file>\tName of output file";

  cerr << "\n\t-O <b|u|z|v>\tOutput file type. b: compressed BCF, u: "
          "uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]";

  cerr
      << "\n\t-g <file>\tFile that indicates which gender each sample is. Only "
         "use for x chromosome.";
  cerr << "\n\t\tExample: NA21522 male";
  cerr << "\n\n";
  exit(1);
}

int main(int argc, char **argv) {
  if (argc < 3)
    hapfuse::document();

  try {
    string outFile;
    string genderFile;
    string inFileDir;
    int opt;
    string mode = "v"; // default is compressed bcf
    bool is_x = false;
    //  size_t numThreads = 1;

    while ((opt = getopt(argc, argv, "d:g:o:O:")) >= 0) {
      switch (opt) {
      case 'g':
        is_x = true;
        genderFile = optarg;
        break;

      case 'o':
        outFile = optarg;
        break;

      case 'O':
        mode = optarg;
        break;

      case 'd':
        inFileDir = optarg;
        break;
      default:
        throw runtime_error("unexpected option: " + std::string(optarg));
      }
    }
    //  omp_set_num_threads(numThreads);

    // processing vcf files
    vector<string> inFiles;

    // also check for file existence
    struct stat buffer;
    for (int index = optind; index < argc; index++) {
      if (stat(argv[index], &buffer) == 0)
        inFiles.push_back(argv[index]);
      else
        throw runtime_error("Input file does not exist [" +
                            string(argv[index]) + "]");
    }

    hapfuse hf(outFile, mode, is_x);

    if (hf.is_x())
      hf.gender(genderFile.c_str());

    if (!inFileDir.empty())
      hf.load_dir(inFileDir.c_str());
    else if (!inFiles.empty())
      hf.load_files(inFiles);
    else {
      cerr << "no input files" << endl;
      hf.document();
    }

    hf.work();
  }
  catch (std::exception &e) {
    cerr << "Error: " << e.what() << endl;
    exit(1);
  }
  return 0;
}
