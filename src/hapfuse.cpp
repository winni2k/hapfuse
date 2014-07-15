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
  char all[2];

  static bool is_x;
  static vector<bool> is_male;

  void write(ofile &fusedVCF);
};

bool Site::is_x;
vector<bool> Site::is_male;

void Site::write(ofile &fusedVCF) {
  fusedVCF << chr.c_str() << "\t" << pos << "\t.\t" << all[0] << "\t" << all[1]
           << "\t100\tPASS\t.\tGT:GP:APP";
  uint in = hap.size() / 2;
  double k = 1.0 / cov;
  bool is_par = (pos >= 60001 && pos <= 2699520) ||
                (pos >= 154931044 && pos <= 155270560);

  for (uint i = 0; i < in; i++) {
    uint a, b;
    double p0 = hap[i * 2] * k, p1 = hap[i * 2 + 1] * k;
    double prr = (1 - p0) * (1 - p1), pra = (1 - p0) * p1 + p0 * (1 - p1),
           paa = p0 * p1;

    if (is_x && is_male[i] && !is_par) {
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

    fusedVCF << "\t" << a << "|" << b << ":" << GPs[0] << "," << GPs[1] << ","
             << GPs[2] << ":" << p0 << "," << p1;
  }

  fusedVCF << "\n";
}

class hapfuse {
private:
  list<Site> site;
  uint32_t in;
  FILE *vcf;
  vector<string> file;
  vector<string> name;
  vector<Site> chunk;
  set<string> male;
  vector<string> currentChunkHead;

  void vcf_head(ofile &fusedVCF);
  bool load_chunk(const char *F);
  std::tuple<float, float> extractGP(float *gp, int &gtA, int &gtB);

public:
  static void document(void);
  bool gender(const char *F);
  bool load_dir(const char *D);
  bool load_files(const vector<string> &inFiles);
  void work(string outputFile);
};

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
    GPs[i] = phred2Prob(*gp);
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
    if (gtA == 1)
      pHap1 += GPs[1];
    else
      pHap2 += GPs[1];
  } else {
    pHap1 += GPs[1] / 2;
    pHap2 += GPs[1] / 2;
  }
  return make_tuple(pHap1, pHap2);
}

bool hapfuse::load_chunk(const char *F) {

  // open F and make sure it opened ok
  string inFile(F);
  htsFile *fp = hts_open(inFile.c_str(), "r");
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf1_t *rec = bcf_init1();

  if (!fp)
    throw myException("Could not open file: " + inFile);

  string temp;
  name.clear();
  chunk.clear();

  // skip headers

  // parse #CHROM header line
  name.reserve(bcf_hdr_nsamples(hdr));
  for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
    name.push_back(hdr->samples[i]);
  in = name.size();

  Site::is_male.resize(in);

  for (uint i = 0; i < in; i++)
    Site::is_male[i] = (male.find(name[i]) != male.end());

  // parsing each line of data
  unsigned cnt_lines = 0;
  while (bcf_read1(fp, hdr, rec) >= 0) {

    // get the first five and sample columns
    bcf_unpack(rec, BCF_UN_STR | BCF_UN_IND);

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
    s.chr = bcf_hdr_id2name(hdr, rec->rid);
    s.pos = (rec->pos + 1);

    // make sure site is biallelic
    assert(rec->n_allele == 2);
    string a1(rec->d.allele[0]);
    string a2(rec->d.allele[1]);

    // make sure site is snp
    assert(a1.size() == 1);
    assert(a2.size() == 1);
    s.all[0] = a1[0];
    s.all[1] = a2[0];

    if (!bcf_get_fmt(hdr, rec, "GT"))
      throw std::runtime_error("expected GT field in VCF");
    if (!bcf_get_fmt(hdr, rec, "GP") && !bcf_get_fmt(hdr, rec, "APP"))
      throw std::runtime_error("expected GP or APP field");

    // read sample specific data
    // get genotype array
    int ngt, *gt_arr = NULL, ngt_arr = 0;
    ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);

    if (ngt != 2 * bcf_hdr_nsamples(hdr))
      throw std::runtime_error("number of genotypes found = " + to_string(ngt) +
                               "; with values: " + to_string(gt_arr[0]) + " " +
                               to_string(gt_arr[1]));

    // get APP array
    int n_arr = 0, m_arr = 0;
    float *arr = NULL;
    int stride = 0;
    if (bcf_get_fmt(hdr, rec, "APP")) {
      stride = 2;
      n_arr = bcf_get_format_float(hdr, rec, "APP", &arr, &m_arr);
      assert(n_arr / stride == bcf_hdr_nsamples(hdr));
    }
    // get GP array
    else if (bcf_get_fmt(hdr, rec, "GP")) {
      stride = 3;
      n_arr = bcf_get_format_float(hdr, rec, "GP", &arr, &m_arr);
    } else {
      free(arr);
      throw std::runtime_error("could not read record");
    }

    // cycle through sample vals
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {

      assert(gt_arr[i * 2] != bcf_gt_missing);
      assert(gt_arr[i * 2] != bcf_int32_vector_end);
      assert(gt_arr[i * 2 + 1] != bcf_gt_missing);
      assert(gt_arr[i * 2 + 1] != bcf_int32_vector_end);

      // this assumes the second gt val carries the phase information
      if (!bcf_gt_is_phased(gt_arr[i * 2 + 1]))
        throw std::runtime_error("Error in GT data, genotype is not phased.");

      //    for (uint tokenColIdx = firstSampIdx; tokenColIdx < tokens.size();
      //         ++tokenColIdx) {

      // parse haps

      int gtA = bcf_gt_allele(gt_arr[i]);
      int gtB = bcf_gt_allele(gt_arr[i + 1]);

      // extract/estimate allelic probabilities

      float pHap1;
      float pHap2;
      // parse APPs
      if (bcf_get_fmt(hdr, rec, "APP")) {
        pHap1 = phred2Prob(arr[i * stride]);
        pHap2 = phred2Prob(arr[i * stride + 1]);
      }
      // parse GPs
      else if (bcf_get_fmt(hdr, rec, "GP")) {
        std::tie(pHap1, pHap2) = extractGP(&arr[i * stride], gtA, gtB);
      } else
        throw std::runtime_error("could not load GP or APP field ");

      assert(pHap1 + pHap2 < 2 + EPSILON);

      // make sure pHap1 and 2 are greater than zero
      if (pHap1 < 0)
        pHap1 = 0;

      if (pHap2 < 0)
        pHap2 = 0;

      s.hap[i * 2] = pHap1;
      s.hap[i * 2 + 1] = pHap2;
    }
    chunk.push_back(s);
    free(arr);
    free(gt_arr);
  }

  // gzclose(f);
  return true;
}

void hapfuse::vcf_head(ofile &fusedVCF) {
  string fileFormat("##fileformat=VCFv4.1\n");
  fusedVCF << fileFormat;

  // print the first chunk's header except for its fileformat
  for (auto chunkHeadLine : currentChunkHead)
    if ((!chunkHeadLine.compare(0, 9, "##source=")) ||
        (!chunkHeadLine.compare(0, 21, "##phaseAndImputeCall=")))
      fusedVCF << chunkHeadLine << "\n";

  fusedVCF << "##source=BCM:SNPTools:hapfuseV" << hapfuse_VERSION_MAJOR << "."
           << hapfuse_VERSION_MINOR << "\n";
  fusedVCF
      << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
  fusedVCF << "##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Phred-scaled "
              "genotype posterior probabilities\">\n";

  fusedVCF << "##FORMAT=<ID=APP,Number=2,Type=Float,Description=\"Phred-"
           << "scaled allelic probability, P(Allele=1|Haplotype)\">\n";
  fusedVCF << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
}

bool hapfuse::load_dir(const char *D) {
  string d = D;

  if (d.size() && d[d.size() - 1] != '/')
    d += '/';

  DIR *dir = opendir(D);

  if (dir == NULL) {
    cerr << "fail to open " << D << "\n";
    return false;
  }

  struct dirent *ptr;

  while ((ptr = readdir(dir)) != NULL) {
    string s = d + ptr->d_name;

    if (s.find(".vcf.gz") != string::npos)
      file.push_back(s);
  }

  closedir(dir);
  sort(file.begin(), file.end());
  return true;
}

bool hapfuse::load_files(const vector<string> &inFiles) {

  for (auto inFile : inFiles)
    file.push_back(inFile);

  return true;
}

void hapfuse::work(string outputFile) {

  ofile fusedVCF(outputFile);
  fusedVCF << std::fixed << std::setprecision(3);
  vector<double> sum;

  for (uint i = 0; i < file.size(); i++) {
    if (!load_chunk(file[i].c_str()))
      return;

    // load header here, so we can add first chunk's header to header line
    if (!i) {
      vcf_head(fusedVCF);
      for (uint i = 0; i < in; i++)
        fusedVCF << "\t" << name[i];

      fusedVCF << "\n";
    }

    // write out any sites that have previously been loaded,
    // but are not in the chunk we just loaded
    while (!site.empty() && site.front().pos < chunk[0].pos) {
      site.front().write(fusedVCF);
      site.pop_front();
    }

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

  for (list<Site>::iterator li = site.begin(); li != site.end(); li++)
    li->write(fusedVCF);

  fusedVCF.close();
}

void hapfuse::document(void) {
  cerr << "\nhapfuse v" << hapfuse_VERSION_MAJOR << "."
       << hapfuse_VERSION_MINOR;
  cerr << "\njoint chunked haplotypes into chromosome wide haplotypes";
  cerr << "\nauthor Warren W Kretzschmar @ Marchini Group @ U of Oxford";
  cerr << "\nbased on code by Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
  cerr << "\nusage hapfuse [-g genderFile] [-d inFileDir] < -o out.vcf > < "
          "vcf "
          "files to process in order >";
  cerr << "\n";
  cerr << "\n-g [file]\tFile that indicates which gender each sample is. Only "
          "use for x chromosome.";
  cerr << "\n\t\tExample: NA21522 male";
  cerr << "\n\n";
  exit(1);
}

int main(int argc, char **argv) {
  if (argc < 3)
    hapfuse::document();

  string outFile;
  string genderFile;
  string inFileDir;
  int opt;

  while ((opt = getopt(argc, argv, "d:g:o:")) >= 0) {
    switch (opt) {
    case 'g':
      Site::is_x = true;
      genderFile = optarg;
      break;

    case 'o':
      outFile = optarg;
      break;

    case 'd':
      inFileDir = optarg;
      break;
    }
  }

  // processing vcf files
  vector<string> inFiles;

  for (int index = optind; index < argc; index++)
    inFiles.push_back(argv[index]);

  hapfuse hf;

  if (Site::is_x)
    hf.gender(genderFile.c_str());

  if (!inFileDir.empty())
    hf.load_dir(inFileDir.c_str());
  else if (!inFiles.empty())
    hf.load_files(inFiles);
  else {
    cerr << "no input files" << endl;
    hf.document();
  }

  hf.work(outFile);
  return 0;
}
