

#include "writer.hpp"

using namespace hf;
using namespace std;

void Writer::init(WriterHelper::init init) {

  m_init = std::move(init);

  if (m_init.sampNames.empty())
    throw runtime_error("No output sample names");
  if (m_init.chrom.empty())
    throw runtime_error("Please specify output chromosome");

  // figure out what output files are
  // and open output file for writing

  // this is a WTCCC output file
  if (m_init.mode == "w") {
    m_outputFileType = HapfuseHelper::fileType::WTCCC;
    string haps;
    string sample;
    if (m_init.outputFiles.size() == 1) {
      haps = m_init.outputFiles[0] + ".hap.gz";
      sample = m_init.outputFiles[0] + ".sample";
    } else if (m_init.outputFiles.size() == 2) {
      haps = m_init.outputFiles[0];
      sample = m_init.outputFiles[1];
    } else
      throw runtime_error("[Writer] Could not interpret output files [" +
                          accumulate(m_init.outputFiles.begin(),
                                     m_init.outputFiles.end(), string(",")) +
                          "]");

    clog << "[Writer] Writing output to:\n\tWTCCC haps [" << haps << "]"
         << endl;
    m_fusedWTCCCHaps.open(haps);

    clog << "\tWTCCC sample [" << sample << "]" << endl;
    m_fusedWTCCCSample.open(sample);
  }
  // BCF output file
  else if (m_init.mode.find_first_of("buzv") != string::npos) {
    m_outputFileType = HapfuseHelper::fileType::BCF;

    if (m_init.outputFiles.size() != 1)
      throw runtime_error("[Writer] Writing to BCF/VCF; Number of output files "
                          "supplied is not 1 [" +
                          to_string(m_init.outputFiles.size()) + "]");

    assert(!m_fusedVCF);

    clog << "Writing output to:\n\tVCF/BCF [" << m_init.outputFiles[0] << "]"
         << endl;

    string cmode = "w" + m_init.mode;
    m_fusedVCF = hts_open(m_init.outputFiles[0].c_str(), cmode.c_str());

  } else
    throw runtime_error("Encountered unexpected output type [" + m_init.mode +
                        "]");

  write_head();

  if (m_init.is_x)
    loadGenderFile();
}

void Writer::loadGenderFile() {

  assert(!m_init.genderFile.empty());
  ifile gfd(m_init.genderFile);

  if (!gfd.good())
    throw runtime_error("Could not open file [" + gfd.name() + "]");

  string s, g;

  for (gfd >> s >> g; !gfd.eof(); gfd >> s >> g)
    if (g == "male")
      m_maleIDs.insert(s);

  gfd.close();
}

void Writer::write_head() {

  if (m_outputFileType == HapfuseHelper::fileType::BCF)
    write_vcf_head();
  else if (m_outputFileType == HapfuseHelper::fileType::WTCCC)
    write_wtccc_sample();
  else
    throw runtime_error("Encountered unexpected output file type");
}

void Writer::write_wtccc_sample() {

  if (!m_fusedWTCCCSample.good())
    throw runtime_error("Output file is not good [" +
                        m_fusedWTCCCSample.name() + "]");

  m_fusedWTCCCSample << "ID_1 ID_2 missing\n0 0 0\n";
  for (auto n : m_init.sampNames)
    m_fusedWTCCCSample << n << " " << n << " 0\n";
  m_fusedWTCCCSample.close();
}

void Writer::write_vcf_head() {

  assert(!m_hdr_out);
  m_hdr_out = bcf_hdr_init("w");

  // populate header with sample names
  for (auto sampName : m_init.sampNames)
    bcf_hdr_add_sample(m_hdr_out, sampName.c_str());

  bcf_hdr_add_sample(m_hdr_out, NULL);

  // add contig
  bcf_hdr_printf(m_hdr_out, "##contig=<ID=%s>", m_init.chrom.c_str());

  assert(m_fusedVCF);

  if (m_init.GT == true)
    if (bcf_hdr_append(
            m_hdr_out,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"))
      throw runtime_error("could not append VCF GT format header");

  if (m_init.APP == true)
    if (bcf_hdr_append(m_hdr_out,
                       "##FORMAT=<ID=APP,Number=2,Type=Float,Description="
                       "\"Phred-scaled allelic probability, "
                       "P(Allele=1|Haplotype)\">"))
      throw runtime_error("Could not append VCF APP format header");

  if (m_init.GP == true)
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

void Writer::write_site(const Site &osite) {

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

    if (m_init.APP == true) {
      lineAPPs.push_back(HapfuseHelper::prob2Phred(p0));
      lineAPPs.push_back(HapfuseHelper::prob2Phred(p1));
    }

    const double prr = (1 - p0) * (1 - p1), pra = (1 - p0) * p1 + p0 * (1 - p1),
                 paa = p0 * p1;

    unsigned a, b;
    if (m_init.is_x &&
        (m_maleIDs.find(m_init.sampNames[i]) != m_maleIDs.end()) && !is_par) {
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
    if (m_init.GT == true) {
      gts.push_back(bcf_gt_phased(a));
      gts.push_back(bcf_gt_phased(b));
    }

    if (m_init.GP == true) {
      GPs[0] = (1 - p0) * (1 - p1);
      GPs[1] = (1 - p0) * p1 + p0 * (1 - p1);
      GPs[2] = p0 * p1;

      for (auto &GP : GPs)
        GP = HapfuseHelper::prob2Phred(GP);

      for (auto g : GPs)
        lineGPs.push_back(g);
    }
  } // end numSamps

  switch (m_outputFileType) {
  case HapfuseHelper::fileType::BCF:
    write_bcf_site(osite, gts, lineGPs, lineAPPs);
    break;
  case HapfuseHelper::fileType::WTCCC:
    write_wtccc_site(osite, gts);
    break;
  default:
    throw runtime_error("No output file available!");
  }
}

void Writer::write_wtccc_site(const Site &osite, const vector<unsigned> &gts) {

  if (!m_fusedWTCCCHaps.good())
    throw runtime_error("Output file is not good");

  // write out the first five columns
  m_fusedWTCCCHaps << osite.chr << " " << osite.id << " " << osite.pos << " "
                   << osite.all[0] << " " << osite.all[1];

  // print out every allele
  // 3 and 5 are 0 or 1, respectively, increased by one and right-shifted by
  // one
  // see htslib/vcf.h bcf_gt_phased
  for (auto d : gts) {
    if (d == 3)
      m_fusedWTCCCHaps << " 0";
    else if (d == 5)
      m_fusedWTCCCHaps << " 1";
    else
      throw runtime_error("encountered unexpected allele: " + to_string(d));
  }
  m_fusedWTCCCHaps << "\n";
}

void Writer::write_bcf_site(const Site &osite, const vector<unsigned> &gts,
                            const vector<float> &lineGPs,
                            const vector<float> &lineAPPs) const {
  assert(m_hdr_out);

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

  if (m_init.GT == true)
    bcf_update_genotypes(m_hdr_out, rec.get(), gts.data(), gts.size());
  if (m_init.GP == true)
    bcf_update_format_float(m_hdr_out, rec.get(), "GP", lineGPs.data(),
                            lineGPs.size());
  if (m_init.APP == true)
    bcf_update_format_float(m_hdr_out, rec.get(), "APP", lineAPPs.data(),
                            lineAPPs.size());

  bcf_write(m_fusedVCF, m_hdr_out, rec.get());
}
