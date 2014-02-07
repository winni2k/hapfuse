/**
@file hapfuse.cpp
@brief joint chunked haplotypes into chromosome wide haplotypes
@author Yi Wang
@date 04/01/2011
*/
#include <iostream>
#include <algorithm>
#include <dirent.h>
#include <stdint.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <string>
#include <vector>
#include <zlib.h>
#include <list>
#include <set>
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <cfloat>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#define BUFFER_SIZE 1000000 // 65536 chars is not big enough for more than ~50,000 samples
#define EPSILON 0.001 // precision of input floats

namespace vcfParse {

//[tutorial_adder_using
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

using qi::double_;
using qi::_1;
using ascii::space;
using phoenix::push_back;


//]

///////////////////////////////////////////////////////////////////////////
//  Our adder
///////////////////////////////////////////////////////////////////////////

//[tutorial_adder
template <typename Iterator>
bool parseProbs(Iterator first, Iterator last, std::vector<double>& probs)
{
    bool r = qi::phrase_parse(first, last,

                              //  Begin grammar
                              (
                                  double_[push_back(phoenix::ref(probs),
                                          _1)] >> *(',' >> double_[push_back(phoenix::ref(probs), _1)])
                              )
                              ,

                              //  End grammar
                              space);

    if (first != last) // fail if we did not get a full match
        return false;

    return r;
}

//]
}

using namespace std;

// converts a probability to phred scale
double prob2Phred(double prob)
{
    assert(prob >= 0 - EPSILON);

    if (prob < 0) prob = 0;

    if (prob == 0) return DBL_MAX_10_EXP * 10;
    else if (prob >= 1) return 2 *
                                   DBL_MIN; // 0 comes out negative zero on my system, not cool!
    else return -10 * log10(prob);
}

//converts a phred scaled probability back to a probability
double phred2Prob(double phred)
{
    assert(phred >= 0);

    if (phred >= DBL_MAX_10_EXP * 10) return DBL_MIN;
    else return pow(10, -phred / 10);
}


struct Site {
    vector<double> hap;
    string chr;
    uint32_t pos;
    uint16_t cov;
    char all[2];

    static bool is_x;
    static vector<bool> is_male;

    void write(FILE *F);
};

bool Site::is_x;
vector<bool> Site::is_male;

void Site::write(FILE *F)
{
    fprintf(F, "%s\t%u\t.\t%c\t%c\t100\tPASS\t.\tGT:GP:APP", chr.c_str(), pos,
            all[0],
            all[1]);
    uint in = hap.size() / 2;
    double k = 1.0 / cov;
    bool is_par = (pos >= 60001 && pos <= 2699520) || (pos >= 154931044
                  && pos <= 155270560);

    for (uint i = 0; i < in; i++) {
        uint a, b;
        double p0 = hap[i * 2] * k, p1 = hap[i * 2 + 1] * k;
        double prr = (1 - p0) * (1 - p1), pra = (1 - p0) * p1 + p0 * (1 - p1),
               paa = p0 * p1;

        if (is_x && is_male[i] && !is_par) {
            if (prr >= paa) a = b = 0;
            else a = b = 1;
        } else {
            if (prr >= pra && prr >= paa) a = b = 0;
            else if (pra >= prr && pra >= paa) {
                if (p0 > p1) {
                    a = 1;
                    b = 0;
                } else {
                    a = 0;
                    b = 1;
                }
            } else a = b = 1;
        }

        vector<double> GPs(3);
        GPs[0] = (1 - p0) * (1 - p1);
        GPs[1] = (1 - p0) * p1 + p0 * (1 - p1);
        GPs[2] = p0 * p1;

        for (auto & GP : GPs) GP = prob2Phred(GP);

        p0 = prob2Phred(p0);
        p1 = prob2Phred(p1);

        fprintf(F, "\t%u|%u:%.3f,%.3f,%.3f:%.3f,%.3f", a, b, GPs[0], GPs[1], GPs[2], p0,
                p1);
    }

    fprintf(F, "\n");
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

    void vcf_head(const char *F);
    bool load_chunk(const char *F);
public:
    static void document(void);
    bool gender(const char *F);
    bool load_dir(const char *D);
    bool load_files(const vector<string> & inFiles);
    void work(const char *F);
};

bool hapfuse::gender(const char *F)
{
    ifstream fi(F);

    if (!fi) {
        cerr << "fail to open " << F << endl;
        return false;
    }

    string s, g;

    for (fi >> s >> g; !fi.eof(); fi >> s >> g) if (g == "male") male.insert(s);

    fi.close();
    return true;
}

bool hapfuse::load_chunk(const char *F)
{
    char buff[BUFFER_SIZE];
    string temp;
    name.clear();
    chunk.clear();
    gzFile f = gzopen(F, "rt");

    if (f == Z_NULL) return false;

    do gzgets(f, buff, BUFFER_SIZE);

    while (strstr(buff, "#CHROM") == NULL);

    {
        istringstream si(buff);

        for (uint i = 0; i < 9; i++) si >> temp;

        temp = "";

        while (!si.eof()) {
            si >> temp;

            if (temp != "") name.push_back(temp);

            temp = "";
        }
    }

    in = name.size();
    Site::is_male.resize(in);

    for (uint i = 0; i < in;
            i++) Site::is_male[i] = (male.find(name[i]) != male.end());

    Site s;
    s.hap.resize(in * 2);
    s.cov = 1;
    string a, b;

    while (!gzeof(f)) {
        gzgets(f, buff, BUFFER_SIZE);
        vector<string> tokens;
        boost::split(tokens, buff, boost::is_any_of("\t"));
        s.chr = tokens[0];
        s.pos = stoi(tokens[1]);
        assert(tokens[3].size() == 1);
        assert(tokens[4].size() == 1);
        s.all[0] = tokens[3][0];
        s.all[1] = tokens[4][0];

        //        cerr << s.chr << ":" << s.pos << endl;

        vector<string> GTFields;
        boost::split(GTFields, tokens[8], boost::is_any_of(":"));
        int GTIdx = -1;
        int GPIdx = -1;
        int APPIdx = -1;

        for (unsigned fieldIdx = 0; fieldIdx < GTFields.size(); ++fieldIdx) {
            if (GTFields[fieldIdx].compare("GT") == 0) GTIdx = fieldIdx;
            else if (GTFields[fieldIdx].compare("GP") == 0) GPIdx = fieldIdx;
            else if (GTFields[fieldIdx].compare("APP") == 0) APPIdx = fieldIdx;
        }

        if (!(GTIdx >= 0 && (GPIdx >= 0 || APPIdx >= 0))) {
            cerr << "expected GT:GP or GT:APP input format in chunk " << F << endl;
            exit(1);
        }


        // read sample specific data
        for (uint tokenColIdx = 9; tokenColIdx < tokens.size(); ++tokenColIdx) {

            //            cerr << " " << tokenColIdx;
            vector<string> sampDat;
            boost::split(sampDat, tokens[tokenColIdx], boost::is_any_of(":"));
            assert(sampDat.size() > 1);

            // parse haps
            string GT = sampDat[GTIdx];

            if (GT.at(1) != '|' || GT.size() != 3) {
                cerr << "Error in GT data, could not split on '|': " << GT << endl;
                exit(1);
            }


            // extract/estimate allelic probabilities
            double pHap1, pHap2;

            if (APPIdx >= 0) {
                vector<double> APPs;

                if (!vcfParse::parseProbs(sampDat[APPIdx].begin(), sampDat[APPIdx].end(),
                                          APPs)) {
                    cerr << "Could not parse: " << sampDat[GPIdx] << endl;
                    exit(1);
                }


                //                boost::split(inDat, , boost::is_any_of(","));
                assert(APPs.size() == 2);

                // convert GPs to probabilities
                double sum = 0;

                for (auto& APP : APPs) {
                    APP = phred2Prob(APP);
                    sum += APP;
                }

                assert(sum < 2 + EPSILON);
                pHap1 = APPs[0];
                pHap2 = APPs[1];
            }

            // parse GPs
            else if (GPIdx >= 0) {

                vector<double> GPs;

                if (!vcfParse::parseProbs(sampDat[GPIdx].begin(), sampDat[GPIdx].end(), GPs)) {
                    cerr << "Could not parse: " << sampDat[GPIdx] << endl;
                    exit(1);
                }

                //                boost::split(inDat, sampDat[GPIdx], boost::is_any_of(","));
                assert(GPs.size() == 3);

                // convert GPs to probabilities
                double sum = 0;

                for (auto& GP : GPs) {
                    GP = phred2Prob(GP);
                    sum += GP;
                }

                assert(fabs(sum - 1) < EPSILON);
                pHap1 = GPs[2];
                pHap2 = GPs[2];

                // swap alleles if evidence exists in GT field
                if (GT.at(0) != GT.at(2)) {
                    if (GT.at(0) == '1')
                        pHap1 += GPs[1];
                    else
                        pHap2 += GPs[1];
                } else {
                    pHap1 += GPs[1] / 2;
                    pHap2 += GPs[1] / 2;
                }
            } else {
                cerr << "could not load GP or APP field: " << tokens[tokenColIdx] << endl;
                exit(1);
            }

            // make sure pHap1 and 2 are greater than zero
            if (pHap1 < 0) pHap1 = 0;

            if (pHap2 < 0) pHap2 = 0;

            // assign allelic probs
            s.hap[(tokenColIdx - 9) * 2] = pHap1;
            s.hap[(tokenColIdx - 9) * 2 + 1] = pHap2;
        }

        chunk.push_back(s);
    }

    gzclose(f);
    return true;
}

void hapfuse::vcf_head(const char *F)
{
    vcf = fopen(F, "wt");
    fprintf(vcf, "##fileformat=VCFv4.1\n");
    fprintf(vcf, "##source=BCM:SNPTools:hapfuse\n");
    fprintf(vcf,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(vcf,
            "##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Phred-scaled genotype posterior probabilities\">\n");

    fprintf(vcf,
            "##FORMAT=<ID=APP,Number=2,Type=Float,Description=\"Phred-scaled allelic probability, P(Allele=1|Haplotype)\">\n");
    fprintf(vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
}

bool hapfuse::load_dir(const char *D)
{
    string d = D;

    if (d.size() && d[d.size() - 1] != '/') d += '/';

    DIR *dir = opendir(D);

    if (dir == NULL) {
        cerr << "fail to open " << D << endl;
        return false;
    }

    struct dirent  *ptr;

    while ((ptr = readdir(dir)) != NULL) {
        string s = d + ptr->d_name;

        if (s.find(".vcf.gz") != string::npos) file.push_back(s);
    }

    closedir(dir);
    sort(file.begin(), file.end());
    return true;
}

bool hapfuse::load_files(const vector<string> & inFiles)
{

    for (auto inFile : inFiles)
        file.push_back(inFile);

    return true;
}

void hapfuse::work(const char *F)
{
    vcf_head(F);
    vector<double> sum;

    for (uint i = 0; i < file.size(); i++) {
        if (!load_chunk(file[i].c_str())) return;

        if (!i) {
            for (uint i = 0; i < in; i++) fprintf(vcf, "\t%s", name[i].c_str());

            fprintf(vcf, "\n");
        }

        while (!site.empty() && site.front().pos < chunk[0].pos) {
            site.front().write(vcf);
            site.pop_front();
        }

        // find the correct phase
        sum.assign(in * 2, 0);

        for (list<Site>::iterator li = site.begin(); li != site.end(); ++ li)
            for (uint m = 0; m < chunk.size(); m++) if (li->pos == chunk[m].pos) {
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

            for (list<Site>::iterator li = site.begin(); li != site.end();
                    li++) if (li->pos == chunk[m].pos) {
                    found = true;
                    li->cov++;
                    double *p = &(li->hap[0]), *q = &(chunk[m].hap[0]);

                    for (uint j = 0; j < in * 2; j++) p[j] += q[j];
                }

            if (!found) site.push_back(chunk[m]);
        }

        cerr << file[i] << endl;
    }

    for (list<Site>::iterator li = site.begin(); li != site.end();
            li++) li->write(vcf);

    fclose(vcf);
}

void hapfuse::document(void)
{
    cerr << "\nhapfuse";
    cerr << "\njoint chunked haplotypes into chromosome wide haplotypes";
    cerr << "\nauthor Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
    cerr << "\nusage hapfuse [-g genderFile] [-d inFileDir] < -o out.vcf > < vcf files to process in order >";
    cerr << "\n";
    cerr << "\n-g [file]\tFile that indicates which gender each sample is. Only use for x chromosome.";
    cerr << "\n\t\tExample: NA21522 male";
    cerr << "\n\n";
    exit(1);
}

int main(int argc, char **argv)
{
    if (argc < 3) hapfuse::document();

    string outFile;
    string genderFile;
    string inFileDir;
    int opt;

    while ((opt = getopt(argc, argv, "d:g:o:")) >= 0) {
        switch (opt) {
        case    'g':
            Site::is_x = true;
            genderFile = optarg;
            break;

        case    'o':
            outFile = optarg;
            break;

        case    'd':
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

    hf.work(outFile.c_str());
    return 0;
}











