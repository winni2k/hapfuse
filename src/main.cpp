/* @(#)main.cpp
 */

#include "hapfuse.hpp"

using namespace std;

void printParameters(const HapfuseHelper::init &init);

int main(int argc, char **argv) {

  ostringstream header;
  header << "hapfuse v" << PACKAGE_VERSION << "\n";

  ostringstream documentation;
  documentation
      << "\nAuthor:  Warren W Kretzschmar @ Marchini Group @ Univ. of Oxford.\n"
         "         Based on code by Yi Wang @ Fuli Yu Group @ BCM-HGSC.\n\n"
         "Usage: hapfuse [options] <-o output_file> <VCF/BCF files to "
         "process>\n\n"
         "    -v --verbosity <integer> [0]\n"
         "        Integer between -1 and 1. The larger the integer, the more\n"
         "        verbose the output\n\n"
         "    -o, --output <file> [] - required argument\n"
         "        Name of output file if output is VCF/BCF\n"
         "        Prefix or comma separated haps,sample files if -Ow\n\n"
         "    -O, --output-type <b|u|z|v|w> [b]\n"
         "        Output file type. b: compressed BCF, u: uncompressed BCF,\n"
         "        z: compressed VCF, v: uncompressed VCF,\n"
         "        w: WTCCC style haps/sample\n\n"
         "    -g, --gender-file <file> []\n"
         "        File that indicates which gender each sample is. Only used "
         "        for x chromosome.\n"
         "            Example: NA21522 male\n\n"
         "    -w --ligation-method <string> ['step']\n"
         "        Ligation method. One of 'step', 'linear', or 'average'\n\n"
         "    -h --wtccc-hap-files <file> []\n"
         "        File containing WTCCC style hap file names, one per line\n\n"
         "    -s --wtccc-sample-files <file> []\n"
         "        File containing WTCCC style sample file names, one per\n"
         "        line\n\n"
         "    -t --out_format_tags <string> [GT,GP,APP]\n"
         "        Comma separated string of output format tags (no spaces).\n"
         "        Possible tags: GT, GP, APP\n\n"
         "    -T --in_format_tags <string> [GT,GP,APP]\n"
         "        Comma separated string of input format tags (no spaces).\n"
         "        Possible tags: GT, GP, APP\n\n"
         "    -m --strand_alignment_map <file> []\n"
         "        Space separated file that contains canonical reference\n"
         "        alignment for alleles.\n"
         "        Sites in this file will be used to align input chunk sites.\n"
         "        Example:\n"
         "            20 50000 A T\n"
         "            20 50555 G T\n\n"
         "    -C --assume_chrom <string> []\n"
         "        Ignore the first column of WTCCC and BCF files and use -C\n"
         "        as the chromosome. Do not check if all chunks are from the\n"
         "        same chromosome. \n\n"
         "For details, see README.md\n";
  if (argc < 3) {
    cerr << header.str();
    cerr << documentation.str();
    exit(1);
  }

  try {
    string genderFile;
    int opt;
    //  size_t numThreads = 1;
    HapfuseHelper::init init;
    vector<string> out_format_tags;
    vector<string> in_format_tags;

    static struct option loptions[] = {
        {"gender-file", required_argument, nullptr, 'g'},
        {"output", required_argument, nullptr, 'o'},
        {"output-type", required_argument, nullptr, 'O'},
        {"ligation-method", required_argument, nullptr, 'w'},
        {"wtccc-hap-files", required_argument, nullptr, 'h'},
        {"wtccc-sample-files", required_argument, nullptr, 's'},
        {"out_format_tags", required_argument, nullptr, 'T'},
        {"in_format_tags", required_argument, nullptr, 't'},
        {"strand_alignment_map", required_argument, nullptr, 'm'},
        {"assume_chrom", required_argument, nullptr, 'C'},
        {"verbosity", required_argument, nullptr, 'v'},
        {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "d:g:o:O:w:h:s:t:T:m:C:v:", loptions,
                              nullptr)) >= 0) {
      switch (opt) {
      case 'g':
        init.is_x = true;
        init.genderFile = optarg;
        throw runtime_error("Gender file is not implemented yet.");
        break;

      case 'o':
        boost::split(init.outputFiles, optarg, boost::is_any_of(","));
        break;

      case 'O':
        init.mode = optarg;
        break;

      case 'w':
        if (string(optarg) == "linear")
          init.ws = HapfuseHelper::WeightingStyle::LINEAR;
        else if (string(optarg) == "step")
          init.ws = HapfuseHelper::WeightingStyle::STEP;
        else if (string(optarg) == "average")
          init.ws = HapfuseHelper::WeightingStyle::AVERAGE;
        else
          throw runtime_error("Unexpected option argument -w " +
                              string(optarg));
        break;
      case 'h':
        init.wtcccHapFilesFile = optarg;
        break;
      case 's':
        init.wtcccSampFilesFile = optarg;
        break;
      case 'T':
        boost::split(out_format_tags, optarg, boost::is_any_of(","));
        break;
      case 't':
        boost::split(in_format_tags, optarg, boost::is_any_of(","));
        break;
      case 'm':
        init.alignMapFile = optarg;
        break;
      case 'C':
        init.assumeChrom = optarg;
        break;
      case 'v':
        init.verbosity = stoi(optarg);
        break;
      default:
        throw runtime_error("unexpected option: " + std::string(optarg));
      }
    }

    // processing vcf files
    // also check for file existence
    struct stat buffer;
    for (int index = optind; index < argc; index++)
      init.cmdLineInputFiles.push_back(argv[index]);

    // set format output tags
    if (in_format_tags.empty()) {
      in_format_tags.push_back("GT");

      // assuming this is input BCF
      if (!init.cmdLineInputFiles.empty()) {
        in_format_tags.push_back("GP");
        in_format_tags.push_back("APP");
      }
    }
    if (out_format_tags.empty())
      for (auto t : in_format_tags)
        out_format_tags.push_back(t);

    // now save format tags
    for (auto tag : in_format_tags) {
      try {
        auto &val = init.in_format_tags.at(tag);
        val = true;
      } catch (std::out_of_range &e) {
        cerr << "Encountered unexpected input tag [" << tag << "]" << endl;
      }
    }
    for (auto tag : out_format_tags) {
      try {
        auto &val = init.out_format_tags.at(tag);
        val = true;
      } catch (std::out_of_range &e) {
        cerr << "Encountered unexpected output tag [" << tag << "]" << endl;
      }
    }

    if (init.verbosity > -1) {
      clog << header.str();
      printParameters(init);
    }

    hapfuse hf(init);

    hf.work();
  } catch (std::exception &e) {
    cerr << "Error: " << e.what() << endl;
    exit(1);
  }
  return 0;
}

void printParameters(const HapfuseHelper::init &init) {

  // print input parameters to stderr
  clog << "Parameters:" << endl;
  clog << "* Fusing method: ";
  switch (init.ws) {
  case HapfuseHelper::WeightingStyle::STEP:
    clog << "step";
    break;
  case HapfuseHelper::WeightingStyle::AVERAGE:
    clog << "average";
    break;
  case HapfuseHelper::WeightingStyle::LINEAR:
    clog << "linear";
    break;
  default:
    throw logic_error("Encountered unexpected weighting style");
  }
  clog << endl;

  clog << "* Input tags:";
  for (auto it : init.in_format_tags)
    if (it.second)
      clog << " " << it.first;
  clog << endl;
  clog << "* Output tags:";
  for (auto it : init.out_format_tags)
    if (it.second)
      clog << " " << it.first;
  clog << endl;

  if (!init.alignMapFile.empty())
    clog << "* Alignment map [" + init.alignMapFile + "]" << endl;

  if (!init.assumeChrom.empty())
    clog << "* Assuming chromosome [" + init.assumeChrom + "]\n"
         << "* Ignoring input chunk chromosomes" << endl;
}
