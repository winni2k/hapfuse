/* @(#)main.cpp
 */

#include "hapfuse.hpp"

using namespace std;

int main(int argc, char **argv) {

  ostringstream documentation;
  documentation
      << "\n"
         "Program: hapfuse\n"
         "Version: " << PACKAGE_VERSION
      << "\n\n"
         "Author:  Warren W Kretzschmar @ Marchini Group @ Univ. of Oxford.\n"
         "         Based on code by Yi Wang @ Fuli Yu Group @ BCM-HGSC.\n\n"
         "Usage: hapfuse [options] <-o output_file> <VCF/BCF files to "
         "process>\n\n"
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
         "For details, see README.md\n";
  if (argc < 3) {
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
    vector<string> output_files;

    static struct option loptions[] = {
        {"gender-file", required_argument, nullptr, 'g'},
        {"output", required_argument, nullptr, 'o'},
        {"output-type", required_argument, nullptr, 'O'},
        {"ligation-method", required_argument, nullptr, 'w'},
        {"wtccc-hap-files", required_argument, nullptr, 'h'},
        {"wtccc-sample-files", required_argument, nullptr, 's'},
        {"out_format_tags", required_argument, nullptr, 'T'},
        {"in_format_tags", required_argument, nullptr, 't'},
        {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "d:g:o:O:w:h:s:t:T:", loptions,
                              nullptr)) >= 0) {
      switch (opt) {
      case 'g':
        init.is_x = true;
        genderFile = optarg;
        break;

      case 'o':
        boost::split(output_files, optarg, boost::is_any_of(","));
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

    cout << "Input tags:";
    for (auto it : init.in_format_tags)
      if (it.second)
        cout << " " << it.first;
    cout << endl;
    cout << "Output tags:";
    for (auto it : init.out_format_tags)
      if (it.second)
        cout << " " << it.first;
    cout << endl;

    // figure out what input files are
    if (init.mode == "w") {
      if (output_files.size() == 1) {
        init.outputWTCCCHapsFile = output_files[0] + ".hap.gz";
        init.outputWTCCCSampleFile = output_files[0] + ".sample";
      } else if (output_files.size() == 2) {
        init.outputWTCCCHapsFile = output_files[0];
        init.outputWTCCCSampleFile = output_files[1];
      } else
        throw runtime_error("Could not interpret -o option");
    } else {
      assert(output_files.size() == 1);
      init.outputBCFFile = output_files[0];
    }

    //  omp_set_num_threads(numThreads);

    hapfuse hf(init);

    if (hf.is_x())
      hf.gender(genderFile.c_str());

    hf.work();
  } catch (std::exception &e) {
    cerr << "Error: " << e.what() << endl;
    exit(1);
  }
  return 0;
}
