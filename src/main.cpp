/* @(#)main.cpp
 */

#include "hapfuse.hpp"

using namespace std;

int main(int argc, char **argv) {

  ostringstream documentation;
  documentation
      << "\nhapfuse v" << PACKAGE_VERSION
      << "\njoint chunked haplotypes into chromosome wide haplotypes"
      << "\nauthor Warren W Kretzschmar @ Marchini Group @ U of Oxford"
      << "\nbased on code by Yi Wang @ Fuli Yu' Group @ BCM-HGSC"
      << "\n\nUsage:\thapfuse [options] <-o out.vcf> <VCF/BCF files to "
         "process in order>"
      << "\n\n\t-o, --output <file>\tName of output file"
      << "\n\t-O, --output-type <b|u|z|v>\tOutput file type. b: compressed "
         "BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]"
      << "\n\t-g, --gender-file <file>\tFile that indicates which gender each "
         "sample is. Only use for x chromosome."
      << "\n\t\tExample: NA21522 male"
      << "\n\t-w --ligation-method <string>\tLigation method"
      << "\n\t-h --wtccc-hap-files <file>\tFile containing WTCCC style hap "
         "file names, one per line"
      << "\n\t-h --wtccc-sample-files <file>\tFile containing WTCCC style "
         "sample file names, one per line"
      << "\n\t-t --out_format_tags <string>\tComma separated string of output "
         "format tags (no spaces). Possible tags: GT, GP, APP"
      << "\n\n";

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

    static struct option loptions[] = {
        {"gender-file", required_argument, nullptr, 'g'},
        {"output", required_argument, nullptr, 'o'},
        {"output-type", required_argument, nullptr, 'O'},
        {"ligation-method", required_argument, nullptr, 'w'},
        {"wtccc-hap-files", required_argument, nullptr, 'h'},
        {"wtccc-sample-files", required_argument, nullptr, 's'},
        {"out_format_tags", required_argument, nullptr, 't'},
        {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "d:g:o:O:w:h:s:t:", loptions,
                              nullptr)) >= 0) {
      switch (opt) {
      case 'g':
        init.is_x = true;
        genderFile = optarg;
        break;

      case 'o':
        init.outputFile = optarg;
        break;

      case 'O':
        init.mode = optarg;
        break;

      case 'w':
        if (string(optarg) == "linear")
          init.ws = HapfuseHelper::WeightingStyle::LINEAR;
        else if (string(optarg) == "step")
          init.ws = HapfuseHelper::WeightingStyle::STEP;
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
      case 't':
        boost::split(out_format_tags, optarg, boost::is_any_of(","));
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
    if (out_format_tags.empty()) {
      init.out_format_tags.at("GT") = true;

      // assuming this is input BCF
      if (!init.cmdLineInputFiles.empty()) {
        init.out_format_tags.at("GP") = true;
        init.out_format_tags.at("APP") = true;
      }
    }
    for (auto tag : out_format_tags) {
      try {
        auto val = init.out_format_tags.at(tag);
        val = true;
      } catch (std::out_of_range &e) {
        cerr << "Encountered unexpected output tag [" << tag << endl;
      }
    }

    cout << "Output tags:";
    for (auto it : init.out_format_tags)
      cout << " " << it.first;
    cout << endl;

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
