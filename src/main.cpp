/* @(#)main.cpp
 */

#include "hapfuse.hpp"

using namespace std;

int main(int argc, char **argv) {
  if (argc < 3)
    hapfuse::document();

  try {
    string genderFile;
    int opt;
    //  size_t numThreads = 1;
    HapfuseHelper::init init;

    static struct option loptions[] = {
      { "gender-file", required_argument, nullptr, 'g' },
      { "output", required_argument, nullptr, 'o' },
      { "output-type", required_argument, nullptr, 'O' },
      { "ligation-method", required_argument, nullptr, 'w' },
      { "wtccc-hap-files", required_argument, nullptr, 'h' },
      { "wtccc-sample-files", required_argument, nullptr, 's' },
      { 0, 0, 0, 0 }
    };

    while ((opt = getopt_long(argc, argv, "d:g:o:O:w:", loptions, nullptr)) >=
           0) {
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
          init.useLinearWeighting = true;
        else
          throw runtime_error("Unexpected option argument -w " +
                              string(optarg));
        break;

      case 'h':
        init.wtcccHapFiles = optarg;
        break;
      case 's':
        init.wtcccSampFiles = optarg;
        break;
      default:
        throw runtime_error("unexpected option: " + std::string(optarg));
      }
    }
    //  omp_set_num_threads(numThreads);

    // processing vcf files
    // also check for file existence
    struct stat buffer;
    for (int index = optind; index < argc; index++) {
      init.inFiles.push_back(argv[index]);
    }

    hapfuse hf(init);

    if (hf.is_x())
      hf.gender(genderFile.c_str());

    hf.work();
  }
  catch (std::exception &e) {
    cerr << "Error: " << e.what() << endl;
    exit(1);
  }
  return 0;
}

