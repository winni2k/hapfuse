//$Id: utils.h 611 2012-07-18 14:14:00Z koskos $

#ifndef _UTILS_H
#define _UTILS_H

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI 3.14159265358979323846

#define LOW_POS_DOUBLE 1e-300
#define BIG_POS_DOUBLE 1e300
#define LOW_NEG_DOUBLE -1e-300
#define BIG_NEG_DOUBLE -1e300
#define LOW_POS_FLOAT 1e-8
#define BIG_POS_FLOAT 1e8
#define LOW_NEG_FLOAT -1e-8
#define BIG_NEG_FLOAT -1e8
#define BIG_POS_INT 1000000000
#define BIG_NEG_INT -1000000000

#include <string>
#include <vector>
#include <queue>
#include <map>
#include <bitset>
#include <list>
#include <tr1/unordered_map>
#include <bitset>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include <exception>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/lexical_cast.hpp>


using namespace std;
namespace bio = boost::iostreams;
namespace bpo = boost::program_options;
namespace bid = boost::uuids;

/******************************************************/
/*                  UTILS STATISTICS                  */
/******************************************************/
namespace putils {
	void initRandom(long s);
	double getRandom();
	string getRandomID();
	int getRandom(int);
	long getSeed();
	void normalise(vector < double > & v);
	int sample(vector< double > & v, double sum);
	double entropy(vector < double > & v);
	double KLdistance(vector < double > & P, vector < double > & Q);
};

/******************************************************/
/*                  UTILS ALGORITHM                   */
/******************************************************/
namespace autils {
	int max(vector < double > & v);
	int max(vector < int > & v);
    void findUniqueSet(vector < bool > & B, vector < int > & U); // ?
    void decompose(int min, vector < vector < int > > & B, vector < vector < vector < int > > > & BB);//?
	int checkDuo (int pa1, int pa2, int ca1, int ca2);
	int checkTrio (int fa1, int fa2, int ma1, int ma2, int ca1, int ca2);
};

/******************************************************/
/*                  UTILS STRING                      */
/******************************************************/
namespace sutils {
	int tokenize(string &, vector < string > &);
	int tokenize(string &, vector < string > &, int);
	string uint2str(unsigned int n);
	string int2str(int n);
	string int2str(vector < int > & v);
	string long2str(long int n);
	string double2str(double n, int prc = 4);
	string double2str(vector < double > &v, int prc = 4);
	string bool2str(vector<bool> & v);
	string date2str(time_t * t, string format);
};

/******************************************************/
/*                  UTILS FILE                        */
/******************************************************/
namespace futils {
	bool isFile(string f);
	bool createFile(string f);
	string extensionFile(string & filename);
	void bool2binary(vector < bool > & V, ostream &fd);
	bool binary2bool(vector < bool > & V, istream & fd);
};


/******************************************************/
/*                  EXCEPTIONS                        */
/******************************************************/
class myException : public exception {
public:
   explicit myException(std::string msg) : msg_(msg) {}

   virtual ~myException() throw() {}

   virtual const char* what() const throw() {
      return msg_.c_str();
   }

private:
   std::string msg_;
};

/******************************************************/
/*                  INPUT FILE                        */
/******************************************************/
class ifile : public bio::filtering_istream {
private:
	string file;
	ifstream fd;
        bool m_isGood = false;

public:
	ifile();
	ifile(string filename , bool binary = false);
	~ifile();
	string name();
	bool open(string filename, bool binary = false);
	bool readString(string &);
	void close();
        bool isGood(){return m_isGood;};
};

/******************************************************/
/*                  OUTPUT FILE                       */
/******************************************************/
class ofile : public bio::filtering_ostream {
private:
	string file;
	ofstream fd;

public:
	ofile();
	ofile(string filename , bool binary = false);
	~ofile();
	string name();
	bool open(string filename, bool binary = false);
	void writeString(string &);
	void close();
};

/******************************************************/
/*                  LOG FILE                          */
/******************************************************/
class lfile {
private:
	string file;
	ofstream fd;
	bool verboseC;
	bool verboseL;

public:
	lfile();
	~lfile();
	string name();
	bool open(string filename = "file.log");
	void close();
	string getPrefix();
	void muteL();
	void unmuteL();
	void muteC();
	void unmuteC();
	void print(string s);
	void printC(string s);
	void printL(string s);
	void println(string s);
	void printlnC(string s);
	void printlnL(string s);
	void warning(string s);
	void error(string s);
};

#endif
