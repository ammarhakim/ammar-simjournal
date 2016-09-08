#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

class NameValuePair {
  public:
/** 
 * @param fname File containing data.
 */
    NameValuePair(const std::string& fname) {
      std::string line;
// open file for reading
      std::ifstream file(fname.c_str());
      if (!file)
        throw std::runtime_error("NameValuePair::NameValuePair: Unable to open input file");
// loop over each line in file
      while (std::getline(file, line)) {
        unsigned commentLoc = line.find('#');
        std::string myLine = line.substr(0, commentLoc);
        if (myLine != "") {
// only work on non-blank line: horribly ugly, but works
          unsigned eqLoc = myLine.find('=');
          std::string lhs, lhsStr, valueStr;
          lhsStr = myLine.substr(0, eqLoc);
          valueStr = myLine.substr(eqLoc+1, std::string::npos);
          std::istringstream lhsToken(lhsStr);
          lhsToken >> lhs;

          double value;
          std::istringstream tokens(valueStr);
          tokens >> value;
// insert it into map
          nvmap[lhs] = value;
        }
      }
    }

/** 
 * @param name Name to lookup
 * @return value associated with name
 */
    double getValue(const std::string& name) const {
      std::map<std::string, double>::const_iterator itr =
        nvmap.find(name);
      if (itr != nvmap.end())
        return itr->second;
      throw std::runtime_error("NameValuePair::getValue: not found");
    }

/**
 * @param name Name to lookup
 * @return true if it exists, false otherwise
 */
    bool hasValue(const std::string& name) const {
      std::map<std::string, double>::const_iterator itr =
        nvmap.find(name);
      if (itr != nvmap.end())
        return true;
      return false;
    }

  private:
    std::map<std::string, double> nvmap;
};

/** Store sim data */
struct QuadData {
    int NDIM; /* Number of dimensions */
    int Np, Nq; /* Number of basis, Number of quadrature points */
    void show() {
      std::cout << "Np: " << Np << " Nq: " << Nq << std::endl;
    }
};

void
vol(int nloop, const QuadData& qd)
{
  int np = qd.Np, nq = qd.Nq;
  int NDIM = qd.NDIM;
  
  Eigen::VectorXd f(np), alpha(np), result(np);
  Eigen::VectorXd fQuad(nq), alphaQuad(nq);
  Eigen::MatrixXd interpMatrix(nq, np);
  Eigen::MatrixXd bigMatrix(np, nq);

  f = Eigen::VectorXd::Random(np);
  alpha = Eigen::VectorXd::Random(np);
  interpMatrix = Eigen::MatrixXd::Random(nq, np);
  bigMatrix = Eigen::MatrixXd::Random(np, nq);

  for (unsigned n=0; n<nloop; ++n)
  {
// interpolate to quadrature nodes    
    fQuad.noalias() = interpMatrix*f; // Np*Nq
    for (unsigned d=0; d<NDIM; ++d)
    {
      alphaQuad.noalias() = interpMatrix*alpha; // Np*Nq
      fQuad.cwiseProduct(alphaQuad); // Nq
// compute updated solution
      result.noalias() = bigMatrix*fQuad; // Np*Nq
    }
  }
}

double cost(const QuadData& qd)
{
  int np = qd.Np, nq = qd.Nq;
  return (double) (3*np*nq+nq);
}

double run(int nloop, const QuadData& qd)
{
  clock_t t1 = clock();
  vol(nloop, qd);
  clock_t t2 = clock();
  return (double) (t2-t1)/CLOCKS_PER_SEC;
}

void printInfo(const QuadData& qd, double tm)
{
  std::cout << qd.NDIM << " " << qd.Np << " " << qd.Nq << " " << tm << std::endl;
}

int
main (int argc, char **argv)
{
  
  if (argc != 2) {
    std::cout << "Usage::" << std::endl;
    std::cout << " dg-int <input-file>" << std::endl;
    std::cout << "  input-file: Name of input file" << std::endl;

    exit(1);
  }
// name of input file and output prefix
  std::string inFile(argv[1]);

// initialize problem state from input file
  NameValuePair nvpair(inFile);

  QuadData forceVolQuad, streamVolQuad;
  forceVolQuad.Np = nvpair.getValue("NpVolForce");
  forceVolQuad.Nq = nvpair.getValue("NqVolForce");
  forceVolQuad.NDIM = nvpair.getValue("VDIM");

  streamVolQuad.Np = nvpair.getValue("NpVolStream");
  streamVolQuad.Nq = nvpair.getValue("NqVolStream");
  streamVolQuad.NDIM = nvpair.getValue("CDIM");
  
  int nloop = nvpair.getValue("nloop");

  std::cout << "# Nloop  " << nloop << std::endl;
  std::cout << "# NDIM | Basis | Quadrature | Time " << std::endl;

// force terms
  double tForce = run(nloop, forceVolQuad);
  printInfo(forceVolQuad, tForce);

// stream terms
  double tStream = run(nloop, streamVolQuad);
  printInfo(streamVolQuad, tStream);

// stats
  std::cout << std::endl;
  std::cout << "# Total time | Time per DOF " << std::endl;  
  double tm = (tForce+tStream)/nloop;
  double tmPerDof = tm/forceVolQuad.Np;
  std::cout << tm << " " << tmPerDof << std::endl;
  
  return 0;
}
