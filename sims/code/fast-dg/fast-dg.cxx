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

class Update1x3vSerendip {
  public:
    void calcDirectStreamVol(double w, const std::vector<double>& f, std::vector<double>& V)
    {
      using namespace std;
      double sq3 = sqrt(3), sq5 = sqrt(5), sq1_5 = 1/sqrt(5);
  
      V[1] += sq3*f[0]*w+f[2];
      V[5] += sq3*f[2]*w+(2*f[12])*sq1_5+f[0];
      V[6] += sq3*f[3]*w+f[7];
      V[8] += sq3*f[4]*w+f[9];
      V[11] += sq3*sq5*f[1]*w+sq5*f[5];
      V[15] += sq3*f[7]*w+(2*f[22])*sq1_5+f[3];
      V[16] += sq3*f[9]*w+(2*f[26])*sq1_5+f[4];
      V[17] += sq3*f[10]*w+f[18];
      V[19] += sq3*sq5*f[5]*w+2*f[20]+sq5*f[1];
      V[20] += sq3*f[12]*w+(2*f[2])*sq1_5;
      V[21] += sq3*sq5*f[6]*w+sq5*f[15];
      V[23] += sq3*f[13]*w+f[24];
      V[25] += sq3*sq5*f[8]*w+sq5*f[16];
      V[28] += sq3*f[14]*w+f[29];
      V[31] += sq3*f[18]*w+(2*f[38])*sq1_5+f[10];
      V[32] += sq3*sq5*f[15]*w+2*f[33]+sq5*f[6];
      V[33] += sq3*f[22]*w+(2*f[7])*sq1_5;
      V[34] += sq3*f[24]*w+f[13];
      V[35] += sq3*sq5*f[16]*w+2*f[36]+sq5*f[8];
      V[36] += sq3*f[26]*w+(2*f[9])*sq1_5;
      V[37] += sq3*sq5*f[17]*w+sq5*f[31];
      V[39] += sq3*f[27]*w+f[40];
      V[41] += sq3*f[29]*w+f[14];
      V[42] += sq3*f[30]*w+f[43];
      V[44] += sq3*sq5*f[31]*w+2*f[45]+sq5*f[17];
      V[45] += sq3*f[38]*w+(2*f[18])*sq1_5;
      V[46] += sq3*f[40]*w+f[27];
      V[47] += sq3*f[43]*w+f[30];
    }

    void calcDirectForceVol(double w, const std::vector<double>& f, const std::vector<double>& EB, std::vector<double>& V)
    { // this is only place holder for now
      register double EB0 = EB[0], EB1 = EB[1], EB2 = EB[2];
  
      V[2] += 1.224744871391589*EB2*f[11]+1.224744871391589*f[1]*EB1
        +1.224744871391589*f[0]*EB0;
      V[5] += 1.095445115010332*EB1*f[11]+1.095445115010332*f[1]*EB2
        +1.224744871391589*f[0]*EB1
        +1.224744871391589*EB0*f[1];
      V[7] += 1.224744871391589*EB2*f[21]+1.224744871391589*EB1*f[6]
        +1.224744871391589*EB0*f[3];
      V[9] += 1.224744871391589*EB2*f[25]+1.224744871391589*EB1*f[8]
        +1.224744871391589*EB0*f[4];
      V[12] += 2.738612787525831*EB2*f[19]+2.738612787525831*EB1*f[5]
        +2.738612787525831*EB0*f[2];
      V[15] += 1.095445115010332*EB1*f[21]+1.095445115010332*EB2*f[6]
        +1.224744871391589*EB0*f[6]
        +1.224744871391589*EB1*f[3];
      V[16] += 1.095445115010332*EB1*f[25]+1.095445115010332*EB2*f[8]
        +1.224744871391589*EB0*f[8]
        +1.224744871391589*EB1*f[4];
      V[18] += 1.224744871391589*EB2*f[37]+1.224744871391589*EB1*f[17]
        +1.224744871391589*EB0*f[10];
      V[19] += 0.7824607964359517*EB2*f[11]+1.224744871391589*EB0*f[11]
        +1.224744871391589*f[0]*EB2
        +1.095445115010332*f[1]*EB1;
      V[20] += 2.449489742783178*EB1*f[19]+2.449489742783178*EB2*f[5]
        +2.738612787525831*EB0*f[5]
        +2.738612787525831*f[2]*EB1;
      V[22] += 2.738612787525831*EB2*f[32]+2.738612787525831*EB1*f[15]
        +2.738612787525831*EB0*f[7];
      V[24] += 1.224744871391589*EB1*f[23]+1.224744871391589*EB0*f[13];
      V[26] += 2.738612787525831*EB2*f[35]+2.738612787525831*EB1*f[16]
        +2.738612787525831*EB0*f[9];
      V[29] += 1.224744871391589*EB1*f[28]+1.224744871391589*EB0*f[14];
      V[31] += 1.095445115010332*EB1*f[37]+1.095445115010332*EB2*f[17]
        +1.224744871391589*EB0*f[17]
        +1.224744871391589*EB1*f[10];
      V[32] += 0.7824607964359517*EB2*f[21]+1.224744871391589*EB0*f[21]
        +1.095445115010332*EB1*f[6]
        +1.224744871391589*f[3]*EB2;
      V[33] += 2.449489742783178*EB1*f[32]+2.449489742783178*EB2*f[15]
        +2.738612787525831*EB0*f[15]
        +2.738612787525831*EB1*f[7];
      V[34] += 1.095445115010332*EB2*f[23]+1.224744871391589*EB0*f[23]
        +1.224744871391589*EB1*f[13];
      V[35] += 0.7824607964359517*EB2*f[25]+1.224744871391589*EB0*f[25]
        +1.095445115010332*EB1*f[8]
        +1.224744871391589*EB2*f[4];
      V[36] += 2.449489742783178*EB1*f[35]+2.449489742783178*EB2*f[16]
        +2.738612787525831*EB0*f[16]
        +2.738612787525831*EB1*f[9];
      V[38] += 2.738612787525831*EB2*f[44]+2.738612787525831*EB1*f[31]
        +2.738612787525831*EB0*f[18];
      V[40] += 1.224744871391589*EB1*f[39]+1.224744871391589*EB0*f[27];
      V[41] += 1.095445115010332*EB2*f[28]+1.224744871391589*EB0*f[28]
        +1.224744871391589*EB1*f[14];
      V[43] += 1.224744871391589*EB1*f[42]+1.224744871391589*EB0*f[30];
      V[44] += 0.7824607964359517*EB2*f[37]+1.224744871391589*EB0*f[37]
        +1.095445115010332*EB1*f[17]
        +1.224744871391589*EB2*f[10];
      V[45] += 2.449489742783178*EB1*f[44]+2.449489742783178*EB2*f[31]
        +2.738612787525831*EB0*f[31]
        +2.738612787525831*EB1*f[18];
      V[46] += 1.095445115010332*EB2*f[39]+1.224744871391589*EB0*f[39]
        +1.224744871391589*EB1*f[27];
      V[47] +=  1.095445115010332*EB2*f[42]+1.224744871391589*EB0*f[42]
        +1.224744871391589*EB1*f[30];
    }

    void vol(int nloop)
    {
      std::vector<double> f(48), F(48), fOut(48), EB(3);
      double w = 0.1;
      for (unsigned i=0; i<nloop; ++i)
      {
        calcDirectStreamVol(w, f, fOut);
        calcDirectForceVol(w, f, EB, fOut);
        calcDirectForceVol(w, f, EB, fOut);
        calcDirectForceVol(w, f, EB, fOut);
      }
    }
};

class Update1x3vMax {
  public:
    void calcDirectStreamVol(double w, const std::vector<double>& f, std::vector<double>& V)
    {
      using namespace std;
      double sq3 = sqrt(3), sq5 = sqrt(5), sq1_5 = 1/sqrt(5);
      V[0] += 0;
      V[1] += sq3*f[0]*w+f[2];
      V[5] += sq3*f[2]*w+(2*f[12])*sq1_5+f[0];
      V[6] += sq3*f[3]*w+f[7];
      V[8] += sq3*f[4]*w+f[9];
      V[11] += sq3*sq5*f[1]*w+sq5*f[5];
    }

    void calcDirectForceVol(double w, const std::vector<double>& f, const std::vector<double>& EB, std::vector<double>& V)
    { // this is only place holder for now
      register double EB0 = EB[0], EB1 = EB[1], EB2 = EB[2];
  
    }

    void vol(int nloop)
    {
      std::vector<double> f(15), F(15), fOut(15), EB(3);
      double w = 0.1;
      for (unsigned i=0; i<nloop; ++i)
      {
        calcDirectStreamVol(w, f, fOut);
        calcDirectForceVol(w, f, EB, fOut);
        calcDirectForceVol(w, f, EB, fOut);
        calcDirectForceVol(w, f, EB, fOut);
      }
    }
};


double run(NameValuePair& nvpair)
{
  int nloop = nvpair.getValue("nloop");
  int dim = nvpair.getValue("dim");
  int basis = nvpair.getValue("basis");
  clock_t t1 = clock();

  Update1x3vSerendip u1x3vs;

  if (dim == 4 && basis == 0) u1x3vs.vol(nloop);

  clock_t t2 = clock();
  return (double) (t2-t1)/CLOCKS_PER_SEC;
}

int
getNumDof(NameValuePair& nvpair) {
  int dim = nvpair.getValue("dim");
  int basis = nvpair.getValue("basis");  
  if (dim == 4 && basis == 0)
    return 48;
  else if (dim == 4 && basis == 1)
    return 15;
  else if (dim == 5 && basis == 0)
    return 112;
  else if (dim == 5 && basis == 1)
    return 21;
}

int
main (int argc, char **argv)
{
  if (argc != 2) {
    std::cout << "Usage::" << std::endl;
    std::cout << " fast-dg <input-file>" << std::endl;
    std::cout << "  input-file: Name of input file" << std::endl;

    exit(1);
  }
  std::string inFile(argv[1]); // input file
  NameValuePair nvpair(inFile);
  double tm = run(nvpair);

  int nloop = nvpair.getValue("nloop");
  double ndof = 1.0*getNumDof(nvpair);
  std::cout << "Update took " << tm << std::endl;
  std::cout << "Cost per DOF " << tm/nloop/ndof << std::endl;
  
  return 1;
}
