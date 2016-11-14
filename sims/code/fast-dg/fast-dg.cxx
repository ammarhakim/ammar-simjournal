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

class Update1x1vSerendip {
  public:
    void calcDirectStreamVol(double w, const std::vector<double>& f, std::vector<double>& V)
    {
      double vFac = 1.5, wFac = 1.5*w;
      
      V[1] += 1.732050807568877*f[0]*wFac+f[2]*vFac; 
      V[3] += 1.732050807568877*f[2]*wFac+0.8944271909999159*f[5]*vFac+f[0]*vFac; 
      V[4] += 3.872983346207417*f[1]*wFac+2.23606797749979*f[3]*vFac; 
      V[6] += 3.872983346207417*f[3]*wFac+2.0*f[7]*vFac+2.23606797749979*f[1]*vFac; 
      V[7] += 1.732050807568877*f[5]*wFac+0.8944271909999159*f[2]*vFac;
    }

    void calcDirectForceVol(double w, const std::vector<double>& f, const std::vector<double>& E, std::vector<double>& V)
    {
      double exFac = 1.5;
      double E0 = E[0]*exFac, E1 = E[1]*exFac, E2 = E[2]*exFac; 
      V[2] += 1.224744871391589*f[4]*E2+1.224744871391589*f[1]*E1+1.224744871391589*f[0]*E0; 
      V[3] += 1.095445115010332*f[1]*E2+1.095445115010332*f[4]*E1+1.224744871391589*f[0]*E1+1.224744871391589*f[1]*E0; 
      V[5] += 2.738612787525831*f[6]*E2+2.738612787525831*f[3]*E1+2.738612787525831*f[2]*E0; 
      V[6] += 0.7824607964359517*f[4]*E2+1.224744871391589*f[0]*E2+1.095445115010332*f[1]*E1+1.224744871391589*f[4]*E0; 
      V[7] += 2.449489742783178*f[3]*E2+2.449489742783178*f[6]*E1+2.738612787525831*f[2]*E1+2.738612787525831*f[3]*E0; 
    }

    void vol(int nloop)
    {
      std::vector<double> f(48), F(48), fOut(48), EB(3);
      double w = 0.1;
      for (unsigned i=0; i<nloop; ++i)
      {
        calcDirectStreamVol(w, f, fOut);
        calcDirectForceVol(w, f, EB, fOut);
      }
    }
};

class Update1x3vSerendip {
  public:
    void calcDirectStreamVol(double w, const std::vector<double>& f, std::vector<double>& V)
    {
      double vFac = 1.5, wFac = 1.5*w;
      V[1] += 1.732050807568877*f[0]*wFac+f[2]*vFac; 
      V[5] += 1.732050807568877*f[2]*wFac+0.8944271909999159*f[12]*vFac+f[0]*vFac; 
      V[6] += 1.732050807568877*f[3]*wFac+f[7]*vFac; 
      V[8] += 1.732050807568877*f[4]*wFac+f[9]*vFac; 
      V[11] += 3.872983346207417*f[1]*wFac+2.23606797749979*f[5]*vFac; 
      V[15] += 1.732050807568877*f[7]*wFac+0.8944271909999159*f[22]*vFac+f[3]*vFac; 
      V[16] += 1.732050807568877*f[9]*wFac+0.8944271909999159*f[26]*vFac+f[4]*vFac; 
      V[17] += 1.732050807568877*f[10]*wFac+f[18]*vFac; 
      V[19] += 3.872983346207417*f[5]*wFac+2.0*f[20]*vFac+2.23606797749979*f[1]*vFac; 
      V[20] += 1.732050807568877*f[12]*wFac+0.8944271909999159*f[2]*vFac; 
      V[21] += 3.872983346207417*f[6]*wFac+2.23606797749979*f[15]*vFac; 
      V[23] += 1.732050807568877*f[13]*wFac+f[24]*vFac; 
      V[25] += 3.872983346207417*f[8]*wFac+2.23606797749979*f[16]*vFac; 
      V[28] += 1.732050807568877*f[14]*wFac+f[29]*vFac; 
      V[31] += 1.732050807568877*f[18]*wFac+0.8944271909999159*f[38]*vFac+f[10]*vFac; 
      V[32] += 3.872983346207417*f[15]*wFac+2.0*f[33]*vFac+2.23606797749979*f[6]*vFac; 
      V[33] += 1.732050807568877*f[22]*wFac+0.8944271909999159*f[7]*vFac; 
      V[34] += 1.732050807568877*f[24]*wFac+f[13]*vFac; 
      V[35] += 3.872983346207417*f[16]*wFac+2.0*f[36]*vFac+2.23606797749979*f[8]*vFac; 
      V[36] += 1.732050807568877*f[26]*wFac+0.8944271909999159*f[9]*vFac; 
      V[37] += 3.872983346207417*f[17]*wFac+2.23606797749979*f[31]*vFac; 
      V[39] += 1.732050807568877*f[27]*wFac+f[40]*vFac; 
      V[41] += 1.732050807568877*f[29]*wFac+f[14]*vFac; 
      V[42] += 1.732050807568877*f[30]*wFac+f[43]*vFac; 
      V[44] += 3.872983346207417*f[31]*wFac+2.0*f[45]*vFac+2.23606797749979*f[17]*vFac; 
      V[45] += 1.732050807568877*f[38]*wFac+0.8944271909999159*f[18]*vFac; 
      V[46] += 1.732050807568877*f[40]*wFac+f[27]*vFac; 
      V[47] += 1.732050807568877*f[43]*wFac+f[30]*vFac;       
    }

    void calcDirectForceVol(double w, const std::vector<double>& f, const std::vector<double>& E, std::vector<double>& V)
    {
      double exFac = 1.5;      
      double E0 = E[0]*exFac, E1 = E[1]*exFac, E2 = E[2]*exFac; 
      V[2] += 1.224744871391589*f[11]*E2+1.224744871391589*f[1]*E1+1.224744871391589*f[0]*E0; 
      V[5] += 1.095445115010332*f[1]*E2+1.095445115010332*f[11]*E1+1.224744871391589*f[0]*E1+1.224744871391589*f[1]*E0; 
      V[7] += 1.224744871391589*f[21]*E2+1.224744871391589*f[6]*E1+1.224744871391589*f[3]*E0; 
      V[9] += 1.224744871391589*f[25]*E2+1.224744871391589*f[8]*E1+1.224744871391589*f[4]*E0; 
      V[12] += 2.738612787525831*f[19]*E2+2.738612787525831*f[5]*E1+2.738612787525831*f[2]*E0; 
      V[15] += 1.095445115010332*f[6]*E2+1.095445115010332*f[21]*E1+1.224744871391589*f[3]*E1+1.224744871391589*f[6]*E0; 
      V[16] += 1.095445115010332*f[8]*E2+1.095445115010332*f[25]*E1+1.224744871391589*f[4]*E1+1.224744871391589*f[8]*E0; 
      V[18] += 1.224744871391589*f[37]*E2+1.224744871391589*f[17]*E1+1.224744871391589*f[10]*E0; 
      V[19] += 0.7824607964359517*f[11]*E2+1.224744871391589*f[0]*E2+1.095445115010332*f[1]*E1+1.224744871391589*f[11]*E0; 
      V[20] += 2.449489742783178*f[5]*E2+2.449489742783178*f[19]*E1+2.738612787525831*f[2]*E1+2.738612787525831*f[5]*E0; 
      V[22] += 2.738612787525831*f[32]*E2+2.738612787525831*f[15]*E1+2.738612787525831*f[7]*E0; 
      V[24] += 1.224744871391589*f[23]*E1+1.224744871391589*f[13]*E0; 
      V[26] += 2.738612787525831*f[35]*E2+2.738612787525831*f[16]*E1+2.738612787525831*f[9]*E0; 
      V[29] += 1.224744871391589*f[28]*E1+1.224744871391589*f[14]*E0; 
      V[31] += 1.095445115010332*f[17]*E2+1.095445115010332*f[37]*E1+1.224744871391589*f[10]*E1+1.224744871391589*f[17]*E0; 
      V[32] += 0.7824607964359517*f[21]*E2+1.224744871391589*f[3]*E2+1.095445115010332*f[6]*E1+1.224744871391589*f[21]*E0; 
      V[33] += 2.449489742783178*f[15]*E2+2.449489742783178*f[32]*E1+2.738612787525831*f[7]*E1+2.738612787525831*f[15]*E0; 
      V[34] += 1.095445115010332*f[23]*E2+1.224744871391589*f[13]*E1+1.224744871391589*f[23]*E0; 
      V[35] += 0.7824607964359517*f[25]*E2+1.224744871391589*f[4]*E2+1.095445115010332*f[8]*E1+1.224744871391589*f[25]*E0; 
      V[36] += 2.449489742783178*f[16]*E2+2.449489742783178*f[35]*E1+2.738612787525831*f[9]*E1+2.738612787525831*f[16]*E0; 
      V[38] += 2.738612787525831*f[44]*E2+2.738612787525831*f[31]*E1+2.738612787525831*f[18]*E0; 
      V[40] += 1.224744871391589*f[39]*E1+1.224744871391589*f[27]*E0; 
      V[41] += 1.095445115010332*f[28]*E2+1.224744871391589*f[14]*E1+1.224744871391589*f[28]*E0; 
      V[43] += 1.224744871391589*f[42]*E1+1.224744871391589*f[30]*E0; 
      V[44] += 0.7824607964359517*f[37]*E2+1.224744871391589*f[10]*E2+1.095445115010332*f[17]*E1+1.224744871391589*f[37]*E0; 
      V[45] += 2.449489742783178*f[31]*E2+2.449489742783178*f[44]*E1+2.738612787525831*f[18]*E1+2.738612787525831*f[31]*E0; 
      V[46] += 1.095445115010332*f[39]*E2+1.224744871391589*f[27]*E1+1.224744871391589*f[39]*E0; 
      V[47] += 1.095445115010332*f[42]*E2+1.224744871391589*f[30]*E1+1.224744871391589*f[42]*E0;       
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
      double vFac = 1.5, wFac = 1.5*w;
      V[1] += 1.732050807568877*f[0]*wFac+f[2]*vFac; 
      V[5] += 1.732050807568877*f[2]*wFac+0.8944271909999159*f[12]*vFac+f[0]*vFac; 
      V[6] += 1.732050807568877*f[3]*wFac+f[7]*vFac; 
      V[8] += 1.732050807568877*f[4]*wFac+f[9]*vFac; 
      V[11] += 3.872983346207417*f[1]*wFac+2.23606797749979*f[5]*vFac;       
    }

    void calcDirectForceVol(double w, const std::vector<double>& f, const std::vector<double>& E, std::vector<double>& V)
    {
      double exFac = 1.5;            
      double E0 = E[0]*exFac, E1 = E[1]*exFac, E2 = E[2]*exFac; 
      V[2] += 1.224744871391589*f[11]*E2+1.224744871391589*f[1]*E1+1.224744871391589*f[0]*E0; 
      V[5] += 1.095445115010332*f[1]*E2+1.095445115010332*f[11]*E1+1.224744871391589*f[0]*E1+1.224744871391589*f[1]*E0; 
      V[7] += 1.224744871391589*f[6]*E1+1.224744871391589*f[3]*E0; 
      V[9] += 1.224744871391589*f[8]*E1+1.224744871391589*f[4]*E0; 
      V[12] += 2.738612787525831*f[5]*E1+2.738612787525831*f[2]*E0; 
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

class Update2x3vMax {
  public:
    void calcDirectStreamVol_X(double w, const std::vector<double>& f, std::vector<double>& V)
    {
      double vFac = 1.5, wFac = 1.5*w;
      V[1] += 1.732050807568877*f[0]*wFac+f[3]*vFac; 
      V[6] += 1.732050807568877*f[2]*wFac+f[8]*vFac; 
      V[7] += 1.732050807568877*f[3]*wFac+0.8944271909999159*f[18]*vFac+f[0]*vFac; 
      V[9] += 1.732050807568877*f[4]*wFac+f[11]*vFac; 
      V[12] += 1.732050807568877*f[5]*wFac+f[14]*vFac; 
      V[16] += 3.872983346207417*f[1]*wFac+2.23606797749979*f[7]*vFac; 
      
    }
    void calcDirectStreamVol_Y(double w, const std::vector<double>& f, std::vector<double>& V)
    {
      double vFac = 1.5, wFac = 1.5*w;
      V[2] += 1.732050807568877*f[0]*wFac+f[4]*vFac; 
      V[6] += 1.732050807568877*f[1]*wFac+f[9]*vFac; 
      V[8] += 1.732050807568877*f[3]*wFac+f[11]*vFac; 
      V[10] += 1.732050807568877*f[4]*wFac+0.8944271909999159*f[19]*vFac+f[0]*vFac; 
      V[13] += 1.732050807568877*f[5]*wFac+f[15]*vFac; 
      V[17] += 3.872983346207417*f[2]*wFac+2.23606797749979*f[10]*vFac;       
    }    

    void calcDirectForceVol_VX(double w, const std::vector<double>& f, const std::vector<double>& E, std::vector<double>& V)
    {
      double exFac = 1.5;
      double E0 = E[0]*exFac, E1 = E[1]*exFac, E2 = E[2]*exFac, E3 = E[3]*exFac, E4 = E[4]*exFac, E5 = E[5]*exFac; 
      V[3] += 0.8660254037844386*f[17]*E5+0.8660254037844386*f[16]*E4+0.8660254037844386*f[6]*E3+0.8660254037844386*f[2]*E2+0.8660254037844386*f[1]*E1+0.8660254037844386*f[0]*E0; 
      V[7] += 0.7745966692414833*f[1]*E4+0.8660254037844386*f[2]*E3+0.8660254037844386*f[6]*E2+0.7745966692414833*f[16]*E1+0.8660254037844386*f[0]*E1+0.8660254037844386*f[1]*E0; 
      V[8] += 0.7745966692414833*f[2]*E5+0.8660254037844386*f[1]*E3+0.7745966692414833*f[17]*E2+0.8660254037844386*f[0]*E2+0.8660254037844386*f[6]*E1+0.8660254037844386*f[2]*E0; 
      V[11] += 0.8660254037844386*f[10]*E2+0.8660254037844386*f[9]*E1+0.8660254037844386*f[4]*E0; 
      V[14] += 0.8660254037844386*f[13]*E2+0.8660254037844386*f[12]*E1+0.8660254037844386*f[5]*E0; 
      V[18] += 1.936491673103709*f[8]*E2+1.936491673103709*f[7]*E1+1.936491673103709*f[3]*E0;       
    }

    void calcDirectForceVol_VY(double w, const std::vector<double>& f, const std::vector<double>& E, std::vector<double>& V)
    {
      double exFac = 1.5;
      double E0 = E[0]*exFac, E1 = E[1]*exFac, E2 = E[2]*exFac, E3 = E[3]*exFac, E4 = E[4]*exFac, E5 = E[5]*exFac; 
      V[4] += 0.8660254037844386*f[17]*E5+0.8660254037844386*f[16]*E4+0.8660254037844386*f[6]*E3+0.8660254037844386*f[2]*E2+0.8660254037844386*f[1]*E1+0.8660254037844386*f[0]*E0; 
      V[9] += 0.7745966692414833*f[1]*E4+0.8660254037844386*f[2]*E3+0.8660254037844386*f[6]*E2+0.7745966692414833*f[16]*E1+0.8660254037844386*f[0]*E1+0.8660254037844386*f[1]*E0; 
      V[10] += 0.7745966692414833*f[2]*E5+0.8660254037844386*f[1]*E3+0.7745966692414833*f[17]*E2+0.8660254037844386*f[0]*E2+0.8660254037844386*f[6]*E1+0.8660254037844386*f[2]*E0; 
      V[11] += 0.8660254037844386*f[8]*E2+0.8660254037844386*f[7]*E1+0.8660254037844386*f[3]*E0; 
      V[15] += 0.8660254037844386*f[13]*E2+0.8660254037844386*f[12]*E1+0.8660254037844386*f[5]*E0; 
      V[19] += 1.936491673103709*f[10]*E2+1.936491673103709*f[9]*E1+1.936491673103709*f[4]*E0; 
    }

    void vol(int nloop)
    {
      std::vector<double> f(21), F(21), fOut(21), EB(6);
      double w = 0.1;
      for (unsigned i=0; i<nloop; ++i)
      {
        calcDirectStreamVol_X(w, f, fOut);
        calcDirectStreamVol_Y(w, f, fOut);
        calcDirectForceVol_VX(w, f, EB, fOut);
        calcDirectForceVol_VY(w, f, EB, fOut);
        calcDirectForceVol_VY(w, f, EB, fOut);
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
  Update1x3vMax u1x3vm;
  Update1x1vSerendip u1x1vs;
  Update2x3vMax u2x3vm;

  if (dim == 2 && basis == 0)
    u1x1vs.vol(nloop);  
  else if (dim == 4 && basis == 0)
    u1x3vs.vol(nloop);
  else if (dim == 4 && basis == 1)
    u1x3vm.vol(nloop);
  else if (dim == 5 && basis == 1)
    u2x3vm.vol(nloop);

  clock_t t2 = clock();
  return (double) (t2-t1)/CLOCKS_PER_SEC;
}

int
getNumDof(NameValuePair& nvpair) {
  int dim = nvpair.getValue("dim");
  int basis = nvpair.getValue("basis");
  if (dim == 2 && basis == 0)
    return 8;
  if (dim == 2 && basis == 1)
    return 6;
  else if (dim == 4 && basis == 0)
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
