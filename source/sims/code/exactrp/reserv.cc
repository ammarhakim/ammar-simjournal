/** Compute exact solution to Euler reservoir problem. */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <limits>

#include "clipp.h"

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

/** Problem state */
class ProblemState {
  public:
/** Density, velocity, pressure, sound-speed on left */
    double dl, ul, pl, cl;
/** Density, velocity, pressure, sound-speed on right */
    double dr, ur, pr, cr;
/** Normalising constant */
    double pscale;
/** Lower and upper bounds */
    double lower, upper;
/** Domain length, time at which solution is needed */
    double domlen, tEnd;
/** Ratio of specific heat */
    double gas_gamma;
/** Location of discontinuity */
    double disLoc;
};

/** Solution of Riemann problem */
class Solution {
  public:
/**
 * @param ncells Number of cells in domain.
 */
    Solution(unsigned ncell)
      : ncell(ncell),
        density(ncell),
        velocity(ncell),
        pressure(ncell),
        internalEnergy(ncell)
    {
    }
    
/** Number of cells */
    unsigned ncell;
/** Density, velocity, pressure and internal-energy */
    std::vector<double> density, velocity, pressure, internalEnergy;
};

void
showInput(const ProblemState& ps) {
  std::cout << "gamma " << ps.gas_gamma << std::endl;

  std::cout << "dl " << ps.dl << std::endl;
  std::cout << "ul " << ps.ul << std::endl;
  std::cout << "pl " << ps.pl << std::endl;

  std::cout << "dr " << ps.dr << std::endl;
  std::cout << "ur " << ps.ur << std::endl;
  std::cout << "pr " << ps.pr << std::endl;
}

void
prefun(const ProblemState& ps, double& F, double& FD,
  double P, double DK, double PK, double CK) {

  double PRATIO, QRT, AK, BK;
  double gas_gamma = ps.gas_gamma;
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;

  if (P<=PK)
  {
// rarefaction wave
    PRATIO = P/PK;
    F = G4*CK*(pow(PRATIO,G1) - 1.0);
    FD = (1.0/(DK*CK))*pow(PRATIO,-G2);
  }
  else
  {
// shock wave
    AK  = G5/DK;
    BK  = G6*PK;
    QRT = sqrt(AK/(BK + P));
    F   = (P - PK)*QRT;
    FD  = (1.0 - 0.5*(P - PK)/(BK + P))*QRT;
  }
}

/**
 * @param ps Problem state
 * @return Guess for pressure in star region
 */
double
guessp(const ProblemState& ps) {
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  double gas_gamma = ps.gas_gamma;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;

  double CUP, GEL, GER, PM, PMAX, PMIN, PPV, PQ,
    PTL, PTR, QMAX, QUSER, UM;

  QUSER = 2.0;

  CUP  = 0.25*(DL + DR)*(CL + CR);
  PPV  = 0.5*(PL + PR) + 0.5*(UL - UR)*CUP;
  PPV  = std::max(0.0, PPV);
  PMIN = std::min(PL,  PR);
  PMAX = std::max(PL,  PR);
  QMAX = PMAX/PMIN;

  if ((QMAX <= QUSER) && ((PMIN<=PPV) && (PPV<=PMAX)))
  {
    //std::cout << "First if" << std::endl;
    PM = PPV;
  }
  else 
  {
    if (PPV<PMIN) 
    {
      //std::cout << "Second if" << std::endl;
      PQ  = pow(PL/PR,G1);
      UM  = (PQ*UL/CL + UR/CR + G4*(PQ - 1.0))/(PQ/CL + 1.0/CR);
      PTL = 1.0 + G7*(UL - UM)/CL;
      PTR = 1.0 + G7*(UM - UR)/CR;
      PM  = 0.5*(PL*pow(PTL,G3) + PR*pow(PTR,G3));
    }
    else 
    {
      //std::cout << "Second else" << std::endl;
      GEL = sqrt((G5/DL)/(G6*PL + PPV));
      GER = sqrt((G5/DR)/(G6*PR + PPV));
      PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER);
    }
  }
  return PM;
}

void
starpu(const ProblemState& ps, double &pm, double& um) {
  double CHANGE, FL, FLD, FR, FRD, P, POLD, 
    PSTART, TOLPRE, U, UDIFF, PSCALE;

  TOLPRE = 1.0e-6;
  unsigned NRITER = 20;

  double gas_gamma = ps.gas_gamma;
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;
  PSCALE = 1.0;

// compute initial guess
  PSTART = guessp(ps);
  //std::cout << "Initial guess " << PSTART << std::endl;

  POLD = PSTART;
  UDIFF = UR-UL;

  bool converged = false;
  for (unsigned i=1; i<NRITER; ++i) {
    prefun(ps, FL, FLD, POLD, DL, PL, CL);
    prefun(ps, FR, FRD, POLD, DR, PR, CR);
    P = POLD - (FL + FR + UDIFF)/(FLD + FRD);
    CHANGE = 2.0*fabs((P - POLD)/(P + POLD));
    //std::cout << "Iteration " << i << " Change " << CHANGE << std::endl;
    if (CHANGE<=TOLPRE) {
      converged = true;
      break;
    }
    if (P<0.0) P = TOLPRE;
    POLD = P;
  }
  if (converged)
    /*std::cout << "Newton iteration converged ..." << std::endl*/;
  else {
    std::cout << "Newton iteration did not converge" << std::endl;
    exit(1);
  }
// compute velocity in star region
  U = 0.5*(UL + UR + FR - FL);
  //std::cout << "Converged pressure " << P << std::endl;
  //std::cout << "Converged velocity " << U << std::endl;
  pm = P;
  um = U;
}

void
sample(const ProblemState& ps, double PM, double UM, double S,
  double& D, double& U, double& P) {
  double gas_gamma = ps.gas_gamma;
// compute constants related to gamma
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;

  double C, CML, CMR, PML, PMR, SHL, SHR, SL, SR, STL, STR;

  //std::cout << "S = " << S << std::endl;
  //std::cout << "UM = " << UM << std::endl;
  if (S<=UM)
  {
// left of contact discontinuity
    //std::cout << "Left of contact .." << std::endl;

    if (PM<=PL)
    {
      //std::cout << "Left rarefaction ...." << std::endl;

// left rarefaction
      SHL = UL-CL;
      if (S<=SHL)
      {
        //std::cout << "Left data state!" << std::endl;

// left data state
        D = DL;
        U = UL;
        P = PL;
      }
      else 
      {
        CML = CL*pow(PM/PL, G1);
        STL = UM-CML;

        if (S>STL)
        {
          //std::cout << "Star left state!" << std::endl;

// Star left state
          D = DL*pow(PM/PL, 1/GAMMA);
          U = UM;
          P = PM;
        } 
        else 
        {
          //std::cout << "Inside left fan!" << std::endl;
// inside left fan
          U = G5*(CL + G7*UL + S);
          C = G5*(CL + G7*(UL - S));
          D = DL*pow(C/CL, G4);
          P = PL*pow(C/CL, G3);
        }
      }
    }
    else
    {
      //std::cout << "Left shock ..." << std::endl;
// left shock
      PML = PM/PL;
      SL = UL - CL*sqrt(G2*PML + G1);

      if (S<=SL)
      {
        //std::cout << "Left data state!" << std::endl;
// point is left data state
        D = DL;
        U = UL;
        P = PL;
      }
      else
      {
        //std::cout << "Star left state!" << std::endl;
// point is star left state
        D = DL*(PML + G6)/(PML*G6 + 1.0);
        U = UM;
        P = PM;
      }
    }
  }
  else
  {
// right of contact discontinuity
    //std::cout << "Right of contact discontinuity ...." << std::endl;
    if (PM>PR)
    {
// right shock
      //std::cout << "Right shock ..." << std::endl;

      PMR = PM/PR;
      SR  = UR + CR*sqrt(G2*PMR + G1);

      if (S>=SR)
      {
        //std::cout << "Right data state!" << std::endl;
// right data state
        D = DR;
        U = UR;
        P = PR;
      }
      else
      {
        //std::cout << "Star right state!" << std::endl;
// right star state
        D = DR*(PMR + G6)/(PMR*G6 + 1.0);
        U = UM;
        P = PM;
      }
    }
    else
    {
      //std::cout << "Right rarefaction ..." << std::endl;
// right rarefaction
      SHR = UR+CR;

      if (S>=SHR)
      {
        //std::cout << "Right data state!" << std::endl;
// right data state
        D = DR;
        U = UR;
        P = PR;
      }
      else
      {
        CMR = CR*pow(PM/PR, G1);
        STR = UM + CMR;

        if (S<=STR)
        {
          //std::cout << "Star right state!" << std::endl;
// star right state
          D = DR*pow(PM/PR, 1.0/GAMMA);
          U = UM;
          P = PM;
        }
        else
        {
          //std::cout << "Inside left fan!" << std::endl;
// inside left fan
          U = G5*(-CR + G7*UR + S);
          C = G5*(CR - G7*(UR - S));
          D = DR*pow(C/CR, G4);
          P = PR*pow(C/CR, G3);
        }
      }
    }
  }
}

void
sampleWithVacuum(const ProblemState& ps, double S, double& D, double& U, double& P) {
  double gas_gamma = ps.gas_gamma;
// compute constants related to gamma
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;

// assume vacuum is to the right (JUST FOR NOW)
  double sstar = UL + 2*CL/G8;
  double relVel = UL-CL;
  if (S<=relVel)
  {
// left of rarefaction
    D = DL;
    U = UL;
    P = PL;
  }
  else
  {
    if ((relVel<S) && (S<=sstar))
    {
// inside rarefaction fan
      D = DL*pow(G5 + G6/CL*(UL-S), G4);
      U = G5*(CL + G7*UL + S);
      P = PL*pow(G5 + G6/CL*(UL-S), G3);
    }
    else
    {
// inside vacuum
      D = 0.0;
      U = 0.0;
      P = 0.0;
    }
  }
}

void
writeSolution(const ProblemState& ps, const std::string& outPrefix, const Solution& solution) {
  double xcell;
  double dx = ps.domlen/solution.ncell;

// write density
  std::string dhFileNm(outPrefix);
  dhFileNm += "-density.txt";
  std::ofstream dhFile(dhFileNm.c_str());
  dhFile << "## X [m] Density [kg/m^3]" << std::endl;
  for (unsigned i=0; i<solution.ncell; ++i) {
    xcell = ps.lower + (i+0.5)*dx;
    dhFile << xcell << " " << solution.density[i] << std::endl;
  }

// write pressure
  std::string prFileNm(outPrefix);
  prFileNm += "-pressure.txt";
  std::ofstream prFile(prFileNm.c_str());
  prFile << "## X [m] Pressure [Pa]" << std::endl;
  for (unsigned i=0; i<solution.ncell; ++i) {
    xcell = ps.lower + (i+0.5)*dx;
    prFile << xcell << " " << solution.pressure[i] << std::endl;
  }

// write velocity
  std::string ulFileNm(outPrefix);
  ulFileNm += "-velocity.txt";
  std::ofstream ulFile(ulFileNm.c_str());
  ulFile << "## X [m] Velocity [m/s]" << std::endl;
  for (unsigned i=0; i<solution.ncell; ++i) {
    xcell = ps.lower + (i+0.5)*dx;
    ulFile << xcell << " " << solution.velocity[i] << std::endl;
  }

// write internal energy
  std::string ieFileNm(outPrefix);
  ieFileNm += "-internal-energy.txt";
  std::ofstream ieFile(ieFileNm.c_str());
  ieFile << "## X [m] Internal-energy [J]" << std::endl;
  for (unsigned i=0; i<solution.ncell; ++i) {
    xcell = ps.lower + (i+0.5)*dx;
    ieFile << xcell << " " << solution.internalEnergy[i] << std::endl;
  }
}

void flux(double d, double u, double p, double g, std::vector<double>& f) {
  f[0] = d*u;
  f[1] = d*u*u + p;
  double E = p/(g-1)+0.5*d*u*u;
  f[2] = (E+p)*u;
}

void
numFlux(const ProblemState& ps, std::vector<double>& fl) {
  double g1, g2, g3, g4, g5, g6, g7, g8;

  double gas_gamma = ps.gas_gamma;
// compute constants related to gamma
  g1 = (gas_gamma - 1)/(2*gas_gamma);
  g2 = (gas_gamma + 1)/(2*gas_gamma);
  g3 = 2*gas_gamma/(gas_gamma - 1);
  g4 = 2/(gas_gamma - 1);
  g5 = 2/(gas_gamma + 1);
  g6 = (gas_gamma - 1)/(gas_gamma + 1);
  g7 = (gas_gamma - 1)/2;
  g8 = gas_gamma - 1;

// check if a vacuum is generated
  if (g4*(ps.cl+ps.cr) <= (ps.ur-ps.ul)) {
    std::cout << "Initial conditions will lead to vacuum. Aborting ..." << std::endl;
    exit(1);
  }

// compute pressure and velocity in "star" region
  double pm, um;
  starpu(ps, pm, um);

  double dsol, usol, psol;
  sample(ps, pm, um, 0.0, dsol, usol, psol);
  flux(dsol, usol, psol, gas_gamma, fl);
}

double
cmpStates(std::vector<double>& q1, std::vector<double>& q2) {
  double deld = (q1[0]-q2[0])/q1[0];
  double dele = (q1[2]-q2[2])/q1[2];
  return std::fabs(deld)+std::fabs(dele);
}

void
exactEulerReserv(const ProblemState& ps, unsigned maxIter) {
  ProblemState psl(ps), psr(ps);
    
  // this is an effective dt/dx computed from 1D stability condition c*dt/dx < 1
  double dtdx = 1/std::max(ps.cl+std::fabs(ps.ul), ps.cr+std::fabs(ps.ur));
  
  double gas_gamma = ps.gas_gamma;
  
  double dm = 0.5*(ps.dl+ps.dr);
  double um = 0.5*(ps.ul+ps.ur);
  double pm = 0.5*(ps.pl+ps.pr);

  // initial state (conserved variables)
  std::vector<double> qstar(3);
  qstar[0] = dm;
  qstar[1] = dm*um;
  qstar[2] = pm/(gas_gamma-1)+0.5*dm*um*um;

  std::vector<double> qnew(3), fl(3), fr(3);

  unsigned count = 0;
  while(1) {
    count = count+1;
    if (count>maxIter) break;
    
    // current intermediate state
    double dsm, usm, psm, csm;
    dsm = qstar[0];
    usm = qstar[1]/qstar[0];
    psm = (qstar[2]-0.5*dsm*usm*usm)*(gas_gamma-1);
    csm = std::sqrt(gas_gamma*psm/dsm);

    // copy into left/right state objects
    psl.dr = dsm;
    psl.ur = usm;
    psl.pr = psm;
    psl.cr = csm;

    psr.dl = dsm;
    psr.ul = usm;
    psr.pl = psm;
    psr.cl = csm;

    // compute numerical fluxes
    numFlux(psl, fl);
    numFlux(psr, fr);

    // update solution
    for (unsigned i=0; i<3; ++i)
      qnew[i] = qstar[i] - dtdx*(fr[i]-fl[i]);

    double err = cmpStates(qnew, qstar);
    for (unsigned i=0; i<3; ++i)
      qstar[i] = qnew[i];

    if (err<1e-6)  break;
  }
  if (count<maxIter)  {
    std::cout << count << " " << qnew[0] << " " << qnew[1] << " " << qnew[2] << std::endl;
  }
  else
    std::cout << "FAILED TO CONVERGE!" << std::endl;
}

void
runWithInputFile(const std::string& fileNm) {
  std::string inFile(fileNm);

  ProblemState ps;

// initialize problem state from input file
  NameValuePair nvpair(inFile);
  ps.dl = nvpair.getValue("dl");
  ps.ul = nvpair.getValue("ul");
  ps.pl = nvpair.getValue("pl");
  ps.dr = nvpair.getValue("dr");
  ps.ur = nvpair.getValue("ur");
  ps.pr = nvpair.getValue("pr");
  ps.pscale = 1.0;
  ps.disLoc = 0.0;
  ps.tEnd = 0.1;
  ps.gas_gamma = nvpair.getValue("gamma");

  ps.lower = -1.0;
  ps.upper = 1.0;
  ps.domlen = ps.upper - ps.lower;

// compute sound speeds in each region
  if (ps.dl != 0)
    ps.cl = std::sqrt(ps.gas_gamma*ps.pl/ps.dl);
  else
    ps.cl = 0.0;

  if (ps.dr != 0)
    ps.cr = std::sqrt(ps.gas_gamma*ps.pr/ps.dr);
  else
    ps.cr = 0.0;

  unsigned maxIter = 1000;
// maximum number of iterations
  if (nvpair.hasValue("maxIter"))
    maxIter = (unsigned) nvpair.getValue("maxIter");

// compute solution
  if ((ps.dl != 0.0) && (ps.dr != 0.0))
    exactEulerReserv(ps, maxIter);
  else if (ps.dr == 0)
    /*exactEulerReservWithVacuum(ps, solution)*/;
  else
    /*exactEulerReservWithVacuum(ps, solution)*/;  
}

void
runWithCmdArgs(unsigned nmax, double gg, double dl, double pl, double dr, double pr) {
  ProblemState ps;
    
  ps.dl = dl;
  ps.ul = 0.0;
  ps.pl = pl;
  ps.dr = dr;
  ps.ur = 0.0;
  ps.pr = pr;
  ps.pscale = 1.0;
  ps.disLoc = 0.0;
  ps.tEnd = 0.1;
  ps.gas_gamma = gg;

  ps.lower = -1.0;
  ps.upper = 1.0;
  ps.domlen = ps.upper - ps.lower;

// compute sound speeds in each region
  if (ps.dl != 0)
    ps.cl = std::sqrt(ps.gas_gamma*ps.pl/ps.dl);
  else
    ps.cl = 0.0;

  if (ps.dr != 0)
    ps.cr = std::sqrt(ps.gas_gamma*ps.pr/ps.dr);
  else
    ps.cr = 0.0;

  unsigned maxIter = nmax;

// compute solution
  if ((ps.dl != 0.0) && (ps.dr != 0.0))
    exactEulerReserv(ps, maxIter);
  else if (ps.dr == 0)
    /*exactEulerReservWithVacuum(ps, solution)*/;
  else
    /*exactEulerReservWithVacuum(ps, solution)*/;    
}

int
main(int argc, char **argv) {
  unsigned nmax = 1000;
  std::string inpFile("");
  double dl, pl, dr, pr, gasGamma;
  
  auto cli = (
    clipp::option("-i") & clipp::opt_value("Input file", inpFile),
    clipp::option("-n").doc("Max number of iteration") & clipp::value("nmax", nmax),
    clipp::option("-gamma").doc("Gas gamma") & clipp::value("gamma", gasGamma),
    clipp::option("-dl").doc("Left density") & clipp::value("dl", dl),
    clipp::option("-pl").doc("Left pressure") & clipp::value("pl", pl),
    clipp::option("-dr").doc("Right density") & clipp::value("dr", dr),
    clipp::option("-pr").doc("Right pressure") & clipp::value("pr", pr)
  );
  
  if (clipp::parse(argc, argv, cli)) {
    if (inpFile.compare("") > 0)
      runWithInputFile(inpFile);
    else
      runWithCmdArgs(nmax, gasGamma, dl, pl, dr, pr);

  } else {
    std::cout << clipp::make_man_page(cli, argv[0]) << std::endl;
  }

  return 1;
}
