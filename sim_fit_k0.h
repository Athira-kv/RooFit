#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <cmath>

#include <TROOT.h>

#include <TCanvas.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>

#include <RooGlobalFunc.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooAddition.h>
#include <RooCategory.h>
#include <RooChebychev.h>
#include <RooConstVar.h>
#include <RooDataHist.h>
#include <RooFFTConvPdf.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooGlobalFunc.h>
//#include <RooMinuit.h>
#include <RooPlot.h>
#include <RooPolynomial.h>
#include <RooRealVar.h>
//#include <RooRelBreitWigner.h>
#include <RooSimultaneous.h>
#include <RooVoigtian.h>

using std::cout;
using std::cerr;
using std::endl;
using std::map;
using std::stringstream;
using std::vector;
using std::pair;
using std::string;
using std::make_pair;
using std::ios;
using std::setw;
using std::setprecision;


using namespace RooFit ;
using namespace std;

const double M_K =  0.493677;
const double M_phi =  1.019461;

