#ifndef HLTv2_h
#define HLTv2_h 1

#include <typeinfo>

// #include "DataFormats/HcalDetId/interface/HcalDetId.h"
//#include "HcalPulseShapes.h"
#include "HcalTimeSlew.h"
#include "PedestalSub.h"
// #include "CalibFormats/HcalObjects/interface/HcalCoder.h"
// #include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
//#include "HybridMinimizer.h"

//#include <TMinuit.h>

//#include <TH1.h>
//#include "Minuit2/FCNBase.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"

class HLTv2 {
 public:
  enum NegStrategy {DoNothing=0, ReqPositive=1};  
  HLTv2();
  ~HLTv2();

  void Init(HcalTimeSlew::ParaSource tsParam, HcalTimeSlew::BiasSetting bias, NegStrategy nStrat, PedestalSub pedSubFxn_);

  // This is the CMSSW Implementation of the apply function
  //void apply(const CaloSamples & cs, const std::vector<int> & capidvec, const HcalCalibrations & calibs, std::vector<double> & correctedOutput) const;
  // This is the edited implementation for our standalone test code

  void apply(const std::vector<double> & inputCharge, const std::vector<double> & inputPedestal, std::vector<double> & HLTOutput) const;
  void applyXM(const std::vector<double> & inputCharge, const std::vector<double> & inputPedestal, std::vector<double> & HLTOutput) const;
  void getLandauFrac(Float_t tStart, Float_t tEnd, Float_t &sum) const;

  //void SolveEquations(Double_t *TS, Double_t *par, Double_t *fit) const;
  double Det2(double *b, double *c) const;
  double Det3(double *a, double *b, double *c) const;
  void PulseFraction(Double_t fC, Double_t *TS46) const;

 private:
  HcalTimeSlew::ParaSource fTimeSlew;
  HcalTimeSlew::BiasSetting fTimeSlewBias;
  NegStrategy fNegStrat;
  PedestalSub fPedestalSubFxn_;
  
};

#endif // HLTAnalyzer_h
