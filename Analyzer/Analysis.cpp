#define Analysis_cxx

#include "Analysis.h"
#include "readparameters/readparameters.h"
#include "TProfile.h"
#include "TProfile2D.h"


#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include <sys/stat.h>
#include "TLine.h"
#include "TLegend.h"

using namespace std;

double round_nplaces(double value, int to){
  int places = 1, whole = value;
  for(int i = 0; i < to; i++) places *= 10;
  value -= whole; //leave decimals
  value *= places; //0.1234 -> 123.4
  value = round(value);//123.4 -> 123
  value /= places; //123 -> .123
  value += whole; //bring the whole value back
  return value;
}

string int2string(int i){
  stringstream ss;
  string ret;
  ss<<i;
  ss>>ret;
  return ret;
}

// Ugly hack to apply energy corrections to some HB- cells
double eCorr(int ieta, int iphi, double energy) {
// return energy correction factor for HBM channels 
// iphi=6 ieta=(-1,-15) and iphi=32 ieta=(-1,-7)
// I.Vodopianov 28 Feb. 2011
  static const float low32[7]  = {0.741,0.721,0.730,0.698,0.708,0.751,0.861};
  static const float high32[7] = {0.973,0.925,0.900,0.897,0.950,0.935,1};
  static const float low6[15]  = {0.635,0.623,0.670,0.633,0.644,0.648,0.600,
				  0.570,0.595,0.554,0.505,0.513,0.515,0.561,0.579};
  static const float high6[15] = {0.875,0.937,0.942,0.900,0.922,0.925,0.901,
				  0.850,0.852,0.818,0.731,0.717,0.782,0.853,0.778};
  
  double slope, mid, en;
  double corr = 1.0;

  if (!(iphi==6 && ieta<0 && ieta>-16) && !(iphi==32 && ieta<0 && ieta>-8)) 
    return corr;

  int jeta = -ieta-1;
  double xeta = (double) ieta;
  if (energy > 0.) en=energy;
  else en = 0.;

  if (iphi == 32) {
    slope = 0.2272;
    mid = 17.14 + 0.7147*xeta;
    if (en > 100.) corr = high32[jeta];
    else corr = low32[jeta]+(high32[jeta]-low32[jeta])/(1.0+exp(-(en-mid)*slope));
  }
  else if (iphi == 6) {
    slope = 0.1956;
    mid = 15.96 + 0.3075*xeta;
    if (en > 100.0) corr = high6[jeta];
    else corr = low6[jeta]+(high6[jeta]-low6[jeta])/(1.0+exp(-(en-mid)*slope));
  }

  return corr;
}

int main(int argc, char **argv)
{
  int ret=0;
  if (argc!=2) {
    cerr<<"Usage: ./Analysis <paramfile>"<<endl;
    ret=1;
  } else {

    readparameters rp(argv[1]);
    TChain* ch = new TChain("ExportTree/HcalTree");

    string filelistname;

    filelistname=rp.get<string>((string("in_filelist")).c_str());
    string line;
    ifstream filelist(filelistname.c_str());
    if (filelist.fail()) { //catch
      cerr << "\nERROR: Could not open " << filelistname << endl;
      exit(1);
    }
    while (getline(filelist,line)) {
      ch->Add(line.c_str());
    }

    Analysis Ana25ns((TTree*) ch);

    Ana25ns.Init(argv[1], (TTree*) ch);
    Ana25ns.Process();
    Ana25ns.Finish();
  }
 return ret;
}

Analysis::~Analysis() {
}

Analysis::Analysis(TTree *tree):analysistree(tree){};

void Analysis::Init(char* paramfile, TTree* tree)
{
  try {
    readparameters rp(paramfile);
    try {Output_File=rp.get<string>("Output_File");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Entries=rp.get<int>("Entries");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Plot_Dir=rp.get<string>("Plot_Dir");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Region=rp.get<int>("Region");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
  } 
  catch (exception& e) {cerr<<e.what()<<endl;} 

  cout << "Running on ";
  if (Entries==-1) cout << "all events." << endl;
  else cout << "first " << Entries << " events. " << endl;

  cout << "Output ROOT file: " << Output_File << endl;

  cout << "Using channels in ";
  if (Region==All) { cout << "all HCAL regions. " << endl; }
  else if (Region==Barrel) { cout << "HCAL barrel. " << endl; }
  else if (Region==Endcap) { cout << "HCAL endcap. " << endl; }

  int check=mkdir(Plot_Dir.c_str(),755);
  if (!check) {
    cout << "Saving files to: " << Plot_Dir << endl;
  }
  else {
    cout << "Double check your plot directory exists! Something funny happened." << endl;
    //exit(1);
  }

  return; 
}

void Analysis::Process() {
  if (fChain == 0) return;

  if (Entries==-1) Entries=fChain->GetEntries();

  fout = new TFile(Output_File.c_str(), "RECREATE");

  tout = new TTree("Events", "Events");
  
  DoHlt();
  //DrawPulses();
}

void Analysis::DrawPulses() {

  float ts1, ts2, ts3, thpd, tpre, wd1, wd2, wd3;

  HcalPulseShapes thePulses_;
  
  Int_t TS[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  Double_t pulse[10];

  tout->Branch("TS",    &TS,    "TS[10]/I");
  tout->Branch("pulse", &pulse, "pulse[10]/D");
  tout->Branch("ts1", &ts1, "ts1/F");
  tout->Branch("ts2", &ts2, "ts2/F");
  tout->Branch("ts3", &ts3, "ts3/F");
  tout->Branch("thpd", &thpd, "thpd/F");
  tout->Branch("tpre", &tpre, "tpre/F");
  tout->Branch("wd1", &wd1, "wd1/F");
  tout->Branch("wd2", &wd2, "wd2/F");
  tout->Branch("wd3", &wd3, "wd3/F");

  //112.5
  ts1=8. ; ts2=19. ; ts3=29.3; tpre=9.0; wd1=2.0; wd2=0.7; wd3=0.32;

  for (Int_t k=0; k<20; k++) {

    thpd=1.0+k*0.5;

    for (UInt_t i=0; i<10; i++) {pulse[i]=0;}
    thePulses_.computeHPDShape(ts1,ts2,ts3,thpd,tpre,wd1,wd2,wd3, thePulses_.hpdShape_);
    for (UInt_t i=0; i<10; i++) {
      for (UInt_t j=0; j<25; j++) {
	pulse[i]+=thePulses_.hpdShape_(25*i+j-92.5);
      }
    }

    tout->Fill();
  }

}

////// Pedestal
/// https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoLocalCalo/HcalRecAlgos/src/PedestalSub.cc

///// reference HLT M3
/// https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoLocalCalo/HcalRecAlgos/src/HcalDeterministicFit.cc

////// configure Method2
///  https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoLocalCalo/HcalRecProducers/python/HcalHitReconstructor_hbhe_cfi.py

////// class that define Method2
///  https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoLocalCalo/HcalRecAlgos/src/PulseShapeFitOOTPileupCorrection.cc

void Analysis::DoHlt() {

  // these are useful to suppress the verbosity of the minimizer
  //  extern int gErrorIgnoreLevel;
  //  gErrorIgnoreLevel = 1001;
  extern int gErrorIgnoreLevel; //miao's
  gErrorIgnoreLevel = 1001;  //miao's

  //=========================================================================      
  // These are the values we should set for Method 2 config with the python files
  // Don't currently have a setup to input with python but we shouldn't have to
  // change these values for our tests
  
  // --------------------------------------------------------------------
/*
  bool iPedestalConstraint = true;
  bool iTimeConstraint = true; // updated
  bool iAddPulseJitter = false;
  bool iUnConstrainedFit = false;
  bool iApplyTimeSlew = true;
  double iTS4Min = 0.; // updated
  double iTS4Max = 100.; // updated 
  double iPulseJitter = 1.;
  double iTimeMean = 0;
  double iTimeSig = 5.;
  double iPedMean = 0.;
  double iPedSig = 0.5;
  double iNoise = 1.;
  double iTMin = -12.5; // updated
  double iTMax = 12.5; // updated
  double its3Chi2 = 5.; // updated , not used
  double its4Chi2 = 15.; // updated 
  double its345Chi2 = 100.; // updated , not used
  double iChargeThreshold = 6.;
  int iFitTimes = 1;
*/

//begin miao's settings
  bool iPedestalConstraint = true;
  bool iTimeConstraint = true; // for M2
  bool iTimeConstraintLAG = false; // for LAG
  bool iAddPulseJitter = false;
  bool iUnConstrainedFit = false;
  bool iApplyTimeSlew = true;  // for M2
  bool iApplyTimeSlewLAG = false;  // for LAG
  double iTS4Min = 0.; // updated
  double iTS4Max = 100.; // updated 
  double iPulseJitter = 1.;
  double iTimeMean = 0;
  double iTimeSig = 5.;  //for M2
  double iTimeSigLAG = -1.0; //for LAG
  double iPedMean = 0.;
  double iPedSig = 0.5;
  double iNoise = 1.;
  double iTMin = -12.5; // updated
  double iTMax = 12.5; // updated
  double its3Chi2 = 5.; // updated , not used
  double its4Chi2 = 15.; // updated 
  double its345Chi2 = 100.; // updated , not used
  double iChargeThreshold = 6.;
  int iFitTimes = 1;
  int iTemptypeM2 = 0;
  int iTemptypeLAG = 1; 
//end miao's setting
   
  //========================================================================
  // Set the Method 2 Parameters here
  psFitOOTpuCorr_->setPUParams(iPedestalConstraint,iTimeConstraint,iAddPulseJitter,iUnConstrainedFit,
			       iApplyTimeSlew,iTS4Min, iTS4Max, iPulseJitter,iTimeMean,iTimeSig,
			       iPedMean,iPedSig,iNoise,iTMin,iTMax,its3Chi2,its4Chi2,its345Chi2,
			       iChargeThreshold,HcalTimeSlew::Medium, iFitTimes, iTemptypeM2);
  
  // Now set the Pulse shape type
  psFitOOTpuCorr_->setPulseShapeTemplate(theHcalPulseShapes_.getShape(105));

//begin miao's LAG para settings
  psFitOOTpuCorrLAG_->setPUParams(iPedestalConstraint,iTimeConstraintLAG,iAddPulseJitter,iUnConstrainedFit,
                               iApplyTimeSlewLAG,iTS4Min, iTS4Max, iPulseJitter,iTimeMean,iTimeSigLAG,
                               iPedMean,iPedSig,iNoise,iTMin,iTMax,its3Chi2,its4Chi2,its345Chi2,
                               iChargeThreshold,HcalTimeSlew::Medium, iFitTimes, iTemptypeLAG);
//end miao's LAG para settings


  Int_t iphi, ieta, depth;
  Int_t TS[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  Double_t Pulse[10];
  Double_t Ped[10];

  Double_t m2Charge, m2Time, m2Ped, m2Chi, m2Status;
  Double_t m2Pulse[10];

  Double_t LAGCharge, LAGTime, LAGPed, LAGChi, LAGStatus;   //miao's
  Double_t LAGPulse[10];                                    //miao's

  Double_t wTime;

  tout->Branch("iphi",  &iphi,  "iphi/I");
  tout->Branch("ieta",  &ieta,  "ieta/I");
  tout->Branch("depth", &depth, "depth/I");
  tout->Branch("TS",    &TS,    "TS[10]/I");
  tout->Branch("pulse", &Pulse, "pulse[10]/D");
  tout->Branch("ped",   &Ped,   "ped[10]/D");

  tout->Branch("m2Charge", &m2Charge, "m2Charge/D");
  tout->Branch("m2Time",   &m2Time,   "m2Time/D");
  tout->Branch("m2Ped",    &m2Ped,    "m2Ped/D");
  tout->Branch("m2Chi",    &m2Chi,    "m2Chi/D");
  tout->Branch("m2Status", &m2Status, "m2Status/D");
  tout->Branch("m2Pulse",  &m2Pulse,  "m2Pulse[10]/D");

//begin miao's
  tout->Branch("LAGCharge", &LAGCharge, "LAGCharge/D");
  tout->Branch("LAGTime",   &LAGTime,   "LAGTime/D");
  tout->Branch("LAGPed",    &LAGPed,    "LAGPed/D");
  tout->Branch("LAGChi",    &LAGChi,    "LAGChi/D");
  tout->Branch("LAGStatus", &LAGStatus, "LAGStatus/D");
  tout->Branch("LAGPulse",  &LAGPulse,  "LAGPulse[10]/D");
//end miao's



  tout->Branch("wTime", &wTime, "wTime/D");

  TProfile* chi2_M2  = new TProfile("Chi2_M2","Chi2_M2",58,15.,605.); //miao's
  TProfile* chi2_LAG = new TProfile("Chi2_LAG","Chi2_LAG",58,15.,605.); //miao's
  TProfile* LAG_M2   = new TProfile("LAG_M2","LAG_M2",100,-0.5,100.5); //miao's
  TProfile2D *LAG_M2_Ch  = new TProfile2D("LAG_M2_Ch","LAG_M2_Ch",100,-0.5,100.5,100,-0.5,100.5,"s"); //miao's

  //Loop over all events
  for (int jentry=0; jentry<Entries;jentry++) {

    fChain->GetEntry(jentry);
    if(jentry%1000==0) cout <<"Analyzed entry "<< jentry <<"/"<< Entries << endl;
    cout <<"Analyzed entry "<< jentry <<"/"<< Entries << endl;
    if(jentry>1000) return;

    for (int j = 0; j < (int)PulseCount; j++) {

      if (IEta[j]>14 && Region==Barrel) continue;
      if (IEta[j]<19 && Region==Endcap) continue;
      if (IEta[j]>27 && Region==Endcap) continue; // remove the bad tower 29 

      //      if (IEta[j]>14 && Region==Barrel) continue;
      //      if (IEta[j]<19 && Region==Endcap) continue;

      std::vector<double> inputCaloSample, inputPedestal, inputGain;                                                           
      std::vector<double> offlineAns, slowAns;
      std::vector<double> LAGAns;
      int status;
      int status1;  

      ieta=IEta[j]; iphi=IPhi[j]; depth=Depth[j];

      wTime=0;
      Double_t temp=0;

      for (int i=0; i<10; i++) {
	inputCaloSample.push_back(Charge[j][i]+Pedestal[j][i]);
	inputPedestal.push_back(Pedestal[j][i]);
	inputGain.push_back(Gain[j][i]);
      	//inputCaloSample.push_back(pulse[j][i]);
      	//inputPedestal.push_back(0);
      	//inputGain.push_back(1.0);
      	
      	Pulse[i]=Charge[j][i];
	Ped[i]=Pedestal[j][i];

	wTime+=i*Pulse[i];
	temp+=Pulse[i];

      }

      wTime=wTime/temp;

      psFitOOTpuCorr_->apply(inputCaloSample,inputPedestal,inputGain,offlineAns, status);

//begin miao's LAG fit
      double totalQ = 0.0;
      double tinyped = 0.0;
      for (int oo=0; oo<3; oo++) { tinyped+=Pulse[oo];}
      for (int oo=0; oo<10; oo++) { totalQ+=Pulse[oo];}

      tinyped/=3;
      totalQ-=10*tinyped;

      if (totalQ<20.0 || totalQ>600.0) continue;
      int ibin= floor( (totalQ-20.0)/10.0 );

      int shapenum;
      if ( abs(ieta)<17 ) {shapenum = 701 + ibin;}
      else {shapenum = 801 + ibin;}
      psFitOOTpuCorrLAG_->setPulseShapeTemplate(theHcalPulseShapes_.getShape(shapenum));

      psFitOOTpuCorrLAG_->apply(inputCaloSample,inputPedestal,inputGain,LAGAns, status1);


      m2Charge=offlineAns[0];
      m2Time=offlineAns[1];
      m2Ped=offlineAns[2];
      m2Chi=offlineAns[3];
      m2Status=status;

//begin miao's
      LAGCharge=LAGAns[0];
      LAGTime=LAGAns[1];
      LAGPed=LAGAns[2];
      LAGChi=LAGAns[3];
      LAGStatus=status1;
//end miao's


      for (uint i=0; i<offlineAns.size(); i++) {
      	if (i>3 && uint(i-4) < 10) {
      	  m2Pulse[uint(i-4)] = offlineAns[i];
      	}
      }

//begin miao's
      double totalEnergy=0.;
      for (uint i=0; i<LAGAns.size(); i++) {
        if (i>3 && uint(i-4) < 10) {
          LAGPulse[uint(i-4)] = LAGAns[i];
          totalEnergy+=LAGPulse[uint(i-4)];
        }
      }
//end miao's

      tout->Fill();

//miao's
      if(totalEnergy>20.0){
          chi2_M2->Fill(ibin*10.+20.,m2Chi);
          if(LAGChi < 1000) chi2_LAG->Fill(ibin*10.+20.,LAGChi);
      }
      if(totalEnergy>20.0 && m2Chi<=100 && LAGChi<=100){
          LAG_M2->Fill(m2Chi,LAGChi);
          LAG_M2_Ch->Fill(m2Chi,LAGChi,ibin*10.+20.);
      }


    }
  }

}

void Analysis::Finish()
{

  fout->cd();
  fout->Write();
  fout->Close();
}
