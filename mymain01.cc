// Main program to generate charged hadron spectra from photon-initiated
// hard processes, by combining sub-runs with direct or resolved photons
// or by generating all with contributions in a single run.

// In case of photon-photon interactions four different contributions are
// present:              ProcessType:
//  - resolved-resolved  1
//  - resolved-direct    2
//  - direct-resolved    3
//  - direct-direct      4
// Events can be generated either with photon beams or with photons emitted
// from lepton beams.

// In case of photon-proton interaction two contributions are present
//  - resolved
//  - direct
// When the photon is from beam A the relevant contributions are
// set with "Photon:ProcessType" values 1 (resolved) and 3 (direct) as the
// convention follows the photon-photon case above.
// Also lepton->photon + proton can be generated.

// Stdlib header file for input and output.
#include <iostream>

// ROOT, for histogramming.
#include "TH1.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TStyle.h"

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {

  // Generator.
  Pythia pythia;

  // Decrease the output.
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Next:numberCount = 0");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Shorthand for some public members of pythia (also static ones).
  Settings& settings  = pythia.settings;
  const Info& info = pythia.info;

  // Photon-proton collisions.
  bool photonProton         = true;

  // Generate photon-photon events in leptonic or photon beams.
  bool photonsFromElectrons = true;

  // Each contributions separately or in a one combined run.
  bool automaticMix         = false;

  // Optionally use different PDFs from LHAPDF for hard process.
  // Requires linkin with LHAPDF5.
  // pythia.readString("PDF:useHard = on");
  // pythia.readString("PDF:GammaHardSet = LHAPDF5:SASG.LHgrid/5");

  // Beam parameters.
  pythia.readString("Beams:eCM = 200.");

  // Set up beam particles for (electron -> photon) + proton.
  if ( photonProton) {
    if ( photonsFromElectrons) {
      pythia.readString("Beams:idA = 11");
      pythia.readString("Beams:idB = 2212");
      pythia.readString("PDF:beamA2gamma = on");

    // Set up beam particles for photon + proton.
    } else {
      pythia.readString("Beams:idA = 22");
      pythia.readString("Beams:idB = 2212");
    }

  // Set up beam particles for photon-photon in e+e-.
  } else if ( photonsFromElectrons) {
    pythia.readString("Beams:idA = -11");
    pythia.readString("Beams:idB =  11");
    pythia.readString("PDF:beamA2gamma = on");
    pythia.readString("PDF:beamB2gamma = on");

  // Set up beam particles for photon-photon.
  } else {
    pythia.readString("Beams:idA = 22");
    pythia.readString("Beams:idB = 22");
  }

  // Cuts on photon virtuality and invariant mass of gamma-gamma/hadron pair.
  if ( photonsFromElectrons) {
    pythia.readString("Photon:Q2max = 1.0");
    pythia.readString("Photon:Wmin  = 10.0");
  }

  // For photon-proton increase pT0Ref (for better agreement with HERA data).
  // Photon-photon has a new default pT0 parametrization tuned to LEP data.
  if ( photonProton)
    pythia.readString("MultipartonInteractions:pT0Ref = 3.00");

  // Limit partonic pThat.
  settings.parm("PhaseSpace:pTHatMin", 5.0);

  // Reset statistics after each subrun.
  pythia.readString("Stat:reset = on");

  // Parameters for histograms.
  double pTmin = 0.0;
  double pTmax = 5.0;
  int nBinsPT  = 100;
  
  TApplication theApp("hist", &argc, argv);
  TFile* outFile = new TFile("photoproduction.root", "RECREATE");
  //TTree* data = new TTree("data","Transverse Momentum Values of all Charged Final State Particles");
  
  //TH1F *mult = new TH1F("mult","charged multiplicity", 100, -0.5, 799.5);
  
  

  // Initialize the histograms.
  Hist pTtot("Total charged hadron pT distribution", nBinsPT, pTmin, pTmax);
  Hist pTresres("Resolved-resolved contribution for pT distribution", nBinsPT, pTmin, pTmax);
  Hist pTresdir("Resolved-direct contribution for pT distribution", nBinsPT, pTmin, pTmax);
  Hist pTdirres("Direct-resolved contribution for pT distribution", nBinsPT, pTmin, pTmax);
  Hist pTdirdir("Direct-direct contribution for pT distribution", nBinsPT, pTmin, pTmax);
  Hist pTiRun("Contribution from Run i for pT distribution", nBinsPT, pTmin, pTmax);
  
  //ROOT Histogram : Attach R at the end of the original histograms
  TH1F *pTtotR = new TH1F("pTtotR","Total charged hadron pT distribution;Transverse Momentum of all Charged Final State Particles ;Events ", nBinsPT, pTmin, pTmax);
  TH1F *pTresresR = new TH1F("pTresresR","Resolved-resolved contribution for pT distribution;Transverse Momentum Res-Res ;Events ", nBinsPT, pTmin, pTmax);
  TH1F *pTresdirR = new TH1F("pTresdirR","Resolved-direct contribution for pT distribution;Transverse Momentum Res-Dir ;Events", nBinsPT, pTmin, pTmax);
  TH1F *pTdirresR = new TH1F("pTdirresR","Direct-resolved contribution for pT distribution;Transverse Momentum Dir-Res ;Events", nBinsPT, pTmin, pTmax);
  TH1F *pTdirdirR = new TH1F("pTdirdirR","Direct-direct contribution for pT distribution;Transverse Momentum Dir-Dir ;Events", nBinsPT, pTmin, pTmax);
  TH1F *pTiRunR = new TH1F("pTiRunR","(Contribution from Run i for pT distribution)", nBinsPT, pTmin, pTmax);
  //double pmin = 0.0;
  //double pmax = 5.0;
  //int nBinsp  = 100;
  TH1F *momentumx = new TH1F("momentumx","Momentum of Final state Charged Particles in x Direction; Px; Events",nBinsPT, pTmin, pTmax);
  TH1F *momentumy = new TH1F("momentumy","Momentum of Final state Charged Particles in y Direction; Py; Events",nBinsPT, pTmin, pTmax);
  TH1F *momentumz = new TH1F("momentumz","Momentum of Final state Charged Particles in z Direction; Pz; Events",nBinsPT, pTmin, pTmax);
  TH1F *energy = new TH1F("energy","Energy of Final state Charged Particles; Energy; Events",100,0,100);
  TH1F *mass = new TH1F("mass","Mass of Final state Charged Particles",100,0,2);

  // Initialize hard QCD processes with 0, 1, or 2 initial photons.
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhotonParton:all = on");
  //pythia.readString("PhotonParton:ggm2qqbar = on");
  //pythia.readString("PhotonParton:ggm2ccbar = on");
  //pythia.readString("PhotonParton:ggm2bbbar = on");
  //pythia.readString("PhotonParton:qgm2qg = on");
  //pythia.readString("PhotonParton:qgm2qgm = on");
  
  if ( !photonProton) {
    pythia.readString("PhotonCollision:gmgm2qqbar = on");
    pythia.readString("PhotonCollision:gmgm2ccbar = on");
    pythia.readString("PhotonCollision:gmgm2bbbar = on");
  }

  // Number of runs.
  int nRuns = photonProton ? 2 : 4;
  if (automaticMix) nRuns = 1;

  // Number of events per run.
  int nEvent = 10000;
  
  //counters for the Beam Particles
  int e_count_scat = 0;
  int p_count_scat = 0;
  int e_count = 0;	//Total Beam Particles from electron beam
  int p_count = 0;	//Total Beam Particles from proton beam
  int beamparticles_c = 0; //Total Beam Particles

  // Loop over relevant processes.
  for ( int iRun = 1; iRun < nRuns + 1; ++iRun) {

    // Turn of MPIs for processes with unresolved photons.
    if (iRun == 2) pythia.readString("PartonLevel:MPI = off");

    // For photon+proton direct contribution with processType = 3.
    if (photonProton && iRun == 2) iRun = 3;

    // Set the type of gamma-gamma process:
    // 0 = mix of all below,
    // 1 = resolved-resolved,
    // 2 = resolved-direct,
    // 3 = direct-resolved,
    // 4 = direct-direct.
    if (automaticMix) settings.mode("Photon:ProcessType", 0);
    else              settings.mode("Photon:ProcessType", iRun);

    // Initialize the generator.
    pythia.init();

    // Clear the histogram.
    pTiRun.null();
    pTiRunR->Reset();

    // Begin event loop. Skip if fails.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate next event.
      if (!pythia.next()) continue;

      // List the first process and event for each run.
      if (iEvent == 0) {
        pythia.process.list();
        pythia.event.list();
      }

      // Possible event weights.
      double weight = info.weight();	//used to Normalize Histograms

      // Loop over event record and find charged final state particles.
      for (int i = 0; i < pythia.event.size(); ++i){
        if ( pythia.event[i].isFinal() && pythia.event[i].isCharged() ) {
			
          // Store the pT value.
          double pTch = pythia.event[i].pT();
          //Store 4 momentum of the final state charged particles in the root file
          double px = pythia.event[i].px(); momentumx->Fill(px);
          double py = pythia.event[i].py(); momentumy->Fill(py);
          double pz = pythia.event[i].py(); momentumz->Fill(pz);
          double e = pythia.event[i].e(); energy->Fill(e);
          double m = pythia.event[i].m(); mass->Fill(m);
          
          if (pythia.event[i].status() >= 11 && pythia.event[i].status() <= 19){
          	beamparticles_c ++;
          	}
          	
          //if (pythia.event[i].status() >= 21 && pythia.event[i].status() <= 29)
          //if (pythia.event[i].status() >= 31 && pythia.event[i].status() <= 39)
          
          /*
          //to check for the scattered electrons :
          if (pythia.event[i].id() == 11)
          {
          	if (pythia.event[i].status() == 14 || pythia.event[i].status() == 15)
          	e_count_scat++;
          	
          	if (pythia.event[i].status() >= 11 && pythia.event[i].status() <= 13)
          	e_count++;
          }
          
          if (pythia.event[i].id() == 2212)
          {
          	if (pythia.event[i].status() == 14 || pythia.event[i].status() == 15)
          	p_count_scat++;
          	
          	if (pythia.event[i].status() >= 11 && pythia.event[i].status() <= 13)
          	p_count++;
          }
          */
          
          

          // Fill the correct histogram depending on the process type.
          
          if (automaticMix) {
            pTtot.fill(pTch, weight);
            pTtotR->Fill(pTch, weight);
            
            if (info.photonMode() == 1) pTresres.fill(pTch, weight);
            if (info.photonMode() == 2) pTresdir.fill(pTch, weight);
            if (info.photonMode() == 3) pTdirres.fill(pTch, weight);
            if (info.photonMode() == 4) pTdirdir.fill(pTch, weight);
            
            if (info.photonMode() == 1) pTresresR->Fill(pTch, weight);
            if (info.photonMode() == 2) pTresdirR->Fill(pTch, weight);
            if (info.photonMode() == 3) pTdirresR->Fill(pTch, weight);
            if (info.photonMode() == 4) pTdirdirR->Fill(pTch, weight);
          } else {
            pTiRun.fill(pTch, weight);
            pTiRunR->Fill(pTch, weight);
          }
        }
      }
    } // End of event loop.

    // Show statistics after each run.
    pythia.stat();

    // Normalize to cross section [mb].
    double sigmaNorm = info.sigmaGen() / info.weightSum();
    double pTBin     = (pTmax - pTmin) / (1. * nBinsPT);
    double val   = sigmaNorm / pTBin;

    // For mix of all contributions normalize with total cross section.
    if (automaticMix) {
      pTtot    *= val;
      pTtotR->SetBinContent(1,val);
      momentumx->SetBinContent(1,val);
      momentumy->SetBinContent(1,val);
      momentumz->SetBinContent(1,val);
      mass->SetBinContent(1,(sigmaNorm/0.05));
      pTresres *= val;	pTresresR->SetBinContent(1,val);
      pTresdir *= val;	pTresdirR->SetBinContent(1,val);
      pTdirres *= val;	pTdirresR->SetBinContent(1,val);
      pTdirdir *= val;	pTdirdirR->SetBinContent(1,val);

    // For each contribution normalize with cross section for the given run.
    } else {
      pTiRun *= val;
      pTiRunR->SetBinContent(1,val);
      momentumx->SetBinContent(1,val);
      momentumy->SetBinContent(1,val);
      momentumz->SetBinContent(1,val);
      
      if (iRun == 1) { pTresres = pTiRun; pTresresR->Add(pTiRunR); }
      if (iRun == 2) { pTresdir = pTiRun; pTresdirR->Add(pTiRunR); }
      if (iRun == 3) { pTdirres = pTiRun; pTdirresR->Add(pTiRunR); }
      if (iRun == 4) { pTdirdir = pTiRun; pTdirdirR->Add(pTiRunR); }
      
      pTtot += pTiRun;
      pTtotR->Add(pTiRunR);
    }

  // End of loop over runs.
  }

  // Print histograms.
  cout << pTresres << pTresdir << pTdirres << pTdirdir << pTtot << endl;
  cout << "------------------------------------------------------------------------------------------------" << endl;
  cout << "Number of Events : " << nEvent << endl;
  cout << "Number of Charged Final State Beam Particles : " << beamparticles_c << endl;
  /*
  cout << "Number of Outgoing Scattered Electrons (DIS or Elasticaly) : " << e_count_scat << endl;
  cout << "Number of Beam Particles (Electrons) : " << e_count << endl;
  cout << "Number of Outgoing Scattered Protons (DIS or Elasticaly) : " << p_count_scat << endl;
  cout << "Number of Beam Particles (Protons) : " << p_count << endl;
  */
  
  
  pTresresR->SetLineWidth(3);	pTresresR->Write(); 
  pTresdirR->SetLineWidth(3);	pTresdirR->Write(); 
  pTdirresR->SetLineWidth(3);	pTdirresR->Write();
  pTdirdirR->SetLineWidth(3);	pTdirdirR->Write();
  //pTiRunR->SetLineWidth(3);		pTiRunR->Write();
  pTtotR->SetLineWidth(3);		pTtotR->Write();
  momentumx->SetLineWidth(3);	momentumx->Write();
  momentumy->SetLineWidth(3);	momentumy->Write();
  momentumz->SetLineWidth(3);	momentumz->Write();
  energy->SetLineWidth(3);		energy->Write();
  mass->SetLineWidth(3);		mass->Write();
  
  gStyle->SetLineWidth(2);
  
  delete outFile;

  // Done.
  return 0;
}
