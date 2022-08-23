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
  // pythia.readString("Init:showChangedSettings = off");
  // pythia.readString("Init:showChangedParticleData = off");
  // pythia.readString("Next:numberCount = 0");
  // pythia.readString("Next:numberShowInfo = 0");
  // pythia.readString("Next:numberShowProcess = 0");
  // pythia.readString("Next:numberShowEvent = 0");

  // Shorthand for some public members of pythia (also static ones).
  Settings& settings  = pythia.settings;
  const Info& info = pythia.info;

  // Photon-proton collisions.
  bool photonProton = true;

  // Generate photon-photon events in leptonic or photon beams.
  bool photonsFromElectrons = true;

  // Each contributions separately or in a one combined run.
  bool automaticMix = false;

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
  int nEvent = 100000;
  
  //counters to check the Beam Particles for res-res and dir-res contributions
  int count_rr[5];// = new int[5];
  int count_dr[5];// = new int[5];
  //{0,0,0,0,0,0,0,0,0,0,0,0};
  int d, d1, d2;
  int dr = 0, rr = 0;

  // Loop over relevant processes. (In photonproton = true case iRun will only have values 1 and 3)
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
      		
      		//********************************************************************************************************************
      		//TO CHECK FOR ELECTRON BEAM REMNANTS rr= resolved-resolved, dr=direct-resolved
      		if (iRun == 1 && pythia.event[i].id() == 11 && pythia.event[i].status() == 63 && pythia.event[i].isFinal()) rr++;
      		if (iRun == 3 && pythia.event[i].id() == 11 && pythia.event[i].status() == 63 && pythia.event[i].isFinal()) dr++;
      		
			//Method 2: By checking the Mother and the daughter particles in each run.
			if (iRun == 1)	//Checking for Resolved-Resolved electrons
			{
         		if(pythia.event[i].id() == 11 && pythia.event[i].mother1() == 0 && pythia.event[i].mother2() == 0) {
          			count_rr[0]++; //Stores the Beam electrons.
          	   	d=i;
		      	d1 = pythia.event[d].daughter1();
		      	d2 = pythia.event[d].daughter2();
		      	
		      	if(pythia.event[d1].id() == 11 && pythia.event[i].isFinal())
		      		count_rr[1]++;	//daughter1 of the beam electron
		      	if(pythia.event[d2].id() == 11 && pythia.event[i].isFinal())
		      		count_rr[2]++;	//daughter2 of the beam electron
		      	if(pythia.event[i].isFinal())
		      		count_rr[3]++;
		     	}
		     }
		     if(iRun == 3)	//Checking for Direct-Resolved electrons
		     {
		     	if(pythia.event[i].id() == 11 && pythia.event[i].mother1() == 0 && pythia.event[i].mother2() == 0) {
          			count_dr[0]++; //Stores the Beam electrons.
          	   	d=i;
		      	d1 = pythia.event[d].daughter1();
		      	d2 = pythia.event[d].daughter2();
		      	
		      	if(pythia.event[d1].id() == 11 && pythia.event[i].isFinal())
		      		count_dr[1]++;	//daughter1 of the beam electron
		      	if(pythia.event[d2].id() == 11 && pythia.event[i].isFinal())
		      		count_dr[2]++;	//daughter2 of the beam electron
		      	if(pythia.event[i].isFinal())
		      		count_dr[3]++;
		     	}
		     }
		     //********************************************************************************************************************
		     
		if ( pythia.event[i].isFinal() && pythia.event[i].isCharged() ) {
			// Store the pT value.
	        double pTch = pythia.event[i].pT();
	        //Store 4 momentum of the final state charged particles in the root file
	        double px = pythia.event[i].px(); momentumx->Fill(px);
	        double py = pythia.event[i].py(); momentumy->Fill(py);
	        double pz = pythia.event[i].pz(); momentumz->Fill(pz);
	        double e = pythia.event[i].e(); energy->Fill(e);
	        double m = pythia.event[i].m(); mass->Fill(m);
					
		// Fill the correct histogram depending on the process type.
        if (automaticMix) {
			pTtotR->Fill(pTch, weight);
            if (info.photonMode() == 1) pTresresR->Fill(pTch, weight);
            if (info.photonMode() == 2) pTresdirR->Fill(pTch, weight);
            if (info.photonMode() == 3) pTdirresR->Fill(pTch, weight);
            if (info.photonMode() == 4) pTdirdirR->Fill(pTch, weight);
          } else {
            pTiRunR->Fill(pTch, weight);
          }
        }
      }
    } // End of event loop for i Run

    // Show statistics after each run.
    pythia.stat();

    // Normalize to cross section [mb].
    double sigmaNorm = info.sigmaGen() / info.weightSum();
    double pTBin     = (pTmax - pTmin) / (1. * nBinsPT);
    double val   = sigmaNorm / pTBin;

    // For mix of all contributions normalize with total cross section.
    if (automaticMix) {
      pTtotR->SetBinContent(1,val);
      momentumx->SetBinContent(1,val);
      momentumy->SetBinContent(1,val);
      momentumz->SetBinContent(1,val);

    // For each contribution normalize with cross section for the given run.
    } else {
      pTiRunR->SetBinContent(1,val);
      momentumx->SetBinContent(1,val);
      momentumy->SetBinContent(1,val);
      momentumz->SetBinContent(1,val);
      
      if (iRun == 1) pTresresR->Add(pTiRunR); 
      if (iRun == 2) pTresdirR->Add(pTiRunR); 
      if (iRun == 3) pTdirresR->Add(pTiRunR); 
      if (iRun == 4) pTdirdirR->Add(pTiRunR); 
      
      pTtotR->Add(pTiRunR);
    }

  // End of loop over runs.
  }
	cout << "" << endl << "Number of Events : " << nEvent << endl;
	//PRINT OUT THE COUNTERS:
	cout << "----------------------------------------------COUNTERS---------------------------------------------" << endl;
	cout << "Method 1 : Checking Status (=63) directly :" << endl;
	cout << "resolved-resolved final state electrons with status codes as 63 :" << rr << endl; 
	cout << "direct-resolved final state electrons with status codes as 63 :" << dr << endl;
	cout << "---------------------------------------------------------------------------------------------------" << endl;
	cout << "Method 2 : Checking Mother-Daughter Status :" << endl;
	cout << "For Res-Res :" << endl;
	cout << "		" << "Beam Electrons (No Mother Particles, m1=0 and m2=0):" << count_rr[0] << endl;
	cout << "		" << "Daughter, d1, of the incoming beam electron is a final state electron :" << count_rr[1] << endl;
	cout << "		" << "Daughter, d2, of the incoming beam electron is a final state electron :" << count_rr[2] << endl;
	cout << "		" << "Final State Beam Electron (Should be zero!)" << count_rr[3] << endl;
	cout << "For Dir-Res :" << endl;
	cout << "		" << "Beam Electrons (No Mother Particles, m1=0 and m2=0):" << count_dr[0] << endl;
	cout << "		" << "Daughter, d1, of the incoming beam electron is a final state electron :" << count_dr[1] << endl;
	cout << "		" << "Daughter, d2, of the incoming beam electron is a final state electron :" << count_dr[2] << endl;
	cout << "		" << "Final State Beam Electron (Should be zero!)" << count_dr[3] << endl;
	
  
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
