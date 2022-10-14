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


#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"
#include "Pythia8/Pythia.h"
#include <iostream>

using namespace Pythia8;

int main() {
	
	// Number of events per run.
  int nEvent = 1000;
  
  // Generator.
  Pythia pythia;

  // Decrease the output.
  pythia.readString("Init:showMultipartonInteractions = on");
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
  bool automaticMix = true;

  // Optionally use different PDFs from LHAPDF for hard process.
  // Requires linking with LHAPDF5.
  pythia.readString("PDF:useHard = on");
  pythia.readString("HardQCD:all = on");
  // pythia.readString("PDF:GammaHardSet = LHAPDF5:SASG.LHgrid/5");

  // Beam parameters.
  pythia.readString("Beams:eCM = 140.");	//CoM Energy ranges from 134–277 GeV at HERA 
  pythia.readString("Beams:eA  = 27.5.");
  pythia.readString("Beams:eB  = 820.");

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
    pythia.readString("Photon:Wmin  = 134.0");
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
  
  //TApplication theApp("hist", &argc, argv);
  //TFile *outFile = new TFile("photoproduction.root","RECREATE");
  
  TFile *data = new TFile("data.root","RECREATE");
	//Create a tree for Combined Events :
  TTree *TC = new TTree("combinedevents", "4-Momentum of all Charged Final State Particles (Combined)");
  std::vector<Float_t> px,py,pz,e,bpx,bpy,bpz,be,gmpx,gmpy,gmpz,gme;
  std::vector<Int_t> dirflag, resflag, events;
  TC->Branch("events", "std::vector<Int_t>",  &events);
  TC->Branch("px", "std::vector<Float_t>",  &px);
  TC->Branch("py", "std::vector<Float_t>",  &py);
  TC->Branch("pz", "std::vector<Float_t>",  &pz);
  TC->Branch("e", "std::vector<Float_t>",  &e);
  //4-Momentum of all Beam electrons
  TC->Branch("bpx", "std::vector<Float_t>",  &bpx);
  TC->Branch("bpy", "std::vector<Float_t>",  &bpy);
  TC->Branch("bpz", "std::vector<Float_t>",  &bpz);
  TC->Branch("be", "std::vector<Float_t>",  &be);
  //4-Momentum of all Beam Photoproduction-Photons 
  TC->Branch("gmpx", "std::vector<Float_t>",  &gmpx);
  TC->Branch("gmpy", "std::vector<Float_t>",  &gmpy);
  TC->Branch("gmpz", "std::vector<Float_t>",  &gmpz);
  TC->Branch("gme", "std::vector<Float_t>",  &gme);
  //Set flags for direct or resolved events
  TC->Branch("dirflag", "std::vector<Int_t>",  &dirflag);
  TC->Branch("resflag", "std::vector<Int_t>",  &resflag);
	
	//Tree to save the quark inititated process, for quark jets
	TTree *qqbar2qqbar = new TTree("qqbar2qqbar", "4-Momentum of all qqbar2qqbar process");
  std::vector<Float_t> qqpx,qqpy,qqpz,qqe;
  qqbar2qqbar->Branch("px", "std::vector<Float_t>",  &qqpx);
  qqbar2qqbar->Branch("py", "std::vector<Float_t>",  &qqpy);
  qqbar2qqbar->Branch("pz", "std::vector<Float_t>",  &qqpz);
  qqbar2qqbar->Branch("e", "std::vector<Float_t>",  &qqe);
  
  //Tree to save the gluon inititated process, for gluon jets
	TTree *gg2gg = new TTree("gg2gg", "4-Momentum of all gg2gg process");
  std::vector<Float_t> ggpx,ggpy,ggpz,gge;
  gg2gg->Branch("px", "std::vector<Float_t>",  &ggpx);
  gg2gg->Branch("py", "std::vector<Float_t>",  &ggpy);
  gg2gg->Branch("pz", "std::vector<Float_t>",  &ggpz);
  gg2gg->Branch("e", "std::vector<Float_t>",  &gge);
  
	//ROOT Histogram : Attach R at the end of the original histograms : COMBINED EVENTS HISTOGRAMS.
  TH1F *pTtotR = new TH1F("pTtotR","Total charged hadron pT distribution; pT of Charged Final State Particles ;Events ", nBinsPT, pTmin, pTmax);
  TH1F *pTresresR = new TH1F("pTresresR","Resolved-resolved contribution for pT distribution; pT  Res-Res ;Events ", nBinsPT, pTmin, pTmax);
  TH1F *pTresdirR = new TH1F("pTresdirR","Resolved-direct contribution for pT distribution; pT Res-Dir ;Events", nBinsPT, pTmin, pTmax);
  TH1F *pTdirresR = new TH1F("pTdirresR","Direct-resolved contribution for pT distribution; pT  Dir-Res ;Events", nBinsPT, pTmin, pTmax);
  TH1F *pTdirdirR = new TH1F("pTdirdirR","Direct-direct contribution for pT distribution; pT  Dir-Dir ;Events", nBinsPT, pTmin, pTmax);
  TH1F *pTiRunR = new TH1F("pTiRunR","(Contribution from Run i for pT distribution)", nBinsPT, pTmin, pTmax);
  

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
  
  //nRuns = 1 for just one one run with Mixing turned on

  int tsubp = 0, rescount = 0, dircount = 0, beame = 0, qq=0, gg=0, qsub=0, gsub=0, gm=0, c=0, chk=0;
  int dsub=0,rsub=0,csub=0;
  double ben, gen, y;

  // Loop over relevant processes. (In photonproton = true case iRun will only have values 1 and 3)
  for ( int iRun = 1; iRun < nRuns + 1; iRun++) {
	
    // Turn off MPIs for processes with unresolved photons.
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
			
			events.clear();
			
			//Reinitializing variables for the next event
			px.clear(); py.clear(); pz.clear(); e.clear();	//Combined Process' 4 momentum
			bpx.clear(); bpy.clear(); bpz.clear(); be.clear();	//Beam Remnant Electrons 4 momentum
			gmpx.clear(); gmpy.clear(); gmpz.clear(); gme.clear();	//Total Photons from beam electrons
			dirflag.clear(); resflag.clear();	//Int Flags to check for direct/resolved events
			qqpx.clear(); qqpy.clear(); qqpz.clear(); qqe.clear();	//qqbar2qqbar subprocess
			ggpx.clear(); ggpy.clear(); ggpz.clear(); gge.clear();	//gg2gg subprocess
			bool dircheck = false, rescheck = false, qqsub = false, ggsub = false;	//Flags

			events.push_back(iEvent);	//Pushing Back the event number
		
      // Generate next event.
      if (!pythia.next()) continue;

      // List the first process and event for each run.
        if (iEvent == 0) {
        pythia.process.list();
        pythia.event.list();
      }
	
      // Possible event weights.
      double weight = info.weight();	//used to Normalize Histograms

		  if (pythia.info.code()!=0){
		  	if (281 <= pythia.info.code() && pythia.info.code() <= 284) {dircheck=true; dircount++;}	
		  	else {rescheck = true; rescount++; }
		  		
		  	if (pythia.info.code() == 114) {qqsub = true; qq++;}
		  	if (pythia.info.code() == 111) {ggsub = true; gg++;}
		  }
		  //Filling the int-flags with code() if true, 0 is false
		  if(dircheck) {dirflag.push_back(pythia.info.code());}
	        	else dirflag.push_back(0);
      if(rescheck) {resflag.push_back(pythia.info.code());}
      	else resflag.push_back(0);
    	
      // Loop over event record and find charged final state particles.
      for (int i = 0; i < pythia.event.size(); ++i){
      	
      	//Condition for storing the 4 momentum of beam particles & skip the code:
      	if (pythia.event[i].status() == -13 || pythia.event[i].status() == -12)
      	if ( pythia.event[i].id() == 22 || pythia.event[i].id() == 11){
      		//store the photon from the beam electron
      		if(pythia.event[i].id() == 22){
				   	gmpx.push_back(pythia.event[i].px());
					  gmpy.push_back(pythia.event[i].py());
					  gmpz.push_back(pythia.event[i].pz());
					  gme.push_back(pythia.event[i].e());
					  gen = pythia.event[i].e();
					  //cout << "incoming photon : " << gen << endl;
					  gm++;
				  }
				  if(pythia.event[i].id() == 11){
				  	bpx.push_back(pythia.event[i].px());
					  bpy.push_back(pythia.event[i].py());
					  bpz.push_back(pythia.event[i].pz());
					  be.push_back(pythia.event[i].e());
					  ben = pythia.event[i].e();
					  //cout << "incoming electron : " <<  ben << endl;
					  beame++;
				  }
				  chk++; continue;
				}		    
		    
		    //Storing 4mom of Final state charged particles
		    if ( pythia.event[i].isFinal() && pythia.event[i].isCharged() ) {
					
					//Quark Subprocess
				  if (pythia.info.code() == 114)
				  {
				  	qqpx.push_back(pythia.event[i].px());
					  qqpy.push_back(pythia.event[i].py());
					  qqpz.push_back(pythia.event[i].pz());
					  qqe.push_back(pythia.event[i].e());
					  qsub++;
				  }
				  //Gluon Subprocess
				  if (pythia.info.code() == 111)
				  {
				  	ggpx.push_back(pythia.event[i].px());
					  ggpy.push_back(pythia.event[i].py());
					  ggpz.push_back(pythia.event[i].pz());
					  gge.push_back(pythia.event[i].e());
					  gsub++;
				  }
		    	        
	        //Store 4-momentum of the final state charged particles (Combined Events)
	        px.push_back(pythia.event[i].px());
	        py.push_back(pythia.event[i].py());
	        pz.push_back(pythia.event[i].pz());
	        e.push_back(pythia.event[i].e());
	        
	        //Checking the number of constituents per event cs=dsub+rsub
	        csub++;
	        if(dircheck) dsub++;
	        if(rescheck) rsub++;
	        					
	        // Store the pT value of the event
	        double pTch = pythia.event[i].pT();
	        
					// Fill the correct histogram depending on the process type.
		      if (automaticMix) {
						pTtotR->Fill(pTch, weight);
	          if (info.photonMode() == 1) pTresresR->Fill(pTch, weight);
	          if (info.photonMode() == 2) pTresdirR->Fill(pTch, weight);
	          if (info.photonMode() == 3) pTdirresR->Fill(pTch, weight);
	          if (info.photonMode() == 4) pTdirdirR->Fill(pTch, weight);
		      } 
		      else pTiRunR->Fill(pTch, weight);
        }	//end-of if condition to check for final state hadrons
      }//end-of for (event size) loop
      TC->Fill(); qqbar2qqbar->Fill(); gg2gg->Fill(); 
    } // End of event loop for i Run; only 1 run in automatic mix

    // Show statistics after each run.
    pythia.stat();

    // Normalize to cross section [mb].
    double sigmaNorm = info.sigmaGen() / info.weightSum();
    double pTBin     = (pTmax - pTmin) / (1. * nBinsPT);
    double val   = sigmaNorm / pTBin;

    // For mix of all contributions normalize with total cross section.
    if (automaticMix) {
      pTtotR->SetBinContent(1,val);

    // For each contribution normalize with cross section for the given run.
    } else {
      pTiRunR->SetBinContent(1,val);
      if (iRun == 1) pTresresR->Add(pTiRunR); 
      if (iRun == 2) pTresdirR->Add(pTiRunR); 
      if (iRun == 3) pTdirresR->Add(pTiRunR); 
      if (iRun == 4) pTdirdirR->Add(pTiRunR); 
      pTtotR->Add(pTiRunR);
    }
  }// End of for loop over runs.
	
  //cout << "Fraction of Resolved and Direct Subprocesses (based on their Subprocess Codes) for "<< nEvent << " Events :" << endl;
  cout << "\nTotal Events : " << nEvent << endl;
  cout << "Total beam remnant electrons (should be same as Events) :" << beame << endl;
  cout << "Total Photons (incoming beam-inside-beam) from beam electrons : " << gm << endl;
  
  cout << "\nTotal Direct Events :" << dircount << endl;
  cout << "Total Direct Constituents :" << dsub << endl;
  cout << "Total Resolved Events :" << rescount << endl;
  cout << "Total Resolved Constituents :" << rsub << endl;
  cout << "Total Combined Events :" << nEvent << endl;
  cout << "Total Combined Constituents :" << csub << endl;
  
  cout << "\nTotal quark initiated events (Code 114):" << qq << endl;
  cout << "Total quark initiated charged final state particles (qqbar2qqbar) :" << qsub << endl;
  cout << "Total gluon initiated events (Code 111):" << gg << endl;
  cout << "Total gluon initiated charged final state particles (gg2gg) :" << gsub << endl;
  cout << c << endl;
    
  gStyle->SetLineWidth(2);
  
  data->Write();
  data->Close();
  delete data;
}
