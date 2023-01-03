//	Event generation for PHOTPRODUCTION process
//	cd /home/siddharth/HEP_Projects/Project1_PhotoProduction/EventGen\ _Photoproduction/
//
// Authors: Ilkka Helenius <ilkka.m.helenius@jyu.fi>
// Siddharth Singh
//
// Keywords: photon beam; photoproduction; photon-photon;
// In case of photon-photon interactions four different contributions are
// present:              ProcessType:
//  - resolved-resolved  1
//  - resolved-direct    2
//  - direct-resolved    3
//  - direct-direct      4
// 
// PYTHIA examples: main70.cc and main69.cc
//***************************************************************************************

#include "TH1.h"
#include "TH3.h"
#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"
#include "TGraph.h"
#include "Pythia8/Pythia.h"
#include <iostream>
using namespace Pythia8;

int main() {
		
	// Number of events per run.
	int nEvent = 100000;

	// Generator.
	Pythia pythia;

	// Decrease the output.
	pythia.readString("Init:showMultipartonInteractions = on");
	pythia.readString("Init:showChangedSettings = on");
	pythia.readString("Init:showChangedParticleData = off");
	pythia.readString("Next:numberCount = 0");
	pythia.readString("Next:numberShowInfo = 0");
	pythia.readString("Next:numberShowProcess = 0");
	pythia.readString("Next:numberShowEvent = 0");

	// Shorthand for some public members of pythia (also static ones).
	Settings& settings  = pythia.settings;
	const Info& info = pythia.info;

	// Beam parameters : 63.2(10x100) 114.9(10x275) 141(18x275)(20x250)
	pythia.readString("Beams:frameType = 2");		// Beams of unequal energies
	pythia.readString("Beams:idA = -11");				// Beam A : Electrons/Positrons
	pythia.readString("Beams:idB = 2212");			// Beam B : Protons
	pythia.readString("Beams:eA = 18");				// electron : HERA=27.5GeV && EIC=18 GeV
	pythia.readString("Beams:eB = 275");				// proton : HERA=820|920Gev && EIC=275 GeV
	pythia.readString("PDF:beamA2gamma = on");	// Set PDF for beam photons

	// pythia.readString("Photon:Wmin  = 134.0");			// invariant mass_min of gamma-hadron pair
	// pythia.readString("Photon:Wmax  = 277.0");			// invariant mass_max of gamma-hadron pair

	// Switch relevant processes on :
	pythia.readString("HardQCD:all = on");									//For All Resolved Process
	// pythia.readString("HardQCD:gg2gg = on");									//gg2gg - Gluon Induced Events (code 111)
	// pythia.readString("HardQCD:gg2qqbar = on");								//gg2qqbar - (code 112)
	// pythia.readString("HardQCD:qg2qg = on");								//qg2qg - (code 113)
	// pythia.readString("HardQCD:qq2qq = on");									//qq2qq - Quark Induced Events (code 114)
	// pythia.readString("HardQCD:qqbar2gg = on");								//qqbar2gg - (code 115)
	// pythia.readString("HardQCD:qqbar2qqbarNew = on");				//qqbar2qqbarNew (code 116)

	// For photon-parton interaction :
	pythia.readString("PhotonParton:all = on");							// For All Direct Process
	// pythia.readString("PhotonParton:ggm2qqbar = on");						// Scattering g gamma → q qbar, (q = u,d,s) Code 271 (281).
	// pythia.readString("PhotonParton:ggm2ccbar = on");						// Scattering g gamma → c cbar. Code 272 (282).
	// pythia.readString("PhotonParton:ggm2bbbar = on");						// Scattering g gamma → b bbar. Code 273 (283).
	// pythia.readString("PhotonParton:qgm2qg = on");						// Scattering q gamma → q g. Code 274 (284).
	// pythia.readString("PhotonParton:qgm2qgm = on");							// Scattering q gamma → q gamma. Code 275 (285).

	// Set Event Settings : 
	settings.mode("Photon:ProcessType", 0);												// Switch on : Automatic Mix
	pythia.readString("Photon:Q2max = 1.0");											// Maximal Q2
	pythia.readString("MultipartonInteractions:pT0Ref = 3");			// Use tuned pT0ref for photon-hadron (for LEP2, EIC = 3)
	pythia.readString("PhaseSpace:pTHatMin = 5");								// Limit partonic pThat.
	pythia.readString("SpaceShower:pTmaxMatch = 1");							// Allow emissions up to the kinematical limit
	// pythia.readString("SpaceShower:dipoleRecoil = on");						// Set dipole recoil on. Necessary for DIS + shower.
	pythia.readString("PartonLevel:MPI = on");										// Master switch for multiparton interactions (default = on)
	pythia.readString("PartonLevel:all= on");											// Master switch for parton-event interactions
	pythia.settings.forceParm("PhaseSpace:pTHatMinDiverge",1.0);	// Extra pT cut to avoid the divergences in the limit pT→0
	// pythia.settings.parm("PhaseSpace:mHatMin",1.0);								// The minimum invariant mass.

	//ROOT for Storing Generated Events :
	TFile *data = new TFile("data.root","recreate");
	std::vector<Float_t> px,py,pz,e,eT,bpx,bpy,bpz,be,beT,inelasticity;			//All combined-events 
	std::vector<Int_t> dirflag, resflag, ggflag, qqflag, eventsid, scattering;					//Flags
	std::vector<Float_t> qpx,qpy,qpz,qe,qeT;																//quark initiated particles (qq)
	std::vector<Float_t> gpx,gpy,gpz,ge,geT;																//gluon initiated particles (gg)
	std::vector<Float_t> qgpx,qgpy,qgpz,qge,qgeT;														//gq

	//Create a tree for all Events :
	TTree *TC = new TTree("combinedevents", "4-Momentum of all Charged Final State Particles (Combined)");
	TC->Branch("eventid", "std::vector<Int_t>",  &eventsid);
	//scattering type: qq, gg or qg
	TC->Branch("scattering", "std::vector<Int_t>",  &scattering);
	TC->Branch("px", "std::vector<Float_t>",  &px);
	TC->Branch("py", "std::vector<Float_t>",  &py);
	TC->Branch("pz", "std::vector<Float_t>",  &pz);
	TC->Branch("e", "std::vector<Float_t>",  &e);
	TC->Branch("eT", "std::vector<Float_t>",  &eT);
	//4-Momentum of all Beam electrons :
	TC->Branch("bpx", "std::vector<Float_t>",  &bpx);
	TC->Branch("bpy", "std::vector<Float_t>",  &bpy);
	TC->Branch("bpz", "std::vector<Float_t>",  &bpz);
	TC->Branch("be", "std::vector<Float_t>",  &be);
	TC->Branch("beT", "std::vector<Float_t>",  &beT);
	TC->Branch("inelasticity", "std::vector<Float_t>",  &inelasticity);
	//Set flags for direct/resolved/gluon/quark processes :
	TC->Branch("dirflag", "std::vector<Int_t>",  &dirflag);
	TC->Branch("resflag", "std::vector<Int_t>",  &resflag);
	TC->Branch("ggflag", "std::vector<Int_t>",  &ggflag);
	TC->Branch("qqflag", "std::vector<Int_t>",  &qqflag);

	//4-Momentum of all quark-initiated final state charged particles (qq):
	TC->Branch("qpx", "std::vector<Float_t>",  &qpx);
	TC->Branch("qpy", "std::vector<Float_t>",  &qpy);
	TC->Branch("qpz", "std::vector<Float_t>",  &qpz);
	TC->Branch("qe", "std::vector<Float_t>",  &qe);
	TC->Branch("qeT", "std::vector<Float_t>",  &qeT);
	//4-Momentum of all gluon-initiated final state charged particles (gg):
	TC->Branch("gpx", "std::vector<Float_t>",  &gpx);
	TC->Branch("gpy", "std::vector<Float_t>",  &gpy);
	TC->Branch("gpz", "std::vector<Float_t>",  &gpz);
	TC->Branch("ge", "std::vector<Float_t>",  &ge);
	TC->Branch("geT", "std::vector<Float_t>", &geT);
	//4-Momentum of all gluon-initiated final state charged particles (gg):
	TC->Branch("qgpx", "std::vector<Float_t>",  &qgpx);
	TC->Branch("qgpy", "std::vector<Float_t>",  &qgpy);
	TC->Branch("qgpz", "std::vector<Float_t>",  &qgpz);
	TC->Branch("qge", "std::vector<Float_t>",  &qge);
	TC->Branch("qgeT", "std::vector<Float_t>", &qgeT);
	

	TGraph *gr_qq= new TGraph(); // subprocesses contribution
	double Xt=0;
	bool qq = false, qg = false, gg = false;	// gg code: 111, 115; qg code: 113,  284; qq code: 112, 281
	int sub[3]={0,0,0};
	// Parameters for histograms.
	double pTmin = 0.0;
	double pTmax = 5.0;
	int nBinsPT  = 100;
	//ROOT Histogram pT of Direct and Resolved Processes :
	TH1F *pTres = new TH1F("pTresresR","Resolved-resolved contribution for pT distribution; pT  Res-Res ;Events ", nBinsPT, pTmin, pTmax);
	TH1F *pTdir = new TH1F("pTdirresR","Direct-resolved contribution for pT distribution; pT  Dir-Res ;Events", nBinsPT, pTmin, pTmax);

	//Counters :
	int ev=0, rescount = 0, dircount = 0, beame = 0, c=0, chk=0, qcc=0, gcc=0, qgcc=0;
	int qcc_fc=0, gcc_fc=0;
	int dsub=0,rsub=0,csub=0;
	int q1=0,q2=0;
	int qq_events=0, gg_events=0, gq_events=0;
	float total_cs=0, qq_cs=0, gg_cs=0, gq_cs=0;
	bool dircheck=false, rescheck=false, ggcheck=false, qqcheck=false, gqcheck=false;
	//y = inelasticity = Beam energy (ben) / Photon Energy (gen)
	double ben, gen, y;

	// Initialize the generator.
	pythia.init();

	for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
		
		// Generate next event.
		if (!pythia.next()) continue;
		
		ev++;
		eventsid.push_back(iEvent);	//Pushing Back the event number

		if (iEvent == 0) {
			pythia.process.list();
			pythia.event.list();
		}

		double weight = info.weight();				//used to Normalize Histograms

		if (pythia.info.code()!=0){
			if (info.photonMode() == 3)								//Direct Process (281-284)
				{dircheck=true; dircount++;}	
			if (info.photonMode() == 1)								//Resolved Process
				{rescheck = true; rescount++;}	

			// Cross Section:
			float cs = info.sigmaGen();
			total_cs+=cs;
			int codechk = pythia.info.code();

			// qq scattering
			if (codechk == 112 || codechk == 114 || codechk == 116 || codechk == 281 || codechk == 282 || codechk == 283 || codechk == 285){
				qqcheck = true; 
				qq_cs+=cs;
				qq_events++;
				scattering.push_back(1);
			}
			// gg scattering
			if (codechk == 111 || codechk == 115){
				ggcheck = true; 
				gg_cs+=cs;
				gg_events++;
				scattering.push_back(2);
			}
			// gq scattering
			if (codechk == 113 || codechk == 284){
				gqcheck = true; 
				gq_cs+=cs;
				gq_events++;
				scattering.push_back(3);
			}
		}
		
		// ----------------------------------------------EVENT LOOP-----------------------------------------------------------
		// loop over all the particles created in an event
		for (int i = 0; i < pythia.event.size(); ++i){

			//Flags for Direct and Combined Events
			{
				if(dircheck) dirflag.push_back(1);
					else dirflag.push_back(0);
				if(rescheck) resflag.push_back(1);
					else resflag.push_back(0);
			}
			//Store incoming beam-photon energy (to calculate inelasticity):
			if (pythia.event[i].id()==22 && pythia.event[i].status()==-13) {
				gen = pythia.event[i].e();
			}
			//Skip for scattered Beam Electron   
			if(pythia.event[i].id()==-11 && pythia.event[i].mother1() == 0 && pythia.event[i].mother2() == 0){
				
				bpx.push_back(pythia.event[i].px());
				bpy.push_back(pythia.event[i].py());
				bpz.push_back(pythia.event[i].pz());
				beT.push_back(pythia.event[i].eT());
				be.push_back(pythia.event[i].e()); 
				
				ben = pythia.event[i].e();
				y=gen/ben;
				// cout << "gen: " << gen << " ben: " << ben << " inelastiity: " << y << endl;
				inelasticity.push_back(y);
				beame++; continue;

			}
			

			//Store 4-momentum of the final state charged particles (Combined Events)
			if ( pythia.event[i].isFinal() && pythia.event[i].isCharged()){ //    && pythia.event[i].isHadron()

				//Checking the number of constituents per event csub=dsub+rsub
				csub++;
				if(dircheck) dsub++;
				if(rescheck) rsub++;

				// Fill in the quark initiated process
				if(qqcheck){
					qcc++;
					qpx.push_back(pythia.event[i].px());
					qpy.push_back(pythia.event[i].py());
					qpz.push_back(pythia.event[i].pz());
					qe.push_back(pythia.event[i].e());
					qeT.push_back(pythia.event[i].eT());
				}

				// Fill in the gluon initiated process
				if(ggcheck){
					gcc++;
					gpx.push_back(pythia.event[i].px());
					gpy.push_back(pythia.event[i].py());
					gpz.push_back(pythia.event[i].pz());
					ge.push_back(pythia.event[i].e());
					geT.push_back(pythia.event[i].eT());
				}

				// Fill in the gluon initiated process
				if(qg){
					qgcc++;
					qgpx.push_back(pythia.event[i].px());
					qgpy.push_back(pythia.event[i].py());
					qgpz.push_back(pythia.event[i].pz());
					qge.push_back(pythia.event[i].e());
					qgeT.push_back(pythia.event[i].eT());
				}
				
				//Save the 4mom of all combined final state particles.
				px.push_back(pythia.event[i].px());
				py.push_back(pythia.event[i].py());
				pz.push_back(pythia.event[i].pz());
				e.push_back(pythia.event[i].e());
				eT.push_back(pythia.event[i].eT());

				//Check how many final state hadrons have Et>10
				if (pythia.event[i].eT() > 10.) 
					c++;

				// Fill Histograms :
				double pTch = pythia.event[i].pT();
				if (info.photonMode() == 1) pTres->Fill(pTch, weight);
				if (info.photonMode() == 3) pTdir->Fill(pTch, weight);

				// Counters to check for unecessary events :
				if (info.photonMode() == 2 || info.photonMode() == 4) chk++;	//Must be 0 for gamma-p collision.
				
			}
		}

		TC->Fill();

		// Clear Variables for next Run
		eventsid.clear(); scattering.clear();
		px.clear(); py.clear(); pz.clear(); e.clear(); eT.clear();	//Combined Process' 4 momentum
		bpx.clear(); bpy.clear(); bpz.clear(); be.clear(); beT.clear();	//Beam Remnant Electrons 4 momentum
		inelasticity.clear();
		qpx.clear();qpy.clear();qpz.clear();qe.clear();qeT.clear();
		gpx.clear();gpy.clear();gpz.clear();ge.clear();geT.clear();
		qgpx.clear();qgpy.clear();qgpz.clear();qge.clear();qgeT.clear();
		dirflag.clear(); resflag.clear();	//Int Flags to check for direct/resolved events
		qqflag.clear(); ggflag.clear();	//Int Flags to check for direct/resolved events
		dircheck = false, rescheck = false, ggcheck = false, qqcheck = false; //flags
		qq = false, qg = false, gg = false; //flags for subprocesses
	}
	// Show statistics after each run.
	pythia.stat();
	//cout << "Fraction of Resolved and Direct Subprocesses (based on their Subprocess Codes) for "<< nEvent << " Events :" << endl;
	cout << "\nTotal Selected Events : " << ev << endl;
	cout << "Total Beam Remnant electrons (should be same as Events) : " << beame << endl;

	cout << "\nTotal Combined (Res+Dir) Events : " << ev << endl;
	cout << "	Combined Constituents : " << csub << endl;
	cout << "Total Direct Events : " << dircount
			<< " | Fraction: " << fixed << (dircount/(double)ev)*100.0 << endl;
	cout << "	Direct Constituents : " << dsub << endl;
	cout << "Total Resolved Events : " << rescount
			<< " | Fraction: " << fixed << (rescount/(double)ev)*100.0 << endl;
	cout << "	Total Resolved Constituents : " << rsub << endl;

	//cout << fixed;	// To remove Scientific Notation.
	cout << "\n\nSelected Events: " << ev << endl;
	cout << fixed << "QQ Scattering: " << qq_events 
			<< "	|	QQ-Fraction of Total Events: " << fixed << (qq_events/(double)ev)*100.0 << endl;
	cout << fixed << "GG Scattering: " << gg_events
			<< "	|	GG-Fraction of Total Events: " << fixed << (gg_events/(double)ev)*100.0 << endl;
	cout << fixed << "GQ Scattering: " << gq_events
			<< "	|	GQ-Fraction of Total Events: " << fixed << (gq_events/(double)ev)*100.0 << endl;

	cout << "\nTotal Cross Section:" << scientific << (total_cs/ev) << endl;
	cout << scientific << "QQ Cross Section: " << (qq_cs/ev) << "	|	QQ-Fraction of Total Events: " 
			<< fixed << (qq_cs/total_cs)*100.0 << endl;
	cout << scientific << "GG Cross Section: " << (gg_cs/ev) << "	|	GG-Fraction of Total Events: "
			<< fixed << (gg_cs/total_cs)*100.0 << endl;
	cout << scientific << "GQ Cross Section: " << (gq_cs/ev) << "	|	GQ-Fraction of Total Events: "
			<< fixed << (gq_cs/total_cs)*100.0 << endl;

	cout << "\n\nCounters : " << endl;
	cout << "photonMode() == 2 or 4 (should be zero!): " << chk << endl;
	cout << "Final state particles with Et > 10Gev : " << c << endl;

	data->Write();
	data->Close();
	delete data;
}