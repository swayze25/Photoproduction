//	Event generation for PHOTPRODUCTION process
//	cd HEP_Projects/Project1_PhotoProduction/EventGen\ _Photoproduction/

// Authors: Ilkka Helenius <ilkka.m.helenius@jyu.fi>.

// Keywords: photon beam; photoproduction; photon-photon;
// In case of photon-photon interactions four different contributions are
// present:              ProcessType:
//  - resolved-resolved  1
//  - resolved-direct    2
//  - direct-resolved    3
//  - direct-direct      4

#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"
#include "Pythia8/Pythia.h"
#include <iostream>

using namespace Pythia8;

vector<int> remove_double(vector<int>&);

int main() {
	
// Number of events per run.
int nEvent = 500;

// Generator.
Pythia pythia;

// Decrease the output.
pythia.readString("Init:showMultipartonInteractions = on");
pythia.readString("Init:showChangedSettings = on");
pythia.readString("Init:showChangedParticleData = pff");
pythia.readString("Next:numberCount = 0");
pythia.readString("Next:numberShowInfo = 0");
pythia.readString("Next:numberShowProcess = 0");
pythia.readString("Next:numberShowEvent = 0");

// Shorthand for some public members of pythia (also static ones).
Settings& settings  = pythia.settings;
const Info& info = pythia.info;

// Beam parameters :
pythia.readString("Beams:frameType = 2");						// Beams of unequal energies
pythia.readString("Beams:idA = -11");							// Beam A : Electrons/Positrons
pythia.readString("Beams:idB = 2212");							// Beam B : Protons
pythia.readString("Beams:eA = 27.5");							// electron_En : 27.5 GeV
pythia.readString("Beams:eB = 820.");							// proton_En : 920 Gev
pythia.readString("Photon:Wmin  = 134.0");						// invariant mass_min of gamma-hadron pair
pythia.readString("Photon:Wmax  = 277.0");						// invariant mass_max of gamma-hadron pair
pythia.readString("PDF:beamA2gamma = on");						// Set PDF for beam photons

// Switch relevant processes on :
pythia.readString("HardQCD:all = on");								//For All Resolved Process
//pythia.readString("HardQCD:qq2qq = on");							//qq2qq - Quark Induced Events
//pythia.readString("HardQCD:gg2gg = on");							//gg2gg - Gluon Induced Events
//pythia.readString("HardQCD:qg2qg = on");							//Resolved qg2qg, mix

//For photon-parton interaction :
pythia.readString("PhotonParton:all = on");						// For Direct Process
//pythia.readString("PhotonParton:ggm2qqbar = on");				// Scattering g gamma → q qbar, (q = u,d,s) Code 271 (281).
//pythia.readString("PhotonParton:ggm2ccbar = on");				// Scattering g gamma → c cbar. Code 272 (282).
//pythia.readString("PhotonParton:ggm2bbbar = on");				// Scattering g gamma → b bbar. Code 273 (283).
//pythia.readString("PhotonParton:qgm2qg = on");				// Scattering q gamma → q g. Code 274 (284).
//pythia.readString("PhotonParton:qgm2qgm = on");				// Scattering q gamma → q gamma. Code 275 (285).

// Set Event Settings :
settings.mode("Photon:ProcessType", 0);							// Switch on : Automatic Mix
pythia.readString("PhaseSpace:pTHatMin = 5.0");					// Limit partonic pThat.
pythia.readString("MultipartonInteractions:pT0Ref = 3.0");		// Use tuned pT0ref for photon-hadron.
pythia.readString("Photon:Q2max = 1.0");						// Maximal Q2
//pythia.readString("SpaceShower:pTmaxMatch = 2");				// Allow emissions up to the kinematical limit


//ROOT for Storing Generated Events :
TFile *data = new TFile("data.root","recreate");
std::vector<Float_t> px,py,pz,e,eT,bpx,bpy,bpz,be,inelasticity;			//All combined-events 
std::vector<Int_t> dirflag, resflag, ggflag, qqflag, eventsid;			//Flags
std::vector<Float_t> qpx,qpy,qpz,qe,qeT;								//quark initiated particles
std::vector<Float_t> gpx,gpy,gpz,ge,geT;								//gluon initiated particles

//Create a tree for all Events :
TTree *TC = new TTree("combinedevents", "4-Momentum of all Charged Final State Particles (Combined)");

TC->Branch("eventid", "std::vector<Int_t>",  &eventsid);
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
TC->Branch("inelasticity", "std::vector<Float_t>",  &inelasticity);
//Set flags for direct/resolved/gluon/quark processes :
TC->Branch("dirflag", "std::vector<Int_t>",  &dirflag);
TC->Branch("resflag", "std::vector<Int_t>",  &resflag);
TC->Branch("ggflag", "std::vector<Int_t>",  &ggflag);
TC->Branch("qqflag", "std::vector<Int_t>",  &qqflag);

//4-Momentum of all quark-initiated final state charged particles :
TC->Branch("qpx", "std::vector<Float_t>",  &qpx);
TC->Branch("qpy", "std::vector<Float_t>",  &qpy);
TC->Branch("qpz", "std::vector<Float_t>",  &qpz);
TC->Branch("qe", "std::vector<Float_t>",  &qe);
TC->Branch("qeT", "std::vector<Float_t>",  &qeT);
//4-Momentum of all gluon-initiated final state charged particles :
TC->Branch("gpx", "std::vector<Float_t>",  &gpx);
TC->Branch("gpy", "std::vector<Float_t>",  &gpy);
TC->Branch("gpz", "std::vector<Float_t>",  &gpz);
TC->Branch("ge", "std::vector<Float_t>",  &ge);
TC->Branch("geT", "std::vector<Float_t>",  &geT);



// Parameters for histograms.
double pTmin = 0.0;
double pTmax = 5.0;
int nBinsPT  = 100;
//ROOT Histogram pT of Direct and Resolved Processes :
TH1F *pTres = new TH1F("pTresresR","Resolved-resolved contribution for pT distribution; pT  Res-Res ;Events ", nBinsPT, pTmin, pTmax);
TH1F *pTdir = new TH1F("pTdirresR","Direct-resolved contribution for pT distribution; pT  Dir-Res ;Events", nBinsPT, pTmin, pTmax);

//Counters :
int ev=0, rescount = 0, dircount = 0, beame = 0, qq=0, gg=0, qsub=0, gsub=0, gm=0, c=0, c1=0, c2=0, chk=0, qcc=0, gcc=0;
int qcc_fc=0, gcc_fc=0;
int dsub=0,rsub=0,csub=0;
bool dircheck = false, rescheck = false, qqcheck = false, ggcheck = false;
//y = inelasticity = Beam energy (ben) / Photon Energy (gen)
double ben, gen, y;
std::vector<int> d_quark, d_gluon;						// To store the initial daughters
std::vector<int> d_quark_all, d_gluon_all;		// To store all the daughters of the quarks
std::vector<int> motherData;									// Temporary vector to store all the data

// Initialize the generator.
pythia.init();

for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
	// Generate next event.
  if (!pythia.next()) continue;
	
	ev++;
	eventsid.push_back(iEvent);	//Pushing Back the event number

	// List the first process and event for each run.
	if (iEvent == 1) {
		pythia.process.list();
		pythia.event.list();
	}

	double weight = info.weight();				//used to Normalize Histograms

	if (pythia.info.code()!=0){
	
		if (info.photonMode() == 3)								//Direct Process (281-284)
			{dircheck=true; dircount++;}	
		if (info.photonMode() == 1)								//Resolved Process
			{rescheck = true; rescount++;}	

		if (pythia.info.code() == 114 || pythia.info.code() == 281)			//Quark inititated process
			{qqcheck = true; qq++;}
		if (pythia.info.code() == 111) 																	//Gluon inititated process
			{ggcheck = true; gg++;}
		
		//Filling the int-flags with 1 if true, 0 is false
		if(dircheck) dirflag.push_back(1);
			else dirflag.push_back(0);
		if(rescheck) resflag.push_back(1);
			else resflag.push_back(0);
		if(qqcheck) qqflag.push_back(1);
			else qqflag.push_back(0);
		if(ggcheck) ggflag.push_back(1);
			else ggflag.push_back(0);
	}
	
	// Loop over event record :
	for (int i = 0; i < pythia.event.size(); ++i){			
		


		//----------------------------------------QUARK INDUCED FINAL STATE CHARGED PARTICLES----------------------------------------

		// Select the "First" split quark (directly from the beam particle)
		if (pythia.event[i].isQuark() && pythia.event[i].mother1() <= 4 && pythia.event[i].mother2() == 0){
			// Find first generation of daughters, store in 'd_quark'
			d_quark = pythia.event.daughterList(i);
			// Recursively add daughters of unstable particles.
			int size = d_quark.size();
			for (int iDau = 0; iDau < size; ++iDau){
			Particle& partNow = pythia.event[d_quark[iDau]];
			if (!partNow.isFinal()){
				vector<int> grandDauVec = partNow.daughterList();
				for (int i = 0; i < int(grandDauVec.size()); ++i)
					d_quark.push_back(grandDauVec[i]);
					size += grandDauVec.size();
				}
			}

			// // Remove the doubly counted daughter particles due to recursive process :
			cout << "\n->Total Daughters of mother quark (added recursively) - " << d_quark.size();
			d_quark = remove_double(d_quark);
			cout << " | After removing double count - " << d_quark.size() << endl;

			cout << "Quark induced final state particles for Event " << iEvent << " (" << pythia.info.name() << ")" << endl;
			cout << "First-beam-split Quark: " << pythia.event[i].nameWithStatus() << "::" 
				<< pythia.event[i].index() << "::"
				<< pythia.event[i].status() << endl;
			
			// Access the list of all daughters and separate the charged and final particles
			for(int c : d_quark){
				Particle& partChk = pythia.event[c];
				if (partChk.isFinal() && partChk.isCharged() ){ // && partChk.status() >= 81 && partChk.status() <= 89){
					qcc++;	//count the number of quark initiated particles & store the 4Mom of the particle
					qpx.push_back(partChk.px());
					qpy.push_back(partChk.py());
					qpz.push_back(partChk.pz());
					qe.push_back(partChk.e());
					qeT.push_back(partChk.eT());

					// Access the 'first' mother of the daughter particle (IS or FS showers):
					motherData = partChk.motherList();
					// int firstmother = motherData.back();
					cout << "\nFinal_Charged Particle (Q): " << pythia.event[c].nameWithStatus() << "::" 
						 << pythia.event[c].index() << "::"
						 << pythia.event[c].status();
					cout << "\nMother List: "; for(int r : motherData) cout<< pythia.event[r].nameWithStatus() << "::"
						 << pythia.event[r].index() << "::"
						 << pythia.event[r].status() << "		";
					motherData.clear();
				}
			}
			d_quark.clear();
		}
		

		//----------------------------------------GLUON INDUCED FINAL STATE CHARGED PARTICLES----------------------------------------
		// Select the "First" quark (for Initial State Showers)
		if (pythia.event[i].isGluon() && pythia.event[i].mother1() <= 4 && pythia.event[i].mother2() == 0){
			// Find first generation of daughters, store in 'd_gluon'
			
			d_gluon = pythia.event.daughterList(i);
			// Recursively add daughters of unstable particles.
			int size = d_gluon.size();
			for (int iDau = 0; iDau < size; ++iDau){
			Particle& partNow = pythia.event[d_gluon[iDau]];
			if (!partNow.isFinal()){
				vector<int> grandDauVec = partNow.daughterList();
				for (int i = 0; i < int(grandDauVec.size()); ++i)
					d_gluon.push_back(grandDauVec[i]);
					size += grandDauVec.size();
				}
			}

			// Remove the doubly counted daughter particles due to recursive process :
			cout << "\n->Total Daughters of mother quark (added recursively) - " << d_gluon.size();
			d_gluon = remove_double(d_gluon);
			cout << " | After removing double count - " << d_gluon.size() << endl;

			// cout << "Gluon induced final state particles for Event " << iEvent << " (" << pythia.info.name() << ")" << endl;
			// cout << "First-beam-split Gluon: " << pythia.event[i].nameWithStatus() << "::" 
			// 	<< pythia.event[i].index() << "::"
			// 	<< pythia.event[i].status() << endl;
			
			// Access the list of daughters and separate the charged and final states
			for(int c : d_gluon){
				Particle& partChk = pythia.event[c];  // checking for all daughters
				if (partChk.isFinal() && partChk.isCharged()){	//&& partChk.status() >= 81 && partChk.status() <= 89){
					gcc++;
					gpx.push_back(partChk.px());
					gpy.push_back(partChk.py());
					gpz.push_back(partChk.pz());
					ge.push_back(partChk.e());
					geT.push_back(partChk.eT());

					// Access the first mother of the daughter particle :
					// motherData = partChk.motherList();
					// int firstmother = motherData.back();
						// cout << "\nFinal_Charged Particle (G): " << pythia.event[c].nameWithStatus() << "::" 
						// 	<< pythia.event[c].index() << "::"
						// 	<< pythia.event[c].status();
						// cout << "\nMother List: "; for(int r : motherData) cout<< pythia.event[r].nameWithStatus() << "::"
						// 	<< pythia.event[r].index() << "::"
						// 	<< pythia.event[r].status()<< "		";
						// motherData.clear();
				}
			}
			d_gluon.clear();
		}
		

		//Store incoming beam-photon energy (to calculate inelasticity):
		if (pythia.event[i].id()==22 && pythia.event[i].status()==-13) {
			gen = pythia.event[i].e();
		}

		//Skip for scattered Beam Electron   
		if(pythia.event[i].id()==-11 && pythia.event[i].mother1() == 1 && pythia.event[i].mother2() == 0){
			
			bpx.push_back(pythia.event[i].px());
			bpy.push_back(pythia.event[i].py());
			bpz.push_back(pythia.event[i].pz());
			be.push_back(pythia.event[i].e()); 
			
			ben = pythia.event[i].e();
			y=gen/ben;
			inelasticity.push_back(y);
			beame++; continue;

		}

		//Counters for number of Direct and Resolved Phontons :
		if (info.photonMode()==3 && pythia.event[i].id()==22 && pythia.event[i].status()==-13)
			c1++;
		if (info.photonMode()==1 && pythia.event[i].id()==22 && pythia.event[i].status()==-13)
			c2++;

		//Store 4-momentum of the final state charged particles (Combined Events)
		if ( pythia.event[i].isFinal() && pythia.event[i].isCharged()){ //    && pythia.event[i].isHadron()
			
			// cout << "\n\n\nfor Event " << iEvent << " (" << pythia.info.name() << ") " << endl;
			// Particle& partChk = pythia.event[i];
			// motherData = partChk.motherList();
			// cout<< "Final_Charged Particle: " << partChk.nameWithStatus() << "::" 
			// 	<< partChk.index() << "::"
			// 	<< partChk.status() << "\nMotherList : ";
			// for (int ck : motherData) {
			// 	cout << pythia.event[ck].nameWithStatus() << "::"
			// 			<< pythia.event[ck].index() << "::"
			// 			<< pythia.event[ck].status() << "	";
			// }
			// int firstmother = motherData.back();
			// if (pythia.event[firstmother].isGluon()){
			// 	gcc_fc++;
			// 	gpx.push_back(partChk.px());
			// 	gpy.push_back(partChk.py());
			// 	gpz.push_back(partChk.pz());
			// 	ge.push_back(partChk.e());
			// 	geT.push_back(partChk.eT());
			// }
			// if (pythia.event[firstmother].isQuark()){
			// 	qcc_fc++;
			// 	qpx.push_back(partChk.px());
			// 	qpy.push_back(partChk.py());
			// 	qpz.push_back(partChk.pz());
			// 	qe.push_back(partChk.e());
			// 	qeT.push_back(partChk.eT());
			// }

			//Save the 4mom of all combined final state particles.
			px.push_back(pythia.event[i].px());
			py.push_back(pythia.event[i].py());
			pz.push_back(pythia.event[i].pz());
			e.push_back(pythia.event[i].e());
			eT.push_back(pythia.event[i].eT());

			//Check how many final state hadrons have Et>10
			if (pythia.event[i].eT() > 10.) 
				c++;
			
			//Checking the number of constituents per event cs=dsub+rsub
			csub++;
			if(dircheck) dsub++;
			if(rescheck) rsub++;
			if(ggcheck) gsub++;
			if(qqcheck) qsub++;

			// Store the pT value of the event
	    	double pTch = pythia.event[i].pT();
			if (info.photonMode() == 1) pTres->Fill(pTch, weight);
			if (info.photonMode() == 3) pTdir->Fill(pTch, weight);
			if (info.photonMode() == 2 || info.photonMode() == 4) chk++;	//Must be 0 for gamma-p collision.
			
		}
	}
	
	TC->Fill();

	eventsid.clear();
	px.clear(); py.clear(); pz.clear(); e.clear(); eT.clear();	//Combined Process' 4 momentum
	bpx.clear(); bpy.clear(); bpz.clear(); be.clear();	//Beam Remnant Electrons 4 momentum
	inelasticity.clear();
	qpx.clear();qpy.clear();qpz.clear();qe.clear();qeT.clear();
	gpx.clear();gpy.clear();gpz.clear();ge.clear();geT.clear();
	dirflag.clear(); resflag.clear();	//Int Flags to check for direct/resolved events
	qqflag.clear(); ggflag.clear();	//Int Flags to check for direct/resolved events
	dircheck = false, rescheck = false, qqcheck = false, ggcheck = false; //flags
}
// Show statistics after each run.
pythia.stat();
//cout << "Fraction of Resolved and Direct Subprocesses (based on their Subprocess Codes) for "<< nEvent << " Events :" << endl;
cout << "\nTotal Selected Events : " << ev << endl;
cout << "Total beam remnant electrons (should be same as Events) : " << beame << endl;

cout << "\nTotal Direct Events : " << dircount << endl;
cout << "Total Direct Constituents : " << dsub << endl;
cout << "Total Resolved Events : " << rescount << endl;
cout << "Total Resolved Constituents : " << rsub << endl;
cout << "Total Combined Events : " << ev << endl;
cout << "Total Combined Constituents : " << csub << endl;

cout << "\nTotal quark initiated events (Code 114) : " << qq << endl;
cout << "Total quark initiated charged final state particles (qqbar2qqbar) : " << qsub << endl;
cout << "Total gluon initiated events (Code 111) : " << gg << endl;
cout << "Total gluon initiated charged final state particles (gg2gg) : " << gsub << endl;

cout << "\nCounters : " << endl;
cout << "photonMode() == 2 or 4 : " << chk << endl;
cout << "Final state particles with Et > 10Gev : " << c << endl;
cout << "Photons when photonMode() == 3 : " << c1 << endl;
cout << "Photons when photonMode() == 1 : " << c2 << endl;

// cout << "\nHow to read event (particleName)::(index)::(status)";
cout << "\nQuark Initiated Process : " << qcc << endl;
cout << "Gluon Initiated Process : " << gcc << endl;
cout << "\nQuark Initiated Process via mothercheck : " << qcc_fc << endl;
cout << "Gluon Initiated Process via mothercheck : " << gcc_fc << endl;



data->Write();
data->Close();
delete data;
}

vector<int> remove_double(vector<int>& vec){
  //cout << "\nbefore :"; for (int i : vec) cout << i << "  ";
std::sort(vec.begin(), vec.end());
vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
  //cout << "\nafter :"; for (int i : vec) cout << i << "  ";
  return vec;
}