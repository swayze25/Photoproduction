//*******************************************
// Photoproduction - Creates Jets and SubJets
//
// Siddharth Singh
// Manjit Kaur, Ritu Aggarwal
//
// Creating Jets/Subjets at HERA energies
//
// cd /home/siddharth/HEP_Projects/Project1_PhotoProduction/Jets_Photoproduction
//*******************************************

#include "fastjet/config.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "TH1.h"
#include <TF2.h>
#include "TH2F.h"
#include "TH3F.h"
#include "THStack.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <TGraph2D.h>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <TLorentzVector.h>
#include <iostream> 
#include <sstream>
#include <vector> 
#include <numeric>
#include <cstdio>
#include <fstream>
using namespace fastjet;
using namespace std;

//Declaring the Functions :
void combined_jets (const vector<fastjet::PseudoJet> &);								// Jet for Combined (Res+Dir) Processes
void direct_jets (const vector<fastjet::PseudoJet> &);									// Jet for Direct Processes
void resolved_jets (const vector<fastjet::PseudoJet> &);								// Jet for Resolved Processes
void quark_jets (const vector<fastjet::PseudoJet> &);										// Jets for QQ scatterings
void gluon_jets (const vector<fastjet::PseudoJet> &);										// Jets for GG scatterings
void print_jets (const vector<fastjet::PseudoJet> &, int flag);					// Printing Jets and Subjets
void jet_visualization (const vector<fastjet::PseudoJet> & , int);			// Visualizing Jets, flag=1 for Jet, 2 for SubJet
double jet_shape_diff (const vector<fastjet::PseudoJet> &);							// Differential Jet Shape
double jet_shape_incl (const vector<fastjet::PseudoJet> &);							// Inclusive Jet Shape, for quark and gluon jets
double x_obs_dijets (const vector<fastjet::PseudoJet> &);								// x_observable for diff res and dir jets		
double average(const vector<double> &);																	// used to count the average vaue of an array

// Boolean variables to know which process is being read
bool qq, gg, qg;

//Set the respective parameters :
double EtMin_C=15.0, R_C=1.0, crapmin=-1.0, crapmax=2., ycut_C = 0.0005, dcut_C;				
double EtMin_D=14.0, R_D=1.0, drapmin=-1.0, drapmax=2.5, ycut_D = 0.0005, dcut_D;				
double EtMin_R=14.0, R_R=1.0, rrapmin=-1.0, rrapmax=2.5, ycut_R = 0.0005, dcut_R;				
double EtMin_Q=14.0, R_Q=1.0, qrapmin=-1.0, qrapmax=2.5, ycut_Q = 0.0005, dcut_Q;
double EtMin_G=14.0, R_G=1.0, grapmin=-1.0, grapmax=2.5, ycut_G = 0.0005, dcut_G;
double EtMin_QG=14.0, R_QG=1.0, qgrapmin=-1.0, qgrapmax=2.5, ycut_QG = 0.0005, dcut_QG;

// differential jet shape
double rr_diff = 0.15;
vector<double> rho_all;

// inclusive jet shape
double rr_incl = 0.3;
vector<double> psi_all;

// Define variables to store the mean subjet multiplicity vs ycut
vector<double> n_C;
vector<double> n_R;
vector<double> n_D;
vector<double> n_Q;
vector<double> n_G;

//Reading  the input root file.
TFile *myFile = new TFile("data.root", "READ");
// For Events :
TTreeReader myReader("combinedevents", myFile);
TTreeReaderArray<Int_t> eventid(myReader, "eventid");
TTreeReaderArray<Float_t> px(myReader, "px");
TTreeReaderArray<Float_t> py(myReader, "py");
TTreeReaderArray<Float_t> pz(myReader, "pz");
TTreeReaderArray<Float_t> e(myReader, "e");
// Reading Beam Electrons & inelasticity:
TTreeReaderArray<Float_t> bpx(myReader, "bpx");
TTreeReaderArray<Float_t> bpy(myReader, "bpy");
TTreeReaderArray<Float_t> bpz(myReader, "bpz");
TTreeReaderArray<Float_t> be(myReader, "be");
TTreeReaderArray<Float_t> beT(myReader, "beT");
TTreeReaderArray<Float_t> inelasticity(myReader, "inelasticity");
// Jet constituents for final state QQ scatterings:
TTreeReaderArray<Float_t> qpx(myReader, "qpx");
TTreeReaderArray<Float_t> qpy(myReader, "qpy");
TTreeReaderArray<Float_t> qpz(myReader, "qpz");
TTreeReaderArray<Float_t> qe(myReader, "qe");
TTreeReaderArray<Float_t> qeT(myReader, "qeT");
// Jet constituents for final state GG scatterings:
TTreeReaderArray<Float_t> gpx(myReader, "gpx");
TTreeReaderArray<Float_t> gpy(myReader, "gpy");
TTreeReaderArray<Float_t> gpz(myReader, "gpz");
TTreeReaderArray<Float_t> ge(myReader, "ge");
TTreeReaderArray<Float_t> geT(myReader, "geT");
// Jet constituents for final state QG scatterings:
TTreeReaderArray<Float_t> qgpx(myReader, "qgpx");
TTreeReaderArray<Float_t> qgpy(myReader, "qgpy");
TTreeReaderArray<Float_t> qgpz(myReader, "qgpz");
TTreeReaderArray<Float_t> qge(myReader, "qge");
TTreeReaderArray<Float_t> qgeT(myReader, "qgeT");
// Reading flags for direct or resolved :
TTreeReaderArray<Int_t> dirflag(myReader, "dirflag");
TTreeReaderArray<Int_t> resflag(myReader, "resflag");
// TTreeReaderArray<Int_t> qqflag(myReader, "qqflag");
// TTreeReaderArray<Int_t> ggflag(myReader, "ggflag");

// Writing an Output ROOT file:
TFile *jets = new TFile("jets.root", "RECREATE");

// COMBINED EVENTS :
std::vector<Float_t> pxcj, pycj, pzcj, ecj, etcj, phicj, etacj, ncjets;
std::vector<Double_t> mean_js_incl, mean_js_diff;
TTree *jetcombined = new TTree("jetcombined", "");
TH1F *ncJets = new TH1F("ncJets", "#Combined Jets; Number of Combined Jets; Events", 10, 0.5, 10);
TH1F *ncsubJets = new TH1F("ncsubJets", "#Combined <mean subJet multiplicity>; Combined nSubjets; Jets", 10, 0.5, 10);


// visualizing jets on the phi-eta plane
int ev1=58494,ev2=74578,ev3=95290; int c1=0,c2=0,c3=0; // Select events for Jet Vis
TH2F *phi_eta = new TH2F("phi_eta", "#eta vs #phi (Res+Dir)", 100, -6, 6, 100, 0, 6.3);
// Store 3 Jets along with their constituents for visualisation.
TH2F *jetvis1 = new TH2F("jetvis1", "", 100, -3, 3, 100, 0, 6.3);
TH2F *jetvis2 = new TH2F("jetvis2", "", 100, -3, 3, 100, 0, 6.3);
TH2F *jetvis3 = new TH2F("jetvis3", "", 100, -3, 3, 100, 0, 6.3);

// RESOLVED EVENTS :
std::vector<Float_t> pxrj, pyrj, pzrj, erj, etrj, phirj, etarj, nrjets;
TTree *jetresolved = new TTree("jetresolved", "");
TH1F *nrJets = new TH1F("nrJets", "#Resolved Jets; Number of Resolved Jets; Events", 10, 0.5, 10);
TH1F *nrsubJets = new TH1F("nrsubJets", "#Resolved <mean subJet multiplicity>; Resolved nSubjets; Jets", 10, 0.5, 10);

// DIRECT EVENTS :
std::vector<Float_t> pxdj, pydj, pzdj, edj, etdj, phidj, etadj, ndjets;
TTree *jetdirect = new TTree("jetdirect", "");
TH1F *ndJets = new TH1F("ndJets", "#Direct Jets; Number of Direct Jets; Events", 10, 0.5, 10);
TH1F *ndsubJets = new TH1F("ndsubJets", "#Direct <mean subJet multiplicity>; Direct nSubjets; Jets", 10, 0.5, 10);

// QUARK INDUCED EVENTS :
std::vector<Float_t> pxqj, pyqj, pzqj, eqj, etqj, phiqj, etaqj, nqjets;
TTree *jetquark = new TTree("jetquark", "");
TH1F *nqJets = new TH1F("nqJets", "#Quark Initiated Jets; Number of Quark Initiated Jets; Events", 10, 0.5, 10);
std::vector<Float_t> pxsubqj, pysubqj, pzsubqj, esubqj, etsubqj, phisubqj, etasubqj, nsubqjets;
TTree *subjetquark = new TTree("subjetquark", "");
TH1F *nqsubJets = new TH1F("nqsubJets", "#quark Initiated <mean subJet multiplicity>; nq_subJets; Jets", 10, 0.5, 10);

//GLUON INDUCED EVENTS :
std::vector<Float_t> pxgj,pygj,pzgj,egj,etgj,phigj,etagj,ngjets;
TTree *jetgluon = new TTree("jetgluon", "");
TH1F *ngJets = new TH1F("ngJets","#Gluon Initiated Jets; Number of Gluon Initiated Jets; Events",10,0.5,10);
std::vector<Float_t> pxsubgj,pysubgj,pzsubgj,esubgj,etsubgj,phisubgj,etasubgj,nsubgjets;
TTree *subjetgluon = new TTree("subjetgluon", "");
TH1F *ngsubJets = new TH1F("ngsubJets","#Gluon Initiated <mean subJet multiplicity>; nq_Initiated subJets; Jets",10,0.5,10);

int ncombjets=0, ndirjets=0, nresjets=0, nquarkjets=0, ngluonjets=0;					// Number of Jets
int ncombsubjets=0, ndirsubjets=0, nressubjets=0, nquarksubjets=0, ngluonsubjets=0;		// Number of subJets
int ev=0;		// Number of Events
double in_y, beam_eT;

int main(){

	vector<fastjet::PseudoJet> input_particles_combined;
	vector<fastjet::PseudoJet> input_particles_direct;
	vector<fastjet::PseudoJet> input_particles_resolved;
	vector<fastjet::PseudoJet> input_particles_quark;
	vector<fastjet::PseudoJet> input_particles_gluon;
	vector<fastjet::PseudoJet> input_particles_qg;
	
	//int ev=0;	
	int c=0, d=0, r=0, cs=0, ds=0, rs=0, g=0, q=0, qs=0, gs=0;
	
	//------------------------------------DECLARING JET PARAMETERS------------------------------------
	//COMBINED ENTRIES:
	jetcombined->Branch("combined_px", "std::vector<Float_t>",  &pxcj);
	jetcombined->Branch("combined_py", "std::vector<Float_t>",  &pycj);
	jetcombined->Branch("combined_pz", "std::vector<Float_t>",  &pzcj);
	jetcombined->Branch("combined_E", "std::vector<Float_t>",  &ecj);
	jetcombined->Branch("combined_Et", "std::vector<Float_t>",  &etcj);
	jetcombined->Branch("combined_phi", "std::vector<Float_t>",  &phicj);
	jetcombined->Branch("combined_eta", "std::vector<Float_t>",  &etacj);
	jetcombined->Branch("combined_njets", "std::vector<Float_t>",  &ncjets);
	jetcombined->Branch("combined_mean_js_incl", "std::vector<Double_t>",  &mean_js_incl);
	jetcombined->Branch("combined_mean_js_diff", "std::vector<Double_t>",  &mean_js_diff);
	
	//DIRECT ENTRIES:
	jetdirect->Branch("direct_px", "std::vector<Float_t>",  &pxdj);
	jetdirect->Branch("direct_py", "std::vector<Float_t>",  &pydj);
	jetdirect->Branch("direct_pz", "std::vector<Float_t>",  &pzdj);
	jetdirect->Branch("direct_E", "std::vector<Float_t>",  &edj);
	jetdirect->Branch("direct_Et", "std::vector<Float_t>",  &etdj);
	jetdirect->Branch("direct_phi", "std::vector<Float_t>",  &phidj);
	jetdirect->Branch("direct_eta", "std::vector<Float_t>",  &etadj);
	jetdirect->Branch("direct_njets", "std::vector<Float_t>",  &ndjets);
	//RESOLVED ENTRIES:
	jetresolved->Branch("resolved_px", "std::vector<Float_t>",  &pxrj);
	jetresolved->Branch("resolved_py", "std::vector<Float_t>",  &pyrj);
	jetresolved->Branch("resolved_pz", "std::vector<Float_t>",  &pzrj);
	jetresolved->Branch("resolved_E", "std::vector<Float_t>",  &erj);
	jetresolved->Branch("resolved_Et", "std::vector<Float_t>",  &etrj);
	jetresolved->Branch("resolved_phi", "std::vector<Float_t>",  &phirj);
	jetresolved->Branch("resolved_eta", "std::vector<Float_t>",  &etarj);
	jetresolved->Branch("resolved_njets", "std::vector<Float_t>",  &nrjets);
	
	
	//QUARK INITIATED ENTRIES QQ:
	jetquark->Branch("q_px", "std::vector<Float_t>",  &pxqj);
	jetquark->Branch("q_py", "std::vector<Float_t>",  &pyqj);
	jetquark->Branch("q_pz", "std::vector<Float_t>",  &pzqj);
	jetquark->Branch("q_E", "std::vector<Float_t>",  &eqj);
	jetquark->Branch("q_Et", "std::vector<Float_t>",  &etqj);
	jetquark->Branch("q_phi", "std::vector<Float_t>",  &phiqj);
	jetquark->Branch("q_eta", "std::vector<Float_t>",  &etaqj);
	jetquark->Branch("q_njets", "std::vector<Float_t>",  &nqjets);
	//SubJets Quarks induced process :
	subjetquark->Branch("px_subjets_quark", "std::vector<Float_t>",  &pxsubqj);
	subjetquark->Branch("py_subjets_quark", "std::vector<Float_t>",  &pysubqj);
	subjetquark->Branch("pz_subjets_quark", "std::vector<Float_t>",  &pzsubqj);
	subjetquark->Branch("E_subjets_quark", "std::vector<Float_t>",  &esubqj);
	subjetquark->Branch("Et_subjets_quark", "std::vector<Float_t>",  &etsubqj);
	subjetquark->Branch("eta_subjets_quark", "std::vector<Float_t>",  &etasubqj);
	subjetquark->Branch("phi_subjets_quark", "std::vector<Float_t>",  &phisubqj);
	subjetquark->Branch("nsjets_quark", "std::vector<Float_t>",  &nsubqjets);
	
	
	//GLUON INITIATED ENTRIES GG:
	jetgluon->Branch("g_px", "std::vector<Float_t>",  &pxgj);
	jetgluon->Branch("g_py", "std::vector<Float_t>",  &pygj);
	jetgluon->Branch("g_pz", "std::vector<Float_t>",  &pzgj);
	jetgluon->Branch("g_E", "std::vector<Float_t>",  &egj);
	jetgluon->Branch("g_Et", "std::vector<Float_t>",  &etgj);
	jetgluon->Branch("g_phi", "std::vector<Float_t>",  &phigj);
	jetgluon->Branch("g_eta", "std::vector<Float_t>",  &etagj);
	jetgluon->Branch("g_njets", "std::vector<Float_t>",  &ngjets);
	//SubJets gluons induced process :
	subjetgluon->Branch("px_subjets_gluon", "std::vector<Float_t>",  &pxsubgj);
	subjetgluon->Branch("py_subjets_gluon", "std::vector<Float_t>",  &pysubgj);
	subjetgluon->Branch("pz_subjets_gluon", "std::vector<Float_t>",  &pzsubgj);
	subjetgluon->Branch("E_subjets_gluon", "std::vector<Float_t>",  &esubgj);
	subjetgluon->Branch("Et_subjets_gluon", "std::vector<Float_t>",  &etsubgj);
	subjetgluon->Branch("eta_subjets_gluon", "std::vector<Float_t>",  &etasubgj);
	subjetgluon->Branch("phi_subjets_gluon", "std::vector<Float_t>",  &phisubgj);
	subjetgluon->Branch("nsjets_gluon", "std::vector<Float_t>",  &nsubgjets);
	//-------------------------------------------------------------------------------------------------
	
	//while loop to read the input particles and call Jets :
	while (myReader.Next()){
	
		for (int i = 0; i < eventid.GetSize(); i++){
			ev++;
			in_y = inelasticity[ev];
			beam_eT = beT[ev];
			qq = false, gg = false, qg = false;

			//------------------------------------COMBINED INPUT PARTICLES------------------------------------
			for (int x = 0; x < px.GetSize(); x++){
				
				input_particles_combined.push_back(fastjet::PseudoJet(px[x],py[x],pz[x],e[x]));
				
				cs++;
			}
			//------------------------------------DIRECT INPUT PARTICLES------------------------------------
			if (dirflag[i]>0){
				for (int x = 0; x < px.GetSize(); x++){
						input_particles_direct.push_back(fastjet::PseudoJet(px[x],py[x],pz[x],e[x]));
						ds++;
					}
				d++;
			}
			//------------------------------------RESOLVED INPUT PARTICLES------------------------------------
			if (resflag[i]>0){
				for (int x = 0; x < px.GetSize(); x++){
						input_particles_resolved.push_back(fastjet::PseudoJet(px[x],py[x],pz[x],e[x]));
						rs++;
					}
				r++;
			}

			//------------------------------------QQ SCATTERED PARTICLES------------------------------------
			for (int x = 0; x < qpx.GetSize(); x++){
				input_particles_quark.push_back(fastjet::PseudoJet(qpx[x],qpy[x],qpz[x],qe[x]));
				qs++;
			}
			//------------------------------------GG SCATTERED PARTICLES------------------------------------
			for (int x = 0; x < gpx.GetSize(); x++){
				input_particles_gluon.push_back(fastjet::PseudoJet(gpx[x],gpy[x],gpz[x],ge[x]));
				gs++;
			}
			//------------------------------------QG SCATTERED PARTICLES------------------------------------
			for (int x = 0; x < qgpx.GetSize(); x++){
				input_particles_qg.push_back(fastjet::PseudoJet(qgpx[x],qgpy[x],qgpz[x],qge[x]));
			}
			
			//Calling Functions for Jets/Subjets :
			combined_jets(input_particles_combined);
			resolved_jets(input_particles_resolved);
			direct_jets(input_particles_direct);
			quark_jets(input_particles_quark);
			gluon_jets(input_particles_gluon);
			
			//Clearing Particles for the next run :
			input_particles_combined.clear();
			input_particles_resolved.clear();
			input_particles_direct.clear();
			input_particles_quark.clear();
			input_particles_gluon.clear();
			input_particles_qg.clear();
		}
  }

  // a "header" for the output
  //cout << "Ran " << jet_def.description() << endl;
	cout << "\n=======================================================================================" << endl;
	cout << "\nNumber of Combined Events Read : " << ev << " (Constituents : " << cs << ")" << endl;
  cout << "Number of Resolved Events : " << r << " (Constituents : " << rs << ")" << endl;
  cout << "Number of Direct Events : " << d << " (Constituents : " << ds << ")" << endl;
  cout << "Number of Quark initiated Events : " << q << " (Constituents : " << qs << ")" << endl;
  cout << "Number of Gluon initiated Events : " << g << " (Constituents : " << gs << ")" << endl << endl;
  //cout << "Strategy adopted by FastJet was "<< clust_seq.strategy_string()<<endl;
  cout << "\n=======================================================================================" << endl;
  
  cout << "\n#Jets formed (with Combined Jet Constituents : " << cs << ") : " << 	ncombjets << endl;
  cout << "		Minimum Etmin Combined: " << EtMin_C << " GeV" << endl;
  cout << "		R-Combined: " << R_C << endl;
  cout << "		PseudoRapidity Range : " << crapmin << " to " << crapmax << endl;
  cout << "		Mean SubJet Multiplicity for Combined Jets : " << ncombsubjets << endl;
  
  cout << "\n#Jets formed (with Resolved Jet Events, Constituents : " << rs << ") : " << 	nresjets << endl;
  cout << "		Minimum Etmin Resolved: " << EtMin_R << " GeV" << endl;
  cout << "		R-Resolved: " << R_R << endl;
  cout << "		PseudoRapidity Range : " << rrapmin << " to " << rrapmax << endl;
	cout << "		Mean SubJet Multiplicity for Resolved Jets : " << nressubjets << endl;
  
  cout << "\n#Jets formed (with Direct Jet Events) Constituents : " << ds << ") : " << 	ndirjets << endl;
  cout << "		Minimum Etmin Direct: " << EtMin_D << " GeV" << endl;
  cout << "		R-Direct: " << R_D << endl;
  cout << "		PseudoRapidity Range : " << drapmin << " to " << drapmax << endl;
	cout << "		Mean SubJet Multiplicity for Direct Jets : " << ndirsubjets << endl;
  cout << "\n---------------------------------------------------------------------------------------" << endl;
  cout << "\n#Jets formed (with quark initiated processes, Constituents : " << qs << ") : " << 	nquarkjets << endl;
  cout << "		Minimum Etmin Quark: " << EtMin_Q << " GeV" << endl;
  cout << "		R-quark: " << R_Q << endl;
  cout << "		PseudoRapidity Range : " << qrapmin << " to " << qrapmax << endl;
  cout << "#Subjets formed (with quark initiated processes) : " << nquarksubjets << endl;

  cout << "\n#Jets formed (with gluon initiated processes, Constituents : " << gs << ") : " << 	ngluonjets << endl;
  cout << "		Minimum Etmin Gluon: " << EtMin_G << " GeV" << endl;
  cout << "		R-gluon: " << R_G << endl;
  cout << "		PseudoRapidity Range : " << grapmin << " to " << grapmax << endl;
	cout << "#Subjets formed (with gluon initiated processes) : " << ngluonsubjets << endl;

	cout << "\n=======================================================================================" << endl;

	// Mean Differential Jet Size 
	double temp = average(rho_all);
	// Mean Subjets 
	double n_C_mean = average(n_C);
	double n_D_mean = average(n_D);
	double n_R_mean = average(n_R);
	double n_Q_mean = average(n_Q);
	double n_G_mean = average(n_G);

	cout << "Differential Jet size for r(" << rr_diff << ") EtJetMin(" <<  EtMin_C << ") :: " << temp << endl;
	cout << "Combined Jets	:	ycut (x-axis) = " << ycut_C << " :: <n_C> (y-axis) = " << n_C_mean << endl;
	cout << "Direct Jets	:	ycut (x-axis) = " << ycut_D << " :: <n_D> (y-axis) = " << n_D_mean << endl;
	cout << "Resolved Jets	:	ycut (x-axis) = " << ycut_R << " :: <n_R> (y-axis) = " << n_R_mean << endl;
	cout << "Quark Jets	:	ycut (x-axis) = " << ycut_Q << " :: <n_Q> (y-axis) = " << n_Q_mean << endl;
	cout << "Gluon Jets	:	ycut (x-axis) = " << ycut_G << " :: <n_G> (y-axis) = " << n_G_mean << endl;
	
	// cout << s << endl;
	// cout << c1 << " " << c2 << " " << c3;
  jets->Write();
	jets->Close();
	delete jets;
}


//------------------------------------COMBINED JETS-------------------------------------------
void combined_jets (const vector<fastjet::PseudoJet> & input_particles_combined){
  
	fastjet::JetDefinition jet_def(kt_algorithm,R_C);
	fastjet::ClusterSequence clust_seq(input_particles_combined, jet_def);
	Selector jet_selector = SelectorEtMin(EtMin_C) && SelectorRapRange(crapmin, crapmax);
	vector<PseudoJet> inclusive_jets_combined = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));

	if(ev == ev1 || ev== ev2 || ev== ev3){
		double ghost_etamax = 3.0;
		double ghost_area    = 0.01;
		int    active_area_repeats = 1;
		fastjet::GhostedAreaSpec ghost_spec(ghost_etamax, active_area_repeats, ghost_area);
		fastjet::AreaDefinition area_def(active_area,ghost_spec);
		fastjet::JetDefinition jet_def(kt_algorithm,R_C);
		fastjet::ClusterSequenceArea clust_seq(input_particles_combined, jet_def, area_def);
		Selector jet_selector = SelectorEtMin(EtMin_C); // && SelectorRapRange(crapmin, crapmax);
		vector<PseudoJet> inclusive_jets_combined = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));
		print_jets(inclusive_jets_combined,1);
		jet_visualization(inclusive_jets_combined, ev);
		for (unsigned int j = 0; j < inclusive_jets_combined.size(); j++){
			// cout << area_def.ghost_spec() << endl;
			cout << inclusive_jets_combined[j].area() << endl;
			cout << inclusive_jets_combined[j].description() << endl;
		}
	}

	if(inclusive_jets_combined.size()!=0){
		// Number of Jets
		ncjets.push_back(inclusive_jets_combined.size());
		ncJets->Fill(inclusive_jets_combined.size());
		ncombjets+=inclusive_jets_combined.size();

		// Jet shape
		jet_shape_incl(inclusive_jets_combined);
		jet_shape_diff(inclusive_jets_combined);

		// extract dijets
		if (inclusive_jets_combined.size() == 2){
			x_obs_dijets(inclusive_jets_combined);
			// print_jets(inclusive_jets_combined,1);
		}
	}
	

	for (int i = 0; i < inclusive_jets_combined.size(); i++){
		pxcj.push_back(inclusive_jets_combined[i].px());
		pycj.push_back(inclusive_jets_combined[i].py());
		pzcj.push_back(inclusive_jets_combined[i].pz());
		ecj.push_back(inclusive_jets_combined[i].e());
		etcj.push_back(inclusive_jets_combined[i].Et());
		phicj.push_back(inclusive_jets_combined[i].phi());
		etacj.push_back(inclusive_jets_combined[i].eta());
		phi_eta->Fill(inclusive_jets_combined[i].eta(), inclusive_jets_combined[i].phi());

		int n_constituents = inclusive_jets_combined[i].constituents().size();
		if(n_constituents > 20)
		print_jets(inclusive_jets_combined,1);

		//--------------------------------------SUBJETS for COMBINED EVENTS--------------------------------------
		dcut_C =  ycut_C*pow(inclusive_jets_combined[i].Et(),2);
		vector<PseudoJet> subJets_comb = inclusive_jets_combined[i].exclusive_subjets(dcut_C);
		
		if (subJets_comb.size() != 0){
			//nsubcjets.push_back(subJets_quark.size());
			ncsubJets->Fill(subJets_comb.size());
			ncombsubjets+=subJets_comb.size();

			double xx=inclusive_jets_combined[i].eta(), yy=subJets_comb.size();
			n_C.push_back(subJets_comb.size());
		}
	}
	
	jetcombined->Fill();
	pxcj.clear(); pycj.clear(); pzcj.clear(); 
	ecj.clear(); etcj.clear(); phicj.clear(); etacj.clear(); ncjets.clear();
	mean_js_incl.clear(); mean_js_diff.clear();
	nsubqjets.clear(); nsubgjets.clear();
}

//---------------------------------------RESOLVED JETS-----------------------------------------
void resolved_jets (const vector<fastjet::PseudoJet> & input_particles_resolved){
  
  fastjet::JetDefinition jet_def(kt_algorithm,R_R);
	fastjet::ClusterSequence clust_seq(input_particles_resolved, jet_def);
	Selector jet_selector = SelectorEtMin(EtMin_R) && SelectorRapRange(rrapmin, rrapmax); 
	vector<PseudoJet> inclusive_jets_resolved = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));
	
	if(inclusive_jets_resolved.size()!=0){
		nrjets.push_back(inclusive_jets_resolved.size());
		nrJets->Fill(inclusive_jets_resolved.size());
		nresjets+=inclusive_jets_resolved.size();
	}

	for (int i = 0; i < inclusive_jets_resolved.size(); i++){
		pxrj.push_back(inclusive_jets_resolved[i].px());
		pyrj.push_back(inclusive_jets_resolved[i].py());
		pzrj.push_back(inclusive_jets_resolved[i].pz());
		erj.push_back(inclusive_jets_resolved[i].e());
		etrj.push_back(inclusive_jets_resolved[i].Et());
		phirj.push_back(inclusive_jets_resolved[i].phi());
		etarj.push_back(inclusive_jets_resolved[i].eta());

		//--------------------------------------SUBJETS for RESOLVED EVENTS--------------------------------------
		dcut_R =  ycut_R*pow(inclusive_jets_resolved[i].Et(),2);
		vector<PseudoJet> subJets_res = inclusive_jets_resolved[i].exclusive_subjets(dcut_R);
		
		if (subJets_res.size() != 0){
			//nsubcjets.push_back(subJets_quark.size());
			nrsubJets->Fill(subJets_res.size());
			nressubjets+=subJets_res.size();
			n_R.push_back(subJets_res.size());
		}
	}
	
	jetresolved->Fill();
	pxrj.clear(); pyrj.clear(); pzrj.clear(); 
	erj.clear(); etrj.clear(); phirj.clear(); etarj.clear(); nrjets.clear();
}

//---------------------------------------DIRECT JETS-------------------------------------------
void direct_jets (const vector<fastjet::PseudoJet> & input_particles_direct){
  
  fastjet::JetDefinition jet_def(kt_algorithm,R_D);
	fastjet::ClusterSequence clust_seq(input_particles_direct, jet_def);
	Selector jet_selector = SelectorEtMin(EtMin_D) && SelectorRapRange(drapmin, drapmax); 
	vector<PseudoJet> inclusive_jets_direct = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));
	
	if(inclusive_jets_direct.size()!=0){
		ndjets.push_back(inclusive_jets_direct.size());
		ndJets->Fill(inclusive_jets_direct.size());
		ndirjets+=inclusive_jets_direct.size();
	}
	
	for (int i = 0; i < inclusive_jets_direct.size(); i++){
		pxdj.push_back(inclusive_jets_direct[i].px());
		pydj.push_back(inclusive_jets_direct[i].py());
		pzdj.push_back(inclusive_jets_direct[i].pz());
		edj.push_back(inclusive_jets_direct[i].e());
		etdj.push_back(inclusive_jets_direct[i].Et());
		phidj.push_back(inclusive_jets_direct[i].phi());
		etadj.push_back(inclusive_jets_direct[i].eta());

		//--------------------------------------SUBJETS for DIRECT EVENTS--------------------------------------
		dcut_D =  ycut_D*pow(inclusive_jets_direct[i].Et(),2);
		vector<PseudoJet> subJets_dir = inclusive_jets_direct[i].exclusive_subjets(dcut_D);
		
		if (subJets_dir.size() != 0){
			//nsubcjets.push_back(subJets_quark.size());
			ndsubJets->Fill(subJets_dir.size());
			ndirsubjets+=subJets_dir.size();
			n_D.push_back(subJets_dir.size());
		}
	}
		
	jetdirect->Fill(); gStyle->SetLineWidth(2);
	pxdj.clear(); pydj.clear(); pzdj.clear(); 
	edj.clear(); etdj.clear(); phidj.clear(); etadj.clear(); ndjets.clear();
}

//---------------------------------------QUARK QQ-SCATTERED JETS-------------------------------
void quark_jets (const vector<fastjet::PseudoJet> & input_particles_quark){
  
	qq = true;

  fastjet::JetDefinition jet_def(kt_algorithm,R_Q);
	fastjet::ClusterSequence clust_seq(input_particles_quark, jet_def);
	Selector jet_selector = SelectorEtMin(EtMin_Q) && SelectorRapRange(qrapmin, qrapmax);  
	vector<PseudoJet> inclusive_jets_quark = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));
	
	if(inclusive_jets_quark.size()!=0){
		nqjets.push_back(inclusive_jets_quark.size());
		nqJets->Fill(inclusive_jets_quark.size());
		nquarkjets+=inclusive_jets_quark.size();
	}

	for (int i = 0; i < inclusive_jets_quark.size(); i++){
		pxqj.push_back(inclusive_jets_quark[i].px());
		pyqj.push_back(inclusive_jets_quark[i].py());
		pzqj.push_back(inclusive_jets_quark[i].pz());
		eqj.push_back(inclusive_jets_quark[i].e());
		etqj.push_back(inclusive_jets_quark[i].Et());
		phiqj.push_back(inclusive_jets_quark[i].phi());
		etaqj.push_back(inclusive_jets_quark[i].eta());
		
		//--------------------------------------SUBJETS for QUARK--------------------------------------
		dcut_Q =  ycut_Q*pow(inclusive_jets_quark[i].Et(),2);
		vector<PseudoJet> subJets_quark = inclusive_jets_quark[i].exclusive_subjets(dcut_Q);
		
		if (subJets_quark.size() != 0){
			nsubqjets.push_back(subJets_quark.size());
			nqsubJets->Fill(subJets_quark.size());			
			nquarksubjets+=subJets_quark.size();
			n_Q.push_back(subJets_quark.size());
		}
		
		for (int j=0; j<subJets_quark.size(); j++){
			pxsubqj.push_back(subJets_quark[j].px());
			pysubqj.push_back(subJets_quark[j].py());
			pzsubqj.push_back(subJets_quark[j].pz());
			esubqj.push_back(subJets_quark[j].e());
			etsubqj.push_back(subJets_quark[i].Et());
			etasubqj.push_back(subJets_quark[i].eta());
			phisubqj.push_back(subJets_quark[i].phi());
		}

	}
	
	jetquark->Fill(); subjetquark->Fill();
	pxqj.clear(); pyqj.clear(); pzqj.clear(); 
	eqj.clear(); etqj.clear(); phiqj.clear(); etaqj.clear(); nqjets.clear();
	pxsubqj.clear(); pysubqj.clear(); pzsubqj.clear(); 
	esubqj.clear(); etsubqj.clear(); etasubqj.clear(); phisubqj.clear(); 
	nsubqjets.clear();
}

//---------------------------------------GLUON GG-SCATTERED JETS-------------------------------
void gluon_jets (const vector<fastjet::PseudoJet> & input_particles_gluon){
  
	gg = true;

  fastjet::JetDefinition jet_def(kt_algorithm,R_G);
	fastjet::ClusterSequence clust_seq(input_particles_gluon, jet_def);
	Selector jet_selector = SelectorEtMin(EtMin_G) && SelectorRapRange(grapmin, grapmax);  
	vector<PseudoJet> inclusive_jets_gluon = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));
	

	if(inclusive_jets_gluon.size()!=0){
		ngjets.push_back(inclusive_jets_gluon.size());
		ngJets->Fill(inclusive_jets_gluon.size());
		ngluonjets+=inclusive_jets_gluon.size();
	}

	for (int i = 0; i < inclusive_jets_gluon.size(); i++){
		pxgj.push_back(inclusive_jets_gluon[i].px());
		pygj.push_back(inclusive_jets_gluon[i].py());
		pzgj.push_back(inclusive_jets_gluon[i].pz());
		egj.push_back(inclusive_jets_gluon[i].e());
		etgj.push_back(inclusive_jets_gluon[i].Et());
		etagj.push_back(inclusive_jets_gluon[i].eta());
		phigj.push_back(inclusive_jets_gluon[i].phi());
		

		//--------------------------------------SUBJETS for GLUON--------------------------------------
		dcut_G =  ycut_G*pow(inclusive_jets_gluon[i].Et(),2);
		vector<PseudoJet> subJets_gluon = inclusive_jets_gluon[i].exclusive_subjets(dcut_G);
		
		if (subJets_gluon.size() != 0){
			nsubgjets.push_back(subJets_gluon.size());
			ngsubJets->Fill(subJets_gluon.size());
			ngluonsubjets+=subJets_gluon.size();
			n_G.push_back(subJets_gluon.size());
		}
		
		for (int j=0; j<subJets_gluon.size(); j++){
			pxsubgj.push_back(subJets_gluon[j].px());
			pysubgj.push_back(subJets_gluon[j].py());
			pzsubgj.push_back(subJets_gluon[j].pz());
			esubgj.push_back(subJets_gluon[j].e());
			etsubgj.push_back(subJets_gluon[i].Et());
			etasubgj.push_back(subJets_gluon[i].eta());
			phisubgj.push_back(subJets_gluon[i].phi());
			
		}	
	}
	jetgluon->Fill(); subjetgluon->Fill();
	pxgj.clear(); pygj.clear(); pzgj.clear(); 
	egj.clear(); etgj.clear(); phigj.clear(); etagj.clear(); ngjets.clear();
	pxsubgj.clear(); pysubgj.clear(); pzsubgj.clear(); 
	esubgj.clear(); etsubgj.clear(); etasubgj.clear(); phisubgj.clear(); 
	nsubgjets.clear();
}

//------------------------------------Jet shape (differential)---------------------------------
double jet_shape_diff (const vector<fastjet::PseudoJet> & inclusive_jets){
	
	double r=rr_diff;								// select r for R_jet = 1
	double delta_r = 0.10;			 		// R2-R1
	double R1 = r-(delta_r/2);			// inner radius
	double R2 = r+(delta_r/2);			// outer radius
	double select_const_r;					// stores the contituent radius
	double select_const_Et = 0;			// stores the contituent Et
	double select_sum = 0;					// summation
	double rho = 0;									// differential 

	// Loop over all inclusive jets
	for (int i = 0; i < inclusive_jets.size(); i++)
	{
		// store the i-jet constituents
		vector <fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();

		// loop over the constituents and select the ones which lie between R1 and R2
		for (unsigned int j = 0; j < constituents.size(); j++){

			select_const_r= pow(inclusive_jets[i].phi()-constituents[j].phi(),2) + pow(inclusive_jets[i].rap()-constituents[j].rap(),2);
			select_const_r= pow(select_const_r,0.5);

			// add Et of constituents in the doughnut of R2-R1
			if( select_const_r > R1 && select_const_r < R2){
				select_const_Et = select_const_Et + constituents[j].Et();
			}
			// divide the total selected jet contituents' Et with the total jet Et 
			select_sum = select_const_Et/inclusive_jets[i].Et();
		}
		// calculate rho, ie the jet shape for the inclusive jets
		rho = rho + (select_sum/delta_r);
		mean_js_diff.push_back(rho);
		rho_all.push_back(rho);

		// reset variable for another r run 
		select_const_Et=0;
	}
	return 0;
}

//------------------------------------Jet shape (inclusive)------------------------------------
double jet_shape_incl (const vector<fastjet::PseudoJet> & inclusive_jets){

	double select_const_r;					// stores the constituent distance from the centre of jet
	double select_const_Et=0;				// stores the constituent Et
	double select_sum;
	double psi = 0;
	double r=rr_incl;								// select r for R_jet = 1
	
	// cout<< "	r :" << rr;

	// Loop over all inclusive jets
	for (int i = 0; i < inclusive_jets.size(); i++)
	{
		// store the i-jet constituents
		vector <fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();

		// loop over the constituents and select the ones which lie < rr_inclusive 
		for (unsigned int j = 0; j < constituents.size(); j++){

			select_const_r= pow(inclusive_jets[i].phi()-constituents[j].phi(),2) + pow(inclusive_jets[i].rap()-constituents[j].rap(),2);
			select_const_r= pow(select_const_r,0.5);

			// add Et of constituents < r_inclusive
			if( select_const_r < r){
				select_const_Et = select_const_Et + constituents[j].Et();
			}
			// divide the total selected jet contituents' Et with the total jet Et 
			select_sum = select_const_Et/inclusive_jets[i].Et();
		}

		// calculate psi, ie the jet shape for the inclusive jet sample 
		psi = psi + select_sum;
		mean_js_incl.push_back(psi);

		// reset variable for another r run 
		select_const_Et=0;
	}
	//cout << " " << psi << endl; 
	return 0;
}

//------------------------------------Jet Visualisation(ev1,ev2,ev3)----------------------------------------
void jet_visualization (const vector<fastjet::PseudoJet> & jets, int ev){
	// Event 1
	if (ev==ev1) {
		for (unsigned int i = 0; i < jets.size(); i++){
			int n_constituents = jets[i].constituents().size();
			vector<PseudoJet> constituents = jets[i].constituents();

			for (unsigned int j = 0; j < n_constituents; j++){
        jetvis1->Fill( constituents[j].rap(), constituents[j].phi()); 
        c1++;
			}
			for (int k=0; k<(int)(jets[i].Et()+0.5); k++)
				{jetvis1->Fill( jets[i].rap(), jets[i].phi() ); c1++;}
			}
	}
	// Event 2
	if (ev==ev2) {
		for (unsigned int i = 0; i < jets.size(); i++){
			int n_constituents = jets[i].constituents().size();
			vector<PseudoJet> constituents = jets[i].constituents();

			for (unsigned int j = 0; j < n_constituents; j++){
				jetvis2->Fill(constituents[j].rap(), constituents[j].phi()); 
        c2++;
			}
			for (int k=0; k<(int)(jets[i].Et()+0.5); k++){
				{jetvis2->Fill( jets[i].rap(), jets[i].phi() ); c2++;}
  		}
		}
	}
	// Event 3
 	if (ev==ev3) {
		for (unsigned int i = 0; i < jets.size(); i++){
			int n_constituents = jets[i].constituents().size();
			vector<PseudoJet> constituents = jets[i].constituents();

			for (unsigned int j = 0; j < n_constituents; j++){
				jetvis3->Fill(constituents[j].rap(), constituents[j].phi()); c3++;
			}
			for (int k=0; k<(int)(jets[i].Et()+0.5); k++)
				{jetvis3->Fill( jets[i].rap(), jets[i].phi() ); c3++;}
		}
	}
}

//------------------------------------PRINT JETS and SUBJETS-----------------------------------------------
void print_jets (const vector<fastjet::PseudoJet> & jets, int flag) {
  // sort jets into increasing pt
  vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(jets);  
	cout << "Event :: " << ev << endl;
  // label the columns
  if (flag == 1)
  printf("%5s %15s %15s %15s %15s %15s\n","Jet", "rapidity", "phi", "pt", "Et", "constituents");
	if (flag == 2)
	printf("%5s %15s %15s %15s %15s %15s\n","sJet", "rapidity", "phi", "pt", "Et", "constituents");
  
  // print out the details for each jet
  for (unsigned int i = 0; i < sorted_jets.size(); i++) {
    int n_constituents = sorted_jets[i].constituents().size();
		vector<PseudoJet> constituents = sorted_jets[i].constituents();

		// print Jets (rapidity, phi_std (−π to π), pT and eT)
    printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n",
	   i, sorted_jets[i].rap(), sorted_jets[i].phi(),
	   sorted_jets[i].pt(), sorted_jets[i].Et(), n_constituents);
		cout << "		Constituents :" << endl;
		// print the constituents (rapidity, phi_std (−π to π), pT and eT)
		for (unsigned int j = 0; j < n_constituents; j++){

				cout<< "						" << constituents[j].rap() << "				"
						<< constituents[j].phi() << "				"
						<< constituents[j].pt() << "					"
						<< constituents[j].Et() << endl;
		}
		cout << "\n";
  }
}

//------------------------------------DIJETS-----------------------------------------------
double x_obs_dijets (const vector<fastjet::PseudoJet> & dijet){

	// cout << in_y << endl;
	// cout << beam_eT << endl;
	// print_jets(dijet,1);
	// double denom = 2*in_y*beam_eT;
	// double num = 0;
	// for (int i = 0; i < dijet.size(); i++){
	// 	num = num + (dijet[i].pt()*pow(exp(dijet[i].eta()),-1));
	// }
	// cout << (num/denom) << endl;
	return 0;
}

// Auxiliary: To print the average of a std::vector<double>
double average(const vector<double> &v){
    if(v.empty())	return 0;
    double avg = accumulate( v.begin(), v.end(), 0.0)/v.size();
    return avg;
}