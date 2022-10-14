//P1
#include "fastjet/ClusterSequence.hh"
#include "TH1.h"
#include <TF2.h>
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <TGraph2D.h>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "Riostream.h"
#include <TLorentzVector.h>
#include <iostream> 
#include <sstream>  // needed for internal io
#include <vector> 
#include <cstdio>

using namespace fastjet;
using namespace std;
using namespace ROOT;

void combined_jets (const vector<fastjet::PseudoJet> &);
double EtMin_C=0.0, R_C=0.01;
void direct_jets (const vector<fastjet::PseudoJet> &);
double EtMin_D=0.0, R_D=0.01;
void resolved_jets (const vector<fastjet::PseudoJet> &);
double EtMin_R=0.0, R_R=0.01;

void print_jets (const vector<fastjet::PseudoJet> &, int flag);


//Reading  the input root file.
TFile *myFile = new TFile("data.root","READ");
//For Combined Events :
TTreeReader myReader("combinedevents", myFile);
TTreeReaderArray<Float_t> px(myReader, "px");
TTreeReaderArray<Float_t> py(myReader, "py");
TTreeReaderArray<Float_t> pz(myReader, "pz");
TTreeReaderArray<Float_t> e(myReader, "e");
//Reading Beam Electrons :
TTreeReaderArray<Float_t> epx(myReader, "bpx");
TTreeReaderArray<Float_t> epy(myReader, "bpy");
TTreeReaderArray<Float_t> epz(myReader, "bpz");
TTreeReaderArray<Float_t> ee(myReader, "be");
//Reading Beam Photons (e->gamma) :
TTreeReaderArray<Float_t> gmpx(myReader, "gmpx");
TTreeReaderArray<Float_t> gmpy(myReader, "gmpy");
TTreeReaderArray<Float_t> gmpz(myReader, "gmpz");
TTreeReaderArray<Float_t> gme(myReader, "gme");
//Reading flags for direct or resolved :
TTreeReaderArray<Int_t> dirflag(myReader, "dirflag");
TTreeReaderArray<Int_t> resflag(myReader, "resflag");
TTreeReaderArray<Int_t> events(myReader, "events");

//Writing an Output ROOT file:
TFile *jets = new TFile("jets.root","RECREATE");

std::vector<Float_t> pxcj,pycj,pzcj,ecj,etcj,phicj,etacj,ncjets;
TTree *jetcombined = new TTree("jetcombined", "");

std::vector<Float_t> pxrj,pyrj,pzrj,erj,etrj,phirj,etarj,nrjets;
TTree *jetresolved = new TTree("jetresolved", "");

std::vector<Float_t> pxdj,pydj,pzdj,edj,etdj,phidj,etadj,ndjets;
TTree *jetdirect = new TTree("jetdirect", "");

int ncombjets=0, ndirjets=0, nresjets=0;

int main(){

	vector<fastjet::PseudoJet> input_particles_combined;
	vector<fastjet::PseudoJet> input_particles_direct;
	vector<fastjet::PseudoJet> input_particles_resolved;
	
	int ev=0;	
	int c=0, d=0, r=0, cs=0, ds=0, rs=0;
	
	//COMBINED ENTRIES:
	jetcombined->Branch("combined_px", "std::vector<Float_t>",  &pxcj);
	jetcombined->Branch("combined_py", "std::vector<Float_t>",  &pycj);
	jetcombined->Branch("combined_pz", "std::vector<Float_t>",  &pzcj);
	jetcombined->Branch("combined_E", "std::vector<Float_t>",  &ecj);
	jetcombined->Branch("combined_Et", "std::vector<Float_t>",  &etcj);
	jetcombined->Branch("combined_phi", "std::vector<Float_t>",  &phicj);
	jetcombined->Branch("combined_eta", "std::vector<Float_t>",  &etacj);
	jetcombined->Branch("combined_njets", "std::vector<Float_t>",  &ncjets);
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
	
	while (myReader.Next()){
		
		for (int i = 0; i < events.GetSize(); i++){
			ev++;
			
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
			
			combined_jets(input_particles_combined);
			resolved_jets(input_particles_resolved);
			direct_jets(input_particles_direct);
			
			input_particles_combined.clear();
			input_particles_resolved.clear();
			input_particles_direct.clear();
		}
  }
  // a "header" for the output
  //cout << "Ran " << jet_def.description() << endl;
	//cout << "Rapidity Range: -1 to 2" << endl;
	//cout << "Showing the subjets with ycut: " << ycut << endl << endl;
	
  cout << "\nNumber of Combined Events Read : " << ev << " (Constituents : " << cs << " )" << endl;
  cout << "Number of Resolved Events : " << r << " (Constituents : " << rs << " )" << endl;
  cout << "Number of Direct Events : " << d << " (Constituents : " << ds << " )" << endl;
  
  cout << "\n#Jets formed (with Combined Jet Constituents : " << cs << ") : " << 	ncombjets << endl;
  cout << "		Minimum Etmin Combined: " << EtMin_C << " GeV" << endl;
  cout << "		R-Combined: " << R_C << endl;
  cout << "\n#Jets formed (with Resolved Jet Constituents : " << rs << ") : " << 	nresjets << endl;
  cout << "		Minimum Etmin Resolved: " << EtMin_R << " GeV" << endl;
  cout << "		R-Resolved: " << R_R << endl;
  cout << "\n#Jets formed (with Direct Jet Constituents : " << ds << ") : " << 	ndirjets << endl;
  cout << "		Minimum Etmin Direct: " << EtMin_D << " GeV" << endl;
  cout << "		R-Direct: " << R_D << endl;
  
  jets->Write();
	jets->Close();
	delete jets;
}



//------------------------------------COMBINED JETS------------------------------------
void combined_jets (const vector<fastjet::PseudoJet> & input_particles_combined){
	
	fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recomb_scheme=fastjet::E_scheme;
  fastjet::JetAlgorithm jet_alg = fastjet::kt_algorithm;
  fastjet::JetDefinition jet_def(jet_alg, R_C, recomb_scheme, strategy);
  
	fastjet::ClusterSequence clust_seq(input_particles_combined, jet_def);
	Selector jet_selector = SelectorEtMin(EtMin_C);
	vector<PseudoJet> inclusive_jets_combined = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));
	
	ncombjets+=inclusive_jets_combined.size();
	
	for (int i = 0; i < inclusive_jets_combined.size(); i++){
		pxcj.push_back(inclusive_jets_combined[i].px());
		pycj.push_back(inclusive_jets_combined[i].py());
		pzcj.push_back(inclusive_jets_combined[i].pz());
		ecj.push_back(inclusive_jets_combined[i].e());
		etcj.push_back(inclusive_jets_combined[i].Et());
		phicj.push_back(inclusive_jets_combined[i].phi());
		etacj.push_back(inclusive_jets_combined[i].eta());
	}
	ncjets.push_back(inclusive_jets_combined.size());
	
	jetcombined->Fill();
	pxcj.clear(); pycj.clear(); pzcj.clear(); 
	ecj.clear(); etcj.clear(); phicj.clear(); etacj.clear(); ncjets.clear();
}



//------------------------------------RESOLVED JETS------------------------------------
void resolved_jets (const vector<fastjet::PseudoJet> & input_particles_resolved){
	//double EtMin=0.0, R=0.1;
	
	fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recomb_scheme=fastjet::E_scheme;
  fastjet::JetAlgorithm jet_alg = fastjet::kt_algorithm;
  fastjet::JetDefinition jet_def(jet_alg, R_R, recomb_scheme, strategy);
  
	fastjet::ClusterSequence clust_seq(input_particles_resolved, jet_def);
	Selector jet_selector = SelectorEtMin(EtMin_R);
	vector<PseudoJet> inclusive_jets_resolved = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));
	
	nresjets+=inclusive_jets_resolved.size();
	for (int i = 0; i < inclusive_jets_resolved.size(); i++){
		pxrj.push_back(inclusive_jets_resolved[i].px());
		pyrj.push_back(inclusive_jets_resolved[i].py());
		pzrj.push_back(inclusive_jets_resolved[i].pz());
		erj.push_back(inclusive_jets_resolved[i].e());
		etrj.push_back(inclusive_jets_resolved[i].Et());
		phirj.push_back(inclusive_jets_resolved[i].phi());
		etarj.push_back(inclusive_jets_resolved[i].eta());
	}
	nrjets.push_back(inclusive_jets_resolved.size());
	
	jetresolved->Fill();
	pxrj.clear(); pyrj.clear(); pzrj.clear(); 
	erj.clear(); etrj.clear(); phirj.clear(); etarj.clear(); nrjets.clear();
}



//------------------------------------DIRECT JETS------------------------------------
void direct_jets (const vector<fastjet::PseudoJet> & input_particles_direct){
	//double EtMin=0.0, R=0.1;
	
	fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recomb_scheme=fastjet::E_scheme;
  fastjet::JetAlgorithm jet_alg = fastjet::kt_algorithm;
  fastjet::JetDefinition jet_def(jet_alg, R_D, recomb_scheme, strategy);
  
	fastjet::ClusterSequence clust_seq(input_particles_direct, jet_def);
	Selector jet_selector = SelectorEtMin(EtMin_R);
	vector<PseudoJet> inclusive_jets_direct = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));
	
	ndirjets+=inclusive_jets_direct.size();
	
	for (int i = 0; i < inclusive_jets_direct.size(); i++){
		pxdj.push_back(inclusive_jets_direct[i].px());
		pydj.push_back(inclusive_jets_direct[i].py());
		pzdj.push_back(inclusive_jets_direct[i].pz());
		edj.push_back(inclusive_jets_direct[i].e());
		etdj.push_back(inclusive_jets_direct[i].Et());
		phidj.push_back(inclusive_jets_direct[i].phi());
		etadj.push_back(inclusive_jets_direct[i].eta());
	}
	ndjets.push_back(inclusive_jets_direct.size());
	
	jetdirect->Fill();
	pxdj.clear(); pydj.clear(); pzdj.clear(); 
	edj.clear(); etdj.clear(); phidj.clear(); etadj.clear(); ndjets.clear();
}




void print_jets (const vector<fastjet::PseudoJet> & jets, int flag) {

  // sort jets into increasing pt
  vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(jets);  

  // label the columns
  if (flag == 1)
  printf("%5s %15s %15s %15s %15s\n","Jet", "rapidity", "phi", "pt", "constituents");
	if (flag == 2)
	printf("%5s %15s %15s %15s %15s\n","sJet", "rapidity", "phi", "pt", "constituents");
  
  // print out the details for each jet
  for (unsigned int i = 0; i < sorted_jets.size(); i++) {
    int n_constituents = sorted_jets[i].constituents().size();
    printf("%5u %15.8f %15.8f %15.8f %8u\n",
	   i, sorted_jets[i].rap(), sorted_jets[i].phi(),
	   sorted_jets[i].perp(), n_constituents);
  }
}




