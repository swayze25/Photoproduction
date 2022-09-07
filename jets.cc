#include "fastjet/ClusterSequence.hh"
#include "TH1.h"
#include <TF2.h>
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <TGraph2D.h>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "Riostream.h"
#include <TLorentzVector.h>
#include <iostream> // needed for io

using namespace fastjet;
using namespace std;

int main(){
	
	//TCanvas *c = new TCanvas("c","-",900,900);
	TH2F *h2 = new TH2F("h2","Jet Constituent",100,-50,50,100,0,100);
	TH2F *jvis = new TH2F("jvis","Jet Visualization in the Rapidity-Azimuth plane",100,-50,50,100,0,100);
	//TGraph2D *dt = new TGraph2D();
	//dt->SetTitle("Heat Map of pT of the Generated Events; X axis : Rapidity; Y axis : Azimuth Angle; Z axis : Transverse Momentum");
	
	//Create a PseudoJet
  vector<fastjet::PseudoJet> input_particles;
  
  //Open root file to read the event entries
  TFile *myFile = new TFile("data.root","READ");
  TTreeReader myReader("events", myFile);
  TTreeReaderArray<Float_t> px(myReader, "px");
  TTreeReaderArray<Float_t> py(myReader, "py");
  TTreeReaderArray<Float_t> pz(myReader, "pz");
  TTreeReaderArray<Float_t> e(myReader, "e");
  
  //Create a root file to save the Jet Parameters
  TFile *jetreco = new TFile("jetsreco.root","RECREATE");
  TTree *jetparam = new TTree("jetparam", "");
  std::vector<Float_t> pt,m,ej;
  jetparam->Branch("energy", "std::vector<Float_t>",  &ej);
  jetparam->Branch("mass", "std::vector<Float_t>",  &m);
  jetparam->Branch("pt", "std::vector<Float_t>",  &pt);
  
  //Fill the PseudoJet with 4 momentum  
  while (myReader.Next()){
  	for (int i = 0; i < px.GetSize(); i++){
  		input_particles.push_back(fastjet::PseudoJet(px[i],py[i],pz[i],e[i]));
  	}
  }
  
  // choose a jet definition
	double R = 1.0;
	JetDefinition jet_def(kt_algorithm, R);
	
	// run the clustering, extract the jets
	ClusterSequence cs(input_particles, jet_def);
	vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
	
	// print out some info
	cout << "Clustered with " << jet_def.description() << endl;
	
	int count = 1;
	
	// print the jets
	cout << " pt y phi" << endl;
	for (int i = 0; i < jets.size(); i++) {
		
		cout << "jet " << i << ": "<< jets[i].pt() << " "	<< jets[i].rap() << " " << jets[i].phi() << endl;
			ej.push_back(jets[i].e());
			pt.push_back(jets[i].pt());
			m.push_back(jets[i].m());
			jetparam->Fill();
			ej.clear(); pt.clear(); m.clear(); 
			h2->Fill(jets[i].rap(), jets[i].phi());
		vector<PseudoJet> constituents = jets[i].constituents();
		//if (count == 1)
		for (int j = 0; j < constituents.size(); j++) {
			//cout << " constituent " << j << "’s pt: "<< constituents[j].pt() << endl;
			//h2->Fill(constituents[j].rap(),constituents[j].phi());
			//dt->SetPoint(j,constituents[j].rap(),constituents[j].phi(),constituents[j].pt());
		}
		jvis->Add(h2);
		h2->Reset();
		//h2->Write();
		//count=2;
	}
	
	jvis->Write();
	jetreco->Write();
	jetreco->Close();
	delete jetreco;
	

	//gStyle->SetPalette(1);
  //dt->Draw("surf1");	
  //c->Show();
  //return c;
}
	
	
	
	
	
