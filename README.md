# Photoproduction
To Study Jets via Photoproduction 

Main program to generate charged hadron spectra from photon-initiated
hard processes, by combining sub-runs with direct or resolved photons
or by generating all with contributions in a single run (automaticmix)

--------------------------------------------------------------------------------------------------------------

Further, I have set the flags for event production: 
bool photonProton = true;
bool photonsFromElectrons = true;
bool automaticMix = true;

I am only considering the beam particles for (electron -> photon) + proton, for:
Resolved - Resolved, and Direct - Resolved contributions only.

## BEAMS and ENERGIES:
  - Beam1 : 11 : electron ("PDF:beamA2gamma = on" ie photons are produced by electrons)
  - Beam2 : 2212 : proton
  - Beams:eCM = 140. : com energy
  
HERA Energies: CoM Energy ranges from 134–277 GeV at HERA 
1. pythia.readString("Beams:eCM = 140.");	 
2. pythia.readString("Beams:eA  = 27.5.");
3. pythia.readString("Beams:eB  = 820.");

EIC Energies: CoM Energy ranges from 30-140 GeV at HERA 
1. pythia.readString("Beams:eCM = 140.");	 
2. pythia.readString("Beams:eA  = 18.");
3. pythia.readString("Beams:eB  = 275.");

### CUTS: (and Initializations on Event Generation)
1. Virtuality : Photon:Q2max = 1.0                 (pythia.readString("Photon:Q2max = 1.0");)
2. Invariant Mass : Photon:Wmin 
3. Limit partonic phat : PhaseSpace:pTHatMin", 5.0
4. MultipartonInteractions:pT0Ref = 3.00
5. HardQCD:all = on
6. PhotonParton:all = on

--------------------------------------------------------------------------------------------------------------

## JETS:

The macro file 'jets_photoproduction.cc' is used for the formation of jets.
Uses Fastjet+ROOT libraries.
