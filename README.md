# Photoproduction
To Study Jets via Photoproduction 

Main program to generate charged hadron spectra from photon-initiated
hard processes, by combining sub-runs with direct or resolved photons
or by generating all with contributions in a single run (automaticmix)

I have modified the main69.ccn(PYTHIA8 example folder) and have made the following changes:
1. The event generator file is labelled as 'eventgen.cc' 
  a. The macro has to be made (requires a makefile) so it is system-specific for now (cant run on other system)
  b. It requires ROOT+Pythia libraries, the events are stored in a log file 'outevents',
     and the events are stored as 4mom in the ROOT file 'data.root'
2. Open the ROOT file 'data.root' to see the 4mom of charged final states hadrons stored in respective branches.

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

### CUTS: (and Initializations on Event Generation)
1. Virtuality : Photon:Q2max = 1.0                 (pythia.readString("Photon:Q2max = 1.0");)
2. Invariant Mass : Photon:Wmin  = 134.0           (pythia.readString("Photon:Wmin  = 134.0");)
3. Limit partonic phat : PhaseSpace:pTHatMin", 5.0
4. MultipartonInteractions:pT0Ref = 3.00
5. HardQCD:all = on
6. PhotonParton:all = on

## HISTOGRAMS:
pTtotR : for the pT values of all the Charged and Final state particles
momentumx, momentumy, momentumz, energy,  mass : 4-momentum and mass of the Charged and FInal state particles
pTresresR, pTresdirR, pTdirresR, pTdirdirR : Various Contributions via Photoproduction interactions

## BRANCHES:
1. combinedevents : stores the 4 momentumm of the Charged and Final state particles (px.py,pz and e)
2. qqbar2qqbar : stores the 4 momentumm of the Charged and Final state particles for quark initiated processes
3. gg2gg : stores the 4 momentumm of the Charged and Final state particles for gluon initiated processes

--------------------------------------------------------------------------------------------------------------

## JETS:

Contains a root file jets.root which stores the jets formed by a standalone FastJet program using kT algorithm.
Read the jets.root for the output


