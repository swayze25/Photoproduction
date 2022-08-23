# Photoproduction
To Study Jets via Photoproduction 

Main program to generate charged hadron spectra from photon-initiated
hard processes, by combining sub-runs with direct or resolved photons
or by generating all with contributions in a single run.

I have modified the main69.cc and have made the following changes:
1. The output Histograms are now stored in a root file named photoproduction.root* instead of being displayed on the output log file.
2. main69.cc have been renamed as photoproduction.cc (Requires a change in makefile)

*(ROOT libraries have been enabled and necessary modifications, check line have been done in the makefile)

--------------------------------------------------------------------------------------------------------------

Further, I have changed the flags : (No Automatic Mixing)
bool photonProton = true;
bool photonsFromElectrons = true;
bool automaticMix = false;

I am only considering the beam particles for (electron -> photon) + proton, for:
Resolved - Resolved, and Direct - Resolved contributions only.

BEAMS: 
Beam1 : 11 : electron ("PDF:beamA2gamma = on" ie photons are produced by electrons)
Beam1 : 2212 : proton

CUTS: (and Initializations)
1. Virtuality : Photon:Q2max = 1.0
2. Invariant Mass : Photon:Wmin  = 10.0 
3. Limit partonic phat : PhaseSpace:pTHatMin", 5.0
4. MultipartonInteractions:pT0Ref = 3.00
5. HardQCD:all = on
6. PhotonParton:all = on

HISTOGRAMS:
pTtotR : for the pT values of all the Charged and FInal state particles
momentumx, momentumy, momentumz, energy,  mass : 4-momentum and mass of the Charged and FInal state particles
pTresresR, pTresdirR, pTdirresR, pTdirdirR : Various Contributions via Photoproduction interactions


