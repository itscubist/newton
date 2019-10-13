# newton
Numeric Event generator in Water for Tens Of MeV Neutrinos
Tell me when you have a better name :)

CURRENT CROSS SECTIONS

All Cross sections are stored in numeric format in xscnData directory. Cross-section calculations for ibd and elastic scattering events are known and replicated from the source material to get numerical data.

The main reason of usage of numerical data is due to 16O cross sections being obtained from digitzed graphs/tables from Haxton and Nakazato papers.

If I had full matrix elements for all the cross sections, then adding these to MARLEY would made more sense, but for now I decided to build this.

IBD (0-200 MeV)  (Strumia - Vissani 2003)
nue-16O (total xscn good between 0-200, but angular and excitation energy dists. are good up to 120 MeV) (Haxton 1987)
nuebar-16O (total xscn good between 0-200, but angular and excitation energy dists. are good up to 120 MeV) (Haxton 1987)
elastic scattering for nue,nuebar,nux,nuxbar (0-200 MeV) (Bahcall et al. 1995)

NC on 16O is not in yet but will implement possibly from Langanke, or Kolbe's papers in the future

INSTALL

One needs to install TALYS separately at https://tendl.web.psi.ch/tendl_2019/talys.
And set $TALYS variable to installed TALYS location (where source, structure, samples directories can be seen)

Then use makefile: make newtonTest


USAGE: 

Currenty not very user friendly...

Control over different xscns, fluxes and detectors are done with card files

Make a detector file for your detector (newton supporst rectangular prism, spherical and cylindrical detectors)
You can copy it from detData/cardSuperK.txt
The material will determine which interactions are included, with what weight, also material's properties will be asked
The years are important for number of events generated...
Probably you do not want to change Avogadro's number and other constants...
You can set detector dimensions at the last part:
  if cylindrical: then set z and r to non-zero values
  if rectangular prism: then set x,y,z to non-zero values
  if spherical: then set r to non-zero values
  others will be ignored.


Input a flux, currently only read from a ROOT file. You can put different flux for 6 type of neutrinos in
TH3D format with axis being (X: zenith from detector top, Y: Azimuth angle, Z:energy, and binContent is the flux)
The number of bins should not matter as long as correct bin numbers are presented in the flux file
Alternative is to provide TH1D histograms with X axis being energy
In that case, you can choose either isotropic or single direction flux
You should also set "scale" and "escale" to convert units of flux and energy axis to ones used in the detector (see available card files for better description)

If you want to use an included cross section in a new material:
  in cross section file add a newline following the convention with WATER X, X being the number of targets per molecule fo water (10 for electrons, 2 for protons,1 for oxygen nuclei)


Run ./newtonTest 
Will output newtonTestRoot.out with some diagnostic graphs/histograms
TALYS decay results will be in talysInterface/talysOutputs/*.root
Number of events generated will be printed out
To decided on energy range interested you have to manually change (will fix that soon)
Also outputs newtonTest.kin in NUANCE format (to be used in simulators)
