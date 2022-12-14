#!/usr/bin/env python3
#####################################################################
# The Belle II code is designed to load the proper geometry of the detector depending on which experiment and run combination is used. Experiment 1003 is reserved as a place holder for run-independent (run=0) MC for early Phase 3 (partial PXD installed). Experiment 0 and run 0 are reserved for run-independent MC for Phase 3 (full PXD installed). (https://confluence.desy.de/display/BI/Experiment+numbering)
#
# We will also start by producing 10000 events. You can change this later.
#
# Note that the global tag used in the script will be the default one for the release, which is appropriate for run-independent MC. If you produce run-dependent MC, you will need to set the appropriate global tag (https://confluence.desy.de/display/BI/Global+Tag+%28GT%29+page).
#
# ####################################################################

# the command line code is: 
# bsub -q l 'basf2 /current/directory/0_Gen_MC.py  D_l  e/mu  1  0.001  -n 10000'

import sys
decaymode = sys.argv[1]
light_lepton = sys.argv[2]

import pdg
M = sys.argv[3]
W = sys.argv[4]
pdg.add_particle(name='N_e', pdgCode=777, mass=float(M), width=float(W), charge=0, spin=0.5, lifetime=0)
pdg.add_particle(name='anti-N_e', pdgCode=-777, mass=float(M), width=float(W), charge=0, spin=0.5, lifetime=0)
#pdg.add_particle(name='N_mu', pdgCode=778, mass=0.005, width=999, charge=0, spin=0.5, lifetime=0)
#pdg.add_particle(name='anti-N_mu', pdgCode=-778, mass=0.005, width=999, charge=0, spin=0.5, lifetime=0)

decfile=f'/home/belle/zhangboy/Danny_HW/B2{decaymode}_N/decfiles/B2D_{light_lepton}_N.dec'
output =f'/home/belle/zhangboy/Danny_HW/B2{decaymode}_N/MC/N_{light_lepton}_mass{M}_width{W}_2.root'

import basf2 as b2
mypath=b2.Path()

# Load the EventInfoSetter module and set the exp/run/evt details
# expList=1003 for early phase 3, 0 for full Belle2 geometry
mypath.add_module("EventInfoSetter", expList=0, runList=0, evtNumList=100)

# Add the generator
import generators as ge
ge.add_evtgen_generator(path=mypath, finalstate='signal', signaldecfile=decfile)

# Simulate the detector response and the L1 trigger
import simulation as si
import background
si.add_simulation(path=mypath, 
                  bkgfiles=background.get_background_files())

# Simulate the L1 trigger
#import L1trigger as l1
#l1.add_tsim(path=mypath)

# Reconstruct the objects
import reconstruction as re
re.add_reconstruction(path=mypath)

# Create the mDST output file
import mdst
mdst.add_mdst_output(path=mypath, filename=output)

# Process the steering path
b2.process(path=mypath)

# Finally, print out some statistics about the modules execution
print(b2.statistics)

###############################################################################
# ### Birks coefficients used in run time
#    G4_POLYSTYRENE     0.07943 mm/MeV     0.00841958 g/cm^2/MeV  massFactor=  101.167 effCharge= 0.027027
#  FLOATING-POINT NUMBERS ASSUMED ACCURATE TO 1e-06
#  FLOATING-POINT NUMBERS ASSUMED ACCURATE TO 1e-06
# EvtGen:In readDecayFile, reading:./B2D_l_N/decfiles/B2D_e_N.dec
# EvtGen:As requested, PHOTOS will be turned on for all decays.
# EvtGen:Redefined decay of Upsilon(4S)
# Myanti-B0 -> MyD+ e- anti-N_e  (BGL)
# EvtGen:BGL did not get the correct daughter spin d=2
# EvtGen:Will terminate execution!
# abort() called, exiting



