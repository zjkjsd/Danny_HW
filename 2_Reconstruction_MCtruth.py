# +
import basf2 as b2
import modularAnalysis as ma
from variables import variables as vm
import variables.collections as vc
import variables.utils as vu
import vertex as vx

# Define the path
main_path = b2.Path()

#input_file = [''] # for grid
# the basf2 command line code is: 
# bsub -q l 'basf2 /current/directory/2_Reconstruction_MCtruth.py  D_tau/Dst_l  MC_e/MC_mu'
import sys
decaymode = 'D_l'#sys.argv[1]
light_lepton = 'e'#sys.argv[2]
input_file = f'./B2{decaymode}_N/MC/B2D_{light_lepton}_N.root'
output_file = f'./B2{decaymode}_N/Ntuples/B2D_{light_lepton}_N_truth.root'

ma.inputMdstList(environmentType='default', filelist=input_file, path=main_path)


ma.fillParticleListFromMC('pi+:MC', cut='', skipNonPrimaryDaughters=True,path=main_path)
ma.fillParticleListFromMC('K-:MC', cut='', skipNonPrimaryDaughters=True,path=main_path)
ma.fillParticleListFromMC('e-:MC', cut='',skipNonPrimaryDaughters=True,path=main_path)
ma.fillParticleListFromMC("gamma:MC", '',skipNonPrimaryDaughters=True, path=main_path)
ma.fillParticleListFromMC('nu_e:MC', '',skipNonPrimaryDaughters=True, path=main_path)
ma.fillParticleListFromMC('mu-:MC', cut='',skipNonPrimaryDaughters=True,path=main_path)
#ma.fillParticleListFromMC('N_e:MC', '',skipNonPrimaryDaughters=True, path=main_path)

# Event Kinematics
ma.buildEventKinematics(fillWithMostLikely=True,path=main_path)
ma.buildEventKinematicsFromMC(path=main_path)



# Reconstruct D
Dcuts = '1.8 < M < 1.9'
ma.reconstructMCDecay('D+:MC -> K-:MC pi+:MC pi+:MC', cut='', path=main_path)
#ma.reconstructMCDecay('N_e:MC -> e-:MC e+:MC nu_e:MC', cut='',path=main_path)

# Reconstruct B
ma.reconstructMCDecay('B0:MC =direct=> D-:MC e+:MC e-:MC e+:MC nu_e:MC', cut='', path=main_path)

# Calculate the distance between vertices De and D+
vm.addAlias('vtxDD', 'vertexDistanceOfDaughter(0, D+)')
vm.addAlias('vtxDDSig', 'vertexDistanceOfDaughterSignificance(0, D+)')

# MC Truth Matching
ma.matchMCTruth('B0:MC', path=main_path)

# generate the decay string
#main_path.add_module('ParticleMCDecayString', listName='anti-B0:MC', fileName=output_hash)
#vm.addAlias('DecayHash','extraInfo(DecayHash)')
#vm.addAlias('DecayHashEx','extraInfo(DecayHashExtended)')





# build the rest of the event
ma.buildRestOfEventFromMC('B0:MC',path=main_path)
roe_mask1 = ('my_mask', '','')
ma.appendROEMasks('B0:MC', [roe_mask1], path=main_path)

# ROE variables
roe_kinematics = ["roeE(my_mask)", "roeP(my_mask)", "roePx(my_mask)",
                  "roePy(my_mask)","roePz(my_mask)","roePt(my_mask)",]
roe_E_Q = ['roeCharge(my_mask)', 'roeNeextra(my_mask)','roeEextra(my_mask)',]

roe_multiplicities = ["nROE_Charged(my_mask)",'nROE_ECLClusters(my_mask)',
                      'nROE_NeutralECLClusters(my_mask)','nROE_KLMClusters',
                      'nROE_NeutralHadrons(my_mask)',"nROE_Photons(my_mask)",
                      'nROE_Tracks(my_mask)','nROE_RemainingTracks(my_mask)',]

vm.addAlias('nROE_e','nROE_Charged(my_mask, 11)')
vm.addAlias('nROE_mu','nROE_Charged(my_mask, 13)')
vm.addAlias('nROE_K','nROE_Charged(my_mask, 321)')
vm.addAlias('nROE_pi','nROE_Charged(my_mask, 211)')
roe_nCharged = ['nROE_e','nROE_mu','nROE_K','nROE_pi',]

# Load missing momentum in the event and use a mask 'cleanMask':
ma.fillParticleListFromROE('nu_e:missing', '', maskName='my_mask',
                           sourceParticleListName='B0:MC', useMissing = True, path=main_path)

# fit B vertex on the tag-side
vx.TagV("B0:MC", fitAlgorithm="Rave", maskName='my_mask', path=main_path)

# Continuum Suppression
ma.buildContinuumSuppression(list_name="B0:MC", roe_mask="my_mask", path=main_path)

CSVariables = ["R2","thrustBm","thrustOm","cosTBTO","cosTBz","isContinuumEvent",]






# Write variables to Ntuples
vm.addAlias('cos_pV','cosAngleBetweenMomentumAndVertexVector')
vm.addAlias('cos_pB','cosThetaBetweenParticleAndNominalB')

# Kinematic variables in CMS
cms_kinematics = vu.create_aliases(vc.kinematics, "useCMSFrame({variable})", "CMS")
roe_cms_kinematics = vu.create_aliases(roe_kinematics, "useCMSFrame({variable})", "CMS")


b_vars = vu.create_aliases_for_selected(
    list_of_variables= cms_kinematics + vc.deltae_mbc + vc.inv_mass +vc.mc_truth
    + ['dM','Q','dQ','missingMass2OfEvent','m2Recoil','cos_pV','cos_pB', 
       "roeMbc(my_mask)", "roeM(my_mask)","roeDeltae(my_mask)",]
    + roe_cms_kinematics + roe_E_Q + roe_multiplicities + roe_nCharged + CSVariables
    + vc.tag_vertex + vc.mc_tag_vertex, # + ft.flavor_tagging
    decay_string='^B0:MC -> D-:MC e+:MC e-:MC e+:MC nu_e:MC',
    prefix=['B0'])

D_vars = vu.create_aliases_for_selected(
    list_of_variables= cms_kinematics + vc.kinematics + vc.inv_mass + vc.mc_truth
    + ['dM',],
    decay_string='B0:MC -> ^D-:MC ^e+:MC ^e-:MC ^e+:MC ^nu_e:MC',
    prefix=['D','e_0','e_1','e_2','nu_e'])

#e_vars = vu.create_aliases_for_selected(
#    list_of_variables= cms_kinematics + vc.kinematics + vc.inv_mass + vc.mc_truth,
#    decay_string='N_e:MC -> ^e-:MC ^e+:MC nu_e:MC',
#    prefix=['e_1','e_2'])

nu_vars = vu.create_aliases_for_selected(
    list_of_variables= cms_kinematics + vc.inv_mass,
    decay_string='^nu_e:missing',
    prefix=['nu_miss'])

candidate_vars = b_vars + D_vars + nu_vars
event_vars=['Ecms', 'IPX', 'IPY', 'IPZ'] + vc.event_kinematics + vc.mc_event_kinematics
ma.variablesToNtuple('B0:MC', event_vars + candidate_vars,# + mc_gen_topo(30),
                     filename=output_file, treename='B0', path=main_path)
#ma.variablesToNtuple('', event_vars, filename=output_file, treename='event', path=main_path)

b2.process(path=main_path)
print(b2.statistics)

# -


