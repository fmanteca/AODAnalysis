import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras


process = cms.Process("Demo", eras.Phase2)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.StandardSequences.IdealGeometry_cff')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")



#process.load('Configuration.StandardSequences.Services_cff')
#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')



# import of standard configurations
process.load("MyAnalysis.AODAnalysis.IFCAAODAnalysis_cff")




#process.load("Configuration.Geometry.GeometrySimTracker_cff")
process.GlobalTag.globaltag = '94X_mc2017_realistic_v10'  # or some other global tag depending on your CMSSW release and sample. 
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/gpfs/users/mantecap/CMSSW_10_6_0_patch2/src/MyAnalysis/AODAnalysis/F0A58BFD-00FF-E443-95A9-FB04741A0179.root'
#        'file:/gpfs/users/mantecap/CMSSW_10_6_0_patch2/src/MyAnalysis/AODAnalysis/SingleMu-Run2018D-ZMu-PromptReco-v2.root'
#        'file:/gpfs/users/mantecap/CMSSW_10_6_0_patch2/src/MyAnalysis/AODAnalysis/GEN-SIM-RECO_RelValZprimeToll_M3000.root'
        'file:/gpfs/projects/cms/fernanpe/MUO-RunIIFall17DRPremix-00025.root'
#        'file:/gpfs/projects/cms/fernanpe/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/ZprimeToMuMu_M-6000_NoPU_106X_upgrade2023_RECO/200207_114013/0000/RECO_Prod_RAW2DIGI_L1Reco_RECO_RECOSIM_inRECOSIM_40.root'
#        'file:/gpfs/projects/cms/fernance/H2ToLLPXToLeptons_MH_400_MX_50_ctau_4mm_TuneCP2_13TeV_pythia8_80X_13082019-1313/400-50-4_HXX_RunIISummer16DR80Premix_step2_230220-1650/200225_115305/0000/EXO-RunIISummer16DR80Premix-03591_308.root'
#        'file:/gpfs/users/mantecap/CMSSW_10_6_0_patch2/src/MyAnalysis/AODAnalysis/FC006F2D-2B75-224B-A5E6-D2D6553A79E7.root'
#        'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring19DR/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/AODSIM/NoPU_106X_upgrade2023_realistic_v3-v2/50000/F0A58BFD-00FF-E443-95A9-FB04741A0179.root'
    )
)





process.p = cms.Path(process.AODanalyzer)


