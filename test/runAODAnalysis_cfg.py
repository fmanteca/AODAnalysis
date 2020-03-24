import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalysis.AODAnalysis.IFCAAODAnalysis_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Choose Tracker Geometry
# process.load("Configuration.Geometry.GeometryReco_cff")
# process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
# process.load("Geometry.TrackerGeometryBuilder.idealForDigiTrackerGeometry_cff")
# process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
# process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
# process.load("Geometry.TrackerNumberingBuilder.trackerTopologyConstants_cfi")

process.load("Configuration.Geometry.GeometrySimTracker_cff")
process.GlobalTag.globaltag = '106X_upgrade2023_realistic_v2'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/gpfs/users/mantecap/CMSSW_10_6_0_patch2/src/MyAnalysis/AODAnalysis/F0A58BFD-00FF-E443-95A9-FB04741A0179.root'
#        'file:/gpfs/users/mantecap/CMSSW_10_6_0_patch2/src/MyAnalysis/AODAnalysis/SingleMu-Run2018D-ZMu-PromptReco-v2.root'
#        'file:/gpfs/users/mantecap/CMSSW_10_6_0_patch2/src/MyAnalysis/AODAnalysis/GEN-SIM-RECO_RelValZprimeToll_M3000.root'
        'file:/gpfs/projects/cms/fernanpe/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/ZprimeToMuMu_M-6000_NoPU_106X_upgrade2023_RECO/200207_114013/0000/RECO_Prod_RAW2DIGI_L1Reco_RECO_RECOSIM_inRECOSIM_40.root'
#        'file:/gpfs/users/mantecap/CMSSW_10_6_0_patch2/src/MyAnalysis/AODAnalysis/FC006F2D-2B75-224B-A5E6-D2D6553A79E7.root'
#        'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring19DR/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/AODSIM/NoPU_106X_upgrade2023_realistic_v3-v2/50000/F0A58BFD-00FF-E443-95A9-FB04741A0179.root'
    )
)





process.p = cms.Path(process.AODanalyzer)


