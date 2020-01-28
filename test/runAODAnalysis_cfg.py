import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalysis.AODAnalysis.IFCAAODAnalysis_cff")
#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '106X_upgrade2023_realistic_v2'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/gpfs/users/mantecap/CMSSW_10_6_0_patch2/src/MyAnalysis/AODAnalysis/F0A58BFD-00FF-E443-95A9-FB04741A0179.root'
#        'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring19DR/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/AODSIM/NoPU_106X_upgrade2023_realistic_v3-v2/50000/F0A58BFD-00FF-E443-95A9-FB04741A0179.root'
    )
)





process.p = cms.Path(process.AODanalyzer)


