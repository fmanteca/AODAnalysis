import FWCore.ParameterSet.Config as cms


AODanalyzer = cms.EDAnalyzer('AODAnalysis',
    nameOfOutput = cms.string('outputAOD.root'),
    MuonCollection = cms.InputTag("muons"),
    # PickyTrackCollection = cms.InputTag("tevMuons", "picky"),
    # TPFMSTrackCollection = cms.InputTag("tevMuons", "firstHit"),
    # dytTrackCollection = cms.InputTag("tevMuons", "dyt"),
    # GlobalTrackCollection = cms.InputTag("tevMuons", "default"),

)


