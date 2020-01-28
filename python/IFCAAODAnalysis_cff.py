import FWCore.ParameterSet.Config as cms


AODanalyzer = cms.EDAnalyzer('AODAnalysis',
    nameOfOutput = cms.string('outputAOD.root'),
    MuonCollection = cms.InputTag("muons"),

)


