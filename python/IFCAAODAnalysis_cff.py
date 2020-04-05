import FWCore.ParameterSet.Config as cms

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.TrackingTools.MuonTrackLoader_cff import *

AODanalyzer = cms.EDAnalyzer('AODAnalysis',
    MuonServiceProxy,
    nameOfOutput = cms.string('outputAOD.root'),
    MuonCollection = cms.InputTag("muons"),
    Propagator = cms.string('SmartPropagatorAny'),
    # PickyTrackCollection = cms.InputTag("tevMuons", "picky"),
    # TPFMSTrackCollection = cms.InputTag("tevMuons", "firstHit"),
    # dytTrackCollection = cms.InputTag("tevMuons", "dyt"),
    # GlobalTrackCollection = cms.InputTag("tevMuons", "default"),

)

