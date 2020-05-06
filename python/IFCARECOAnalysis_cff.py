import FWCore.ParameterSet.Config as cms

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.TrackingTools.MuonTrackLoader_cff import *

RECOanalyzer = cms.EDAnalyzer('RECOAnalysis',
    MuonServiceProxy,
    nameOfOutput = cms.string('tree.root'),
    MuonCollection = cms.InputTag("muons"),
    Propagator = cms.string('SmartPropagatorAny'),
    segmentsDt = cms.InputTag('dt4DSegments'),
    segmentsCSC = cms.InputTag('cscSegments'),
    genParticleCollection = cms.InputTag("genParticles"),
    inputMuonShowerInformationValueMap = cms.InputTag("muons","muonShowerInformation"),
    # PickyTrackCollection = cms.InputTag("tevMuons", "picky"),
    # TPFMSTrackCollection = cms.InputTag("tevMuons", "firstHit"),
    # dytTrackCollection = cms.InputTag("tevMuons", "dyt"),
)

