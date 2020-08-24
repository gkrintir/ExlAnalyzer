import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('ExlAnalyzer',
    recoCaloTower      = cms.InputTag("towerMaker")
)
