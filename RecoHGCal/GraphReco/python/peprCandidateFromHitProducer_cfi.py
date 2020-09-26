# coding: utf-8

"""
Initialization file for the WindowInference module.
"""


__all__ = ["peprCandidateFromHitProducer"]


import math

import FWCore.ParameterSet.Config as cms


peprCandidateFromHitProducer = cms.EDProducer("peprCandidateFromHitProducer",
    # the collections of rechits to use
    recHitCollections=cms.VInputTag(
        cms.InputTag("HGCalRecHit", "HGCEERecHits"),
        cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
        cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    ),
    tracks = cms.InputTag("generalTracks"),
    simClusters = cms.InputTag("mix", "MergedCaloTruth"),
    # whether or not the model in the graph expects a batch dimension
    batchedModel=cms.bool(True),
    # graph to the trained model
    #graphPath=cms.string("graph.pb"),
    tritonPath=cms.string("/afs/cern.ch/work/g/gvonsem/public/HGCAL/ML/HGCalML/triton/"),
    inpipeName=cms.string("arrayspipe"),
    outpipeName=cms.string("arrayspipe_pred"),
    minEta=cms.double(1.6),
    maxEta=cms.double(3.0),
    # # window size in phi and eta
    # etaFrameWidth=cms.double(0.1),
    # phiFrameWidth=cms.double(0.1),
    # # overlap in phi and eta
    # nEtaSegments=cms.uint32(3),
    # nPhiSegments=cms.uint32(1),
    # window size in phi and eta
    etaFrameWidth=cms.double(100),
    phiFrameWidth=cms.double(100),
    # overlap in phi and eta
    nEtaSegments=cms.uint32(1),
    nPhiSegments=cms.uint32(1),   
)
