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
    #tritonPath=cms.string("/afs/cern.ch/work/g/gvonsem/public/HGCAL/ML/HGCalML/triton/"),
    tritonPath=cms.string("RecoHGCal/GraphReco/test/"),
    inpipeName=cms.string("arrayspipe"),
    outpipeName=cms.string("arrayspipe_pred"),
    minEta=cms.double(1.6),
    maxEta=cms.double(3.0)  
)
