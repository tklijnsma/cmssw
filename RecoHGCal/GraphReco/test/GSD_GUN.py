# coding: utf-8

import os
import math

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from reco_prodtools.templates.GSD_fragment import process

# option parsing
options = VarParsing('python')
options.setDefault('outputFile', 'file:partGun_PDGid22_x96_Pt1.0To100.0_GSD_1.root')
options.setDefault('maxEvents', 1)
options.register("pileup", 0, VarParsing.multiplicity.singleton, VarParsing.varType.int,
    "pileup")
options.register("seed", 1, VarParsing.multiplicity.singleton, VarParsing.varType.int,
    "random seed")
options.register(
    'EminFineTrack',
    10000.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum energy in MeV for which secondary tracks in Calo will be saved"
    )
options.register(
    'EminFinePhoton',
    500.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum energy in MeV for which secondary photons in Calo will be saved"
    )
options.register("doFineCalo", 0, VarParsing.multiplicity.singleton, VarParsing.varType.int,
    "turn do fineCalo on/off")
options.register("storeHGCBoundaryCross", 1, VarParsing.multiplicity.singleton, VarParsing.varType.int,
    "turn do StoreHGCBoundarCross on/off")
options.parseArguments()

process.maxEvents.input = cms.untracked.int32(options.maxEvents)

seed = int(options.seed)+1
# random seeds
process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(seed)
process.RandomNumberGeneratorService.VtxSmeared.initialSeed = cms.untracked.uint32(seed)
process.RandomNumberGeneratorService.mix.initialSeed = cms.untracked.uint32(seed)

# Input source
process.source.firstLuminosityBlock = cms.untracked.uint32(seed)

# Output definition
process.FEVTDEBUGoutput.fileName = cms.untracked.string(
    options.__getattr__("outputFile", noTags=True))

process.FEVTDEBUGoutput.outputCommands.append("keep *_*G4*_*_*")

# helper
def calculate_rho(z, eta):
    return z * math.tan(2 * math.atan(math.exp(-eta)))


process.generator = cms.EDProducer("FlatEtaRangeGunProducer",
    # particle ids
    #particleIDs=cms.vint32(22,22,11,-11,211,-211,13,-13),
    particleIDs=cms.vint32(22),
    # max number of particles to shoot at a time
    nParticles=cms.int32(1),
    # shoot exactly the particles defined in particleIDs in that order
    exactShoot=cms.bool(False),
    # randomly shoot [1, nParticles] particles, each time randomly drawn from particleIDs
    randomShoot=cms.bool(False),
    # energy range
    eMin=cms.double(3.0),
    eMax=cms.double(100.0),
    # phi range
    phiMin=cms.double(-math.pi),
    phiMax=cms.double(math.pi),
    # eta range
    etaMin=cms.double(1.52),
    etaMax=cms.double(3.00),

    debug=cms.untracked.bool(True),
)

# Options for saving fine hits
process.g4SimHits.CaloSD.StoreHGCBoundaryCross = cms.bool(bool(options.storeHGCBoundaryCross))
# Seems to be an interplay between hgcboundaryCross and fineCaloID
process.g4SimHits.CaloSD.UseFineCaloID = cms.bool(bool(options.storeHGCBoundaryCross+options.doFineCalo))
process.g4SimHits.CaloTrkProcessing.DoFineCalo = cms.bool(bool(options.doFineCalo))
process.g4SimHits.CaloTrkProcessing.EminFineTrack = cms.double(options.EminFineTrack)
process.g4SimHits.CaloTrkProcessing.EminFinePhoton = cms.double(options.EminFinePhoton)


#load and configure the appropriate pileup modules
if options.pileup > 0:
    process.load("SimGeneral.MixingModule.mix_POISSON_average_cfi")
    process.mix.input.nbPileupEvents.averageNumber = cms.double(options.pileup)
    # process.mix.input.fileNames = cms.untracked.vstring(["/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/F7FE3FE9-565B-544A-855E-902BA4E3C5FD.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/82584FBA-A1E6-DF48-99BA-B1759C3A190F.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/F806295A-492F-EF4F-9D91-15DA8769DD72.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/6FCA2E1D-D1E2-514B-8ABA-5B71A2C1E1B3.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/287275CC-953A-0C4C-B352-E39EC2D571F0.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/657065A5-F35B-3147-AED9-E4ACA915C982.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/2C56BC73-5687-674C-8684-6C785A88DB78.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/B96F4064-156C-5E47-90A0-07475310157A.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/2564B36D-A0DB-6C42-9105-B1CFF44F311D.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v1/10000/2CB8C960-47C0-1A40-A9F7-0B62987097E0.root"])  # noqa: E501
    local_pu_dir = "/eos/user/m/mrieger/data/hgc/RelvalMinBias14"
    process.mix.input.fileNames = cms.untracked.vstring([
        "file://" + os.path.abspath(os.path.join(local_pu_dir, elem))
        for elem in os.listdir(local_pu_dir)
        if elem.endswith(".root")
    ])
else:
    process.load("SimGeneral.MixingModule.mixNoPU_cfi")
    process.mix.digitizers = cms.PSet(process.theDigitizersValid)
    # I don't think this matters, but just to be safe...
    process.mix.bunchspace = cms.int32(25)
    process.mix.minBunch = cms.int32(-3)
    process.mix.maxBunch = cms.int32(3)

