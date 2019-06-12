import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.MessageLogger.cerr = cms.untracked.PSet(threshold = cms.untracked.string("ERROR"))

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v8')
process.GlobalTag = GlobalTag(process.GlobalTag, '90X_upgrade2017_realistic_v20')
process.load("Configuration.StandardSequences.MagneticField_cff")

## Geometry
##process.load("Configuration.Geometry.Geometry_cff") #old
##process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
##process.load('Configuration.Geometry.GeometryIdeal_cff')
#
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
##process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
#
## initialize magnetic field
#process.load("Configuration.StandardSequences.MagneticField_cff")
##process.load("MagneticField.Engine.autoMagneticFieldProducer_cfi")
#
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(
    'file:SingleNeutrino_flatpileup28to42_AOD_1rootfile.root',
#'/store/mc/PhaseISpring17DR/SingleNeutrino/AODSIM/FlatPU28to62_90X_upgrade2017_realistic_v20-v1/60000/00EC273F-682E-E711-A96C-549F3525DFE8.root',
#'/store/mc/PhaseIFall16DR/MinBias_TuneCUETP8M1_13TeV-pythia8/AODSIM/NoPUNZS_90X_upgrade2017_realistic_v6_C1_ext1-v1/120000/0AEFBB63-6C0A-E711-8737-02163E01A74F.root',
#'/store/relval/CMSSW_9_2_2/RelValProdMinBias/GEN-SIM-RECO/91X_mcRun1_realistic_v2-v2/10000/7473BBA6-B14D-E711-8057-0CC47A7452D8.root',
#'/store/relval/CMSSW_9_2_2/RelValProdMinBias/GEN-SIM-RECO/91X_mcRun1_realistic_v2-v2/10000/C0124DA2-B14D-E711-8D0B-0CC47A4D769A.root',
#'file:/afs/cern.ch/user/s/sdutt/work/24F580CC-D0F3-E611-959E-0CC47A4D7634.root',
#'file:/afs/cern.ch/user/s/sdutt/work/30C0FB99-CFF3-E611-B72F-0025905A6090.root',
#'file:/afs/cern.ch/user/s/sdutt/work/381B0902-CFF3-E611-BDA3-0CC47A4D768C.root',
#'file:/afs/cern.ch/user/s/sdutt/work/50EEC11D-CEF3-E611-A3C3-0CC47A4D76C0.root',
#'file:/afs/cern.ch/user/s/sdutt/work/76B0F623-D3F3-E611-9CE2-0025905A610A.root',
#'file:/afs/cern.ch/user/s/sdutt/work/86DE3B1E-CEF3-E611-8355-0CC47A4D768E.root',
#'file:/afs/cern.ch/user/s/sdutt/work/960A6721-D2F3-E611-B94E-0CC47A4C8EC8.root',
#'file:/afs/cern.ch/user/s/sdutt/work/96EFEEBA-D0F3-E611-B6C3-0CC47A78A360.root',
#'file:/afs/cern.ch/user/s/sdutt/work/C2274CF0-D4F3-E611-AA56-0025905B856E.root',
#'file:/afs/cern.ch/user/s/sdutt/work/DAC11501-D5F3-E611-9B89-0CC47A78A3E8.root',

#'/store/relval/CMSSW_8_0_26/RelValTTbar_13/GEN-SIM-RECO/80X_mcRun2_asymptotic_2016_TrancheIV_v8_width1p9mum_BS1p9-v1/00000/24F580CC-D0F3-E611-959E-0CC47A4D7634.root',
    #'/store/relval/CMSSW_3_5_0/RelValMinBias/GEN-SIM-RECO/MC_3XY_V21-v1/0013/F6E3FFAC-6113-DF11-AF06-001BFCDBD15E.root',
    #'/store/relval/CMSSW_3_5_0/RelValMinBias/GEN-SIM-RECO/MC_3XY_V21-v1/0012/BA183BD4-3813-DF11-B015-00304867902C.root',
    #'/store/relval/CMSSW_3_5_0/RelValMinBias/GEN-SIM-RECO/MC_3XY_V21-v1/0012/A6D53A27-3813-DF11-B460-0026189437FA.root'
     
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-6)
)


# filter on PhysDeclared bit
process.skimming = cms.EDFilter("PhysDecl",
    applyfilter = cms.untracked.bool(True)
)

# filter on bit 40 || 41 nad !(bit36 || bit37 || bit38 || bit39)
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

# FilterOutScraping
process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )

# Good Vertex Filter
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNumberOfTracks = cms.uint32(3) ,
                                           maxAbsZ = cms.double(15),	
                                           maxd0 = cms.double(2)	
                                           )



process.load("Validation.EcalValidation.ecalvalidationAOD_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('EcalValidation_MC_AOD.root')
)

process.p = cms.Path(
    #process.skimming*
    #process.hltLevel1GTSeed*
    #process.noscraping*
    #process.primaryVertexFilter*
     process.ecalvalidation
    )

