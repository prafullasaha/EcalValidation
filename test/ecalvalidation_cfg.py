import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.MessageLogger.cerr = cms.untracked.PSet(threshold = cms.untracked.string("ERROR"))

# Geometry
#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
#process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
#
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff" )
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v16')

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(
#'/store/relval/CMSSW_9_0_0_pre4/RelValMinBias_13/GEN-SIM-RECO/90X_upgrade2017_realistic_v6_C1-v1/00000/4A007E08-21FD-E611-A97E-0CC47A7C35D2.root',
#'/store/mc/RunIISpring15DR74/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM-RECO/NoPURawReco_castor_MCRUN2_74_V8B-v1/00000/0228FAE3-200E-E511-BAA9-485B3919F0A3.root',
#'/store/mc/RunIISpring15DR74/MinBias_TuneZ2star_13TeV-pythia6/GEN-SIM-RECO/NoPU0TRawReco_magnetOff_MCRUN2_740TV0-v1/50000/026AC35A-5B10-E511-9479-0002C90B395A.root',
'/store/data/Run2017A/ZeroBias1/RECO/PromptReco-v3/000/296/663/00000/0C40FCDC-5C58-E711-8393-02163E01A50C.root',
#'/store/relval/CMSSW_8_0_27/ZeroBias/RECO/2017_04_18_18_14_PRnewco_80X_dataRun2_2016LegacyRepro_v3-v1/00000/0050CF28-9124-E711-B363-0025905B85D6.root',
#'/store/relval/CMSSW_8_0_27/ZeroBias/RECO/2017_04_18_18_14_PRnewco_80X_dataRun2_2016LegacyRepro_v3-v1/00000/04318A50-8D24-E711-874D-0025905B858E.root',
#'/store/relval/CMSSW_9_2_2/RelValProdMinBias/GEN-SIM-RECO/91X_mcRun1_realistic_v2-v2/10000/7473BBA6-B14D-E711-8057-0CC47A7452D8.root',
#' /store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00099D43-77ED-E611-8889-5065F381E1A1.root',
#'/store/relval/CMSSW_8_0_0/SingleElectron/RECO/80X_dataRun2_relval_v0_RelVal_sigEl2015B-v1/10000/062034D6-1CDA-E511-B192-0025905A6138.root',
#    '/store/relval/CMSSW_3_5_0/RelValMinBias/GEN-SIM-RECO/MC_3XY_V21-v1/0013/F6E3FFAC-6113-DF11-AF06-001BFCDBD15E.root',
#    '/store/relval/CMSSW_3_5_0/RelValMinBias/GEN-SIM-RECO/MC_3XY_V21-v1/0012/BA183BD4-3813-DF11-B015-00304867902C.root',
#    '/store/relval/CMSSW_3_5_0/RelValMinBias/GEN-SIM-RECO/MC_3XY_V21-v1/0012/A6D53A27-3813-DF11-B460-0026189437FA.root'
#    
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20)
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


process.load("Validation.EcalValidation.ecalvalidation_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('EcalValidation_data_test.root')
)
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
#process.Timing = cms.Service("Timing",
#  summaryOnly = cms.untracked.bool(False),
#  useJobReport = cms.untracked.bool(True)
#)
process.p = cms.Path(
#     process.MessageLogger*
    #process.skimming*
    #process.hltLevel1GTSeed*
    #process.noscraping*
    #process.primaryVertexFilter*
    process.ecalvalidation
    )

