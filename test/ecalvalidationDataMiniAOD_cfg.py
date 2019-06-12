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
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_forECALFR_v2')

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(

'/store/data/Run2017C/ZeroBias/MINIAOD/12Sep2017-v1/00000/0230A704-61AA-E711-8E83-A0369FC5FBA4.root',
#'/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/0005AD9F-64ED-E611-A952-0CC47A78A42C.root',
#'/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/00415FAC-B5EC-E611-A1C9-00266CF3E130.root',
#'/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/00763EB4-7FED-E611-AE1F-0CC47A7FC74A.root',
#'/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/009E15D2-E0EC-E611-B01B-0025904B5F8C.root',
#'/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/00A0D9E5-0BED-E611-A89F-0CC47A7E6A4C.root',
#'/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/00F5584D-FDEC-E611-A137-00266CF3DDB4.root',
#'/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/028AC83D-9CED-E611-87E4-0025905B8572.root',
#'/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/02901D3E-87ED-E611-BA53-3417EBE706ED.root',
#'/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/02E7AEDB-1AED-E611-9731-00266CF3E3C4.root',

#'/store/relval/CMSSW_8_0_0/SingleElectron/RECO/80X_dataRun2_relval_v0_RelVal_sigEl2015B-v1/10000/062034D6-1CDA-E511-B192-0025905A6138.root',
#    '/store/relval/CMSSW_3_5_0/RelValMinBias/GEN-SIM-RECO/MC_3XY_V21-v1/0013/F6E3FFAC-6113-DF11-AF06-001BFCDBD15E.root',
#    '/store/relval/CMSSW_3_5_0/RelValMinBias/GEN-SIM-RECO/MC_3XY_V21-v1/0012/BA183BD4-3813-DF11-B015-00304867902C.root',
#    '/store/relval/CMSSW_3_5_0/RelValMinBias/GEN-SIM-RECO/MC_3XY_V21-v1/0012/A6D53A27-3813-DF11-B460-0026189437FA.root'
#    
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


#process.load("Validation.EcalValidation.ecalvalidation_cfi")
process.load("Validation.EcalValidation.ecalvalidationMiniAOD_cfi")
process.ecalvalidation.isMC = cms.bool(False)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('EcalValidation_MiniAODdata_new_cut.root')
)

process.p = cms.Path(
    #process.skimming*
    #process.hltLevel1GTSeed*
    #process.noscraping*
    #process.primaryVertexFilter*
    process.ecalvalidation
    )

