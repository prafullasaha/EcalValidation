import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.MessageLogger.cerr = cms.untracked.PSet(threshold = cms.untracked.string("ERROR"))

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v16')
#process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v5')
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_2017Repro_v4')
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_Prompt_Candidate_forTkDPG_v2')
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
    #'file:ZeroBias_Run2017B-04Jul2017-v2_AOD_1rootfile.root'
#'/store/data/Run2017C/ZeroBias/AOD/12Sep2017-v1/00000/000CD359-A3A8-E711-B551-D8D385FF319D.root'
#'/store/data/Run2017B/ZeroBias/AOD/06Jul2017-v2/130000/184D4248-BA62-E711-8524-0CC47A7C3404.root',
#'/store/data/Run2017C/ZeroBias/AOD/12Sep2017-v1/70000/EAA6AC00-8FA3-E711-B69E-1866DAEA8230.root',
 '/store/data/Run2017F/ZeroBias/AOD/06Nov2017-v1/150000/00DE8807-03C3-E711-A5A0-D067E5F914D3.root',    
 
   ),
    inputCommands = cms.untracked.vstring(
                  'keep *',
                  'drop CTPPSPixelClusteredmDetSetVector_*_*_*',
                  #'keep FEDRawDataCollection_*_*_*'
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
process.ecalvalidation.isMC = cms.bool(False)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('EcalValidation_Data_AOD.root')
)

process.p = cms.Path(
    #process.skimming*
    #process.hltLevel1GTSeed*
    process.noscraping*
    #process.primaryVertexFilter*
     process.ecalvalidation
    )

