import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.MessageLogger.cerr = cms.untracked.PSet(threshold = cms.untracked.string("ERROR"))

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v8')
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_forECALFR_v2')
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
#    'file:/tmp/malberti/901E9460-4E0E-E011-9D0D-0015170AE6E4.root'
#'/store/data/Run2017A/ZeroBias1/AOD/PromptReco-v3/000/296/663/00000/421AE213-A452-E711-A0BD-02163E0135E0.root'
#'/store/mc/RunIISpring15DR74/MinBias_TuneZ2star_13TeV-pythia6/MINIAODSIM/NoPURaw_castor_MCRUN2_74_V8-v1/50000/04F68838-31FA-E411-B993-00259073E398.root'
#'/store/relval/CMSSW_9_2_15/RelValNuGun/MINIAODSIM/PU25ns_92X_upgrade2017_realistic_forECALFR_v2_HS_PF16-v1/20000/56DC1B16-3AE1-E711-89BD-4C79BA3203F7.root',
' /store/relval/CMSSW_9_2_15/RelValNuGun/MINIAODSIM/PU25ns_92X_upgrade2017_realistic_forECALFR_v2_HS_PF16-v1/20000/56DC1B16-3AE1-E711-89BD-4C79BA3203F7.root',
#'/store/mc/RunIISpring15MiniAODv2/MinBias_TuneMBR_13TeV-pythia8/MINIAODSIM/castor_74X_mcRun2_asymptotic_v2-v1/10000/0AF6B2DF-E17B-E511-AE58-0CC47A13CEF4.root',
#'/store/data/Run2017A/MinimumBias/AOD/PromptReco-v3/000/296/662/00000/D4ED9F54-A452-E711-A3A7-02163E014141.root',
#'/store/relval/CMSSW_8_0_0/RelValMinBias_13/MINIAODSIM/80X_mcRun2_asymptotic_v4-v1/10000/9AF3A324-15DA-E511-90E3-0025905A6082.root',
#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00099D43-77ED-E611-8889-5065F381E1A1.root',
#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/008B4E7F-8FED-E611-897F-5065F381F291.root',
#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00BF45FC-5AED-E611-986A-A0000420FE80.root',
#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/02A4AE2C-71ED-E611-B8B1-24BE05C38CA1.root',
#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/02CC705C-70ED-E611-B671-5065F37DD491.root',
#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/08D13E45-77ED-E611-BD0B-5065F381A251.root',
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



process.load("Validation.EcalValidation.ecalvalidationMiniAOD_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('EcalValidationMiniAOD_new_cut.root')
)

process.p = cms.Path(
    #process.skimming*
    #process.hltLevel1GTSeed*
    #process.noscraping*
    #process.primaryVertexFilter*
     process.ecalvalidation
    )

