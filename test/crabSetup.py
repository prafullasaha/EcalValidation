from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName  = 'Analysis'
#config.JobType.outputFiles = ['noise_histos.root']


###for data
config.General.requestName = 'Run2017CMiniAOD_18_new_Cut'
#config.JobType.psetName = 'ecalvalidationDATAAOD_cfg.py'
config.JobType.psetName = 'ecalvalidationDataMiniAOD_cfg.py'

#config.Data.useParent = True
#config.JobType.maxMemoryMB = 4000

config.Data.inputDataset = '/ZeroBias/Run2017C-12Sep2017-v1/MINIAOD'



#config.Data.inputDataset = '/ZeroBias/Run2017B-06Jul2017-v2/AOD'
#config.Data.inputDataset = '/SingleNeutrino/PhaseISpring17DR-FlatPU28to62_90X_upgrade2017_realistic_v20-v1/AODSIM'
#config.Data.secondaryInputDataset = '/SingleNeutrino/PhaseISpring17DR-FlatPU28to62_90X_upgrade2017_realistic_v20-v1/GEN_SIM-RECO'
#config.Data.inputDataset = '/MinBias_TuneCUETP8M1_13TeV-pythia8/PhaseIFall16DR-NoPUNZS_90X_upgrade2017_realistic_v6_C1_ext1-v1/AODSIM'

#
#
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
#config.Data.runRange = '283946'
###config.Data.runRange = '275601-275603'


config.Data.inputDBS = 'global'
#config.JobType.inputFiles = ['Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK8PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFchs.txt']

config.Data.outputDatasetTag = 'Jan14'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False

config.Site.storageSite ='T2_IN_TIFR' 
