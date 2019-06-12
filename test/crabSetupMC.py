from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName  = 'Analysis'
#config.JobType.outputFiles = ['noise_histos.root']


###for data
config.General.requestName = 'NeutrinoMCAOD_12nov_2017'
#config.JobType.psetName = 'ecalvalidation_cfg.py'
#config.JobType.psetName = 'ecalvalidationDATAAOD_cfg.py'
config.JobType.psetName = 'ecalvalidationMCAOD_cfg.py'
#config.Data.useParent = True
#config.JobType.maxMemoryMB = 4000

#config.Data.inputDataset = '/ZeroBias/CMSSW_8_0_27-2017_04_18_18_14_PRnewco_80X_dataRun2_2016LegacyRepro_v3-v1/RECO'
#config.Data.inputDataset = '/ZeroBias/CMSSW_8_0_25-2016_12_21_07_52_PRref_80X_dataRun2_Prompt_v15-v1/RECO'
#config.Data.inputDataset = '/ZeroBias1/Run2017A-PromptReco-v3/RECO'
#config.Data.inputDataset = '/ZeroBias/Run2017A-PromptReco-v1/AOD'
#config.Data.inputDataset = '/MinimumBias/Run2017A-PromptReco-v3/AOD'
config.Data.inputDataset = '/SingleNeutrino/PhaseISpring17DR-FlatPU28to62_90X_upgrade2017_realistic_v20-v1/AODSIM'
#config.Data.secondaryInputDataset = '/SingleNeutrino/PhaseISpring17DR-FlatPU28to62_90X_upgrade2017_realistic_v20-v1/GEN_SIM-RECO'
#config.Data.inputDataset = '/MinBias_TuneCUETP8M1_13TeV-pythia8/PhaseIFall16DR-NoPUNZS_90X_upgrade2017_realistic_v6_C1_ext1-v1/AODSIM'

#
#
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.runRange = '283946'
###config.Data.runRange = '275601-275603'


####for MC
#config.General.requestName = 'ZZTo4L_13TeV_powheg_pythia8'
#config.JobType.psetName = 'ggAnalysis/ggNtuplizer/test/run_mc_80X.py'
#config.Data.inputDataset = '/TTTo2L2Nu_13TeV-powheg/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/MINIAODSIM'
#config.Data.inputDataset = '/WWTo2L2Nu_13TeV-powheg/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
###config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
###config.Data.inputDataset = '/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall15MiniAODv1-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'

#config.Data.inputDataset = '/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
#config.Data.inputDataset = '/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
#config.Data.inputDataset = '/ZZTo4L_13TeV_powheg_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'

#config.Data.inputDataset = '/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'

#config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
#config.JobType.inputFiles = ['Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK8PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK8PFchs.txt', 'Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK8PFchs.txt', 'Summer16_23Sep2016EFV4_DATA_L2Relative_AK8PFchs.txt', 'Summer16_23Sep2016EFV4_DATA_L3Absolute_AK8PFchs.txt', 'Summer16_23Sep2016GV4_DATA_L2L3Residual_AK8PFchs.txt', 'Summer16_23Sep2016GV4_DATA_L2Relative_AK8PFchs.txt', 'Summer16_23Sep2016GV4_DATA_L3Absolute_AK8PFchs.txt', 'Summer16_23Sep2016HV4_DATA_L2L3Residual_AK8PFchs.txt', 'Summer16_23Sep2016HV4_DATA_L2Relative_AK8PFchs.txt', 'Summer16_23Sep2016HV4_DATA_L3Absolute_AK8PFchs.txt', 'Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt', 'Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt', 'Summer16_23Sep2016AllV4_DATA.db', 'Summer16_23Sep2016V4_MC.db']

config.Data.outputDatasetTag = 'Reco_prompt_reprocessing'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 6
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False

config.Site.storageSite ='T2_IN_TIFR' 
