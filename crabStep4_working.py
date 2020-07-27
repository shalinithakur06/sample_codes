from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Run2016C_test'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ggAnalysis/ggNtuplizer/test/run_data_80X.py'
#config.JobType.psetName = 'ggAnalysis/ggNtuplizer/test/run_mc_80X.py'

#config.Data.inputDataset = '/DoubleMuon/Run2016B-PromptReco-v1/MINIAOD'
config.Data.inputDataset = '/DoubleMuon/Run2016C-PromptReco-v2/MINIAOD'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall15MiniAODv1-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
config.Data.inputDBS = 'global'


config.JobType.inputFiles = ['Spring16_25nsV6_DATA.db','Spring16_25nsV6_DATA_L3Absolute_AK8PFchs.txt','Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt','Spring16_25nsV6_DATA_L2L3Residual_AK8PFchs.txt','Spring16_25nsV6_MC.db','Spring16_25nsV6_DATA_L2Relative_AK8PFchs.txt','Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt']


config.Data.outputDatasetTag = 'Run2016C_test'
#config.Data.runRange = '275601-275603'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 15
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

config.Site.storageSite ='T2_IN_TIFR' 
