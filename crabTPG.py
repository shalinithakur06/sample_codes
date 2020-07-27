from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DataTP2018D_ZeroBias_Run2018D-v1_RAW_Run324725_15Feb19'
config.General.workArea = 'DataTP2018D_ZeroBias_Run2018D-v1_RAW_Run324725_15Feb19'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ecalTPGAnalysisTEMPLATE_weekly.py' 
config.JobType.maxMemoryMB = 2000
config.JobType.outputFiles = ['ECALTPGtree.root']
config.Data.inputDataset = '/ZeroBias8/Run2018D-v1/RAW'
config.Data.inputDBS = 'global'
config.Data.runRange = '324725'
#===== RUNNING OVER ALL RUNS ==========
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
#config.Data.splitting = 'EventBased'
##config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#=====================================
#config.Data.splitting = 'EventBased'
#config.Data.unitsPerJob = 10 
#NJOBS = 1
config.Data.publication = False

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.ignoreLocality = True                                                                                                                                                                       
config.section_("Site")
config.Site.storageSite = 'T2_IN_TIFR'
