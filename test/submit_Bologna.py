from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.transferLogs = True
config.General.requestName = 'SingleMuMinusPt20to2500_NoPU_ntupler' 
config.General.workArea = 'crab_projects'


config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName = 'runRECOAnalysis_cfg.py'
#config.JobType.maxMemoryMB = 5000
#config.JobType.numCores = 8
config.JobType.outputFiles = ['tree.root']

config.section_("Data")
config.Data.splitting   = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.outputDatasetTag = 'SingleMuMinusPt20to2500_NoPU_ntupler'

config.Data.inputDBS = 'phys03'
config.Data.outLFNDirBase = '/store/user/fernanpe/' 
config.Data.publication = False
config.Data.outputDatasetTag = 'SingleMuMinusPt20to2500_NoPU_ntupler'

config.Data.inputDataset = '/SingleMuMinusPt20to2500_NoPU/federica-step3_CMSSW_11_0_X_RECOSIMoutput-b9d107c608087746cef6295f322015f7/USER'


config.section_("Site")
config.Site.storageSite = 'T2_ES_IFCA'
