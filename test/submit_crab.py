from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.transferLogs = True
config.General.requestName = 'ZprimeToMuMu_M-5000_ntupler' 
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

config.Data.outputDatasetTag = 'ZprimeToMuMu_M-5000_ntupler'

config.Data.inputDBS = 'phys03'
config.Data.outLFNDirBase = '/store/user/fernanpe/' 
config.Data.publication = False
config.Data.outputDatasetTag = 'ZprimeToMuMu_M-5000_ntupler'

config.Data.inputDataset = '/ZprimeToMuMu_M-5000_TuneCP5_13TeV-madgraphMLM-pythia8/fernanpe-ZprimeToMuMu_M-5000_step2-900573cea688b54324a7e4af95ae2c49/USER'


config.section_("Site")
config.Site.storageSite = 'T2_ES_IFCA'
