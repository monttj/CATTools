from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'pfpaper_isolation'
config.General.workArea = 'crab_cmspfa_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'isolationAnalyzer_cfg.py'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
#config.Data.outLFNDirBase = '/store/user/%s/pf/' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '/store/user/tjkim/pf/' 
config.Data.publication =False
config.Data.outputDatasetTag = 'pf_iso'

config.Site.storageSite = 'T2_CH_CERN'

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    for dataset in ['/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-AsymptNoPURawReco_MCRUN2_74_V9A-v4/MINIAODSIM', '/QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-AsymptNoPUReco_MCRUN2_74_V9A-v2/MINIAODSIM','/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-AsymptNoPUReco_MCRUN2_74_V9A-v1/MINIAODSIM']:
        config.Data.inputDataset = dataset
        config.General.requestName = dataset.split('/')[1] + "_v4" 
        crabCommand('submit', config = config)
