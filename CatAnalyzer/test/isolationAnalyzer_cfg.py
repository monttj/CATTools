import FWCore.ParameterSet.Config as cms
process = cms.Process("IsoAna")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for PHYS14 scenario PU4bx50 : global tag is ???
#    for PHYS14 scenario PU20bx25: global tag is PHYS14_25_V1
#  as a rule, find the global tag in the DAS under the Configs for given dataset
#process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/AsymptNoPURawReco_MCRUN2_74_V9A-v4/10000/F43B9337-811A-E511-B426-008CFA56D870.root'
    )
)

#######################################################################
## for egamma pid temp 
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_74X
from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat,switchOnVIDElectronIdProducer,setupAllVIDIdsInModule,setupVIDElectronSelection

useMiniAOD = True

if not useMiniAOD :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.IsoAna = cms.EDAnalyzer('IsolationAnalyzer',
                                           genLabel      = cms.InputTag("prunedGenParticles"),
                                           vertexLabel   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                           muonLabel     = cms.InputTag("slimmedMuons"),
                                           electronLabel = cms.InputTag("slimmedElectrons"),
                                           # ID decisions (common to all formats)
                                           eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                                           eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                                           )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('cmspfa.root')
                                   )


process.p = cms.Path(process.egmGsfElectronIDSequence*process.IsoAna)
