import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
# NOTE: the pick the right global tag!
# for PHYS14 scenario PU20bx25: global tag is PHYS14_25_V1
# as a rule, find the global tag in the DAS under the Configs for given dataset
#process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
process.GlobalTag.globaltag = '74X_dataRun2_reMiniAOD_v0'
#74X_mcRun2_asymptotic_v2'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#      'file:/eos/uscms/store/user/lammel/r2.2015/chi250_stau225_lsp200/evts_7aod_10k.root',
#       'file:/eos/uscms/store/user/lammel/mg5_13tev_chi125_stau110_lsp095/events_7aod.root' 
        'root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/00759690-D16E-E511-B29E-00261894382D.root'
)
)


# START ELECTRON ID SECTION
#
# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
# overwrite a default parameter: for miniAOD, the collection name is a slimmed one
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')
                             
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']
for idmod in my_id_modules:  
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
 


process.TFileService = cms.Service("TFileService",
fileName = cms.string("OutTree.root")
)

process.TNT = cms.EDAnalyzer("TriggerEff",
    # boolean variables
                             triggerResults      = cms.InputTag( 'TriggerResults', '', 'HLT' ),
                             HLTPath1 =  cms.string( "HLT_IsoMu18_v1"),  #HLT_IsoMu24_eta2p1_v1" ),
                             HLTFilter1a= cms.string( "hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09"), #hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09" ),
                             HLTFilter1b= cms.string( "hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09"), # hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09" ),
#HLTPath1 =  cms.string( "HLT_IsoMu27_v3" ),
#                             HLTFilter1a= cms.string( "hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09" ),
#                             HLTFilter1b= cms.string( "hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09" ),
                             HLTPath2 =  cms.string( "HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v3"),

                             HLTFilter2a = cms.string( "hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09" ),
                             HLTFilter2b = cms.string( "hltOverlapFilterIsoMu17MediumIsoPFTau40Reg" ),
                             electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                             electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                             electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                             electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                             eleHEEPIdMap        = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
                             isOSCharge = cms.bool(True),     
			     l1min = cms.InputTag("patTrigger","l1min"),
			     l1max = cms.InputTag("patTrigger","l1max"),           
			     TauPtCut = cms.double(20.0),
                             TauEtaCut = cms.double(2.1),
                             TauDMF = cms.string("decayModeFinding"),
                             TauEleVeto = cms.string("againstElectronVLooseMVA5"),
                             TauMuVeto = cms.string("againstMuonTight3"),
                             TauIsoString = cms.string("byLooseCombinedIsolationDeltaBetaCorr3Hits"),
                             DYOthersBG = cms.bool(False),
			     TauIsoCutMax = cms.double(99999),
                             TauIsoCutMin = cms.double(0.5),
                             MuonPtCut = cms.double(19.0),
                             MuonEtaCut = cms.double(2.1),
                             IsoMuonMax = cms.double(0.1),
                             MotherpdgID  = cms.double(23),
			     isMC = cms.bool(False),
			     isZtau  = cms.bool(True),
			     isZprime  = cms.bool(False),
			     GenReq = cms.bool(False),
                             )

process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")
                              
process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)                             
                              
process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
   reverseDecision = cms.bool(False)
)                             
 
process.p = cms.Path(         
#                     process.HBHENoiseFilterResultProducer * 
#                     process.ApplyBaselineHBHENoiseFilter *    
                     process.egmGsfElectronIDSequence *
                     process.TNT)                 

