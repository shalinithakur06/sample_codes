import FWCore.ParameterSet.Config as cms
from MiniTree.Selection.LocalRunSkeleton_cff import *
from MiniTree.Selection.LocalSources_cff import toPrint

#-----------------------------
#INPUT & OUTPUT
#-----------------------------
isData=False
#Data
#inFile = "/store/data/Run2016B/SingleMuon/MINIAOD/23Sep2016-v3/00000/00AE0629-1F98-E611-921A-008CFA1112CC.root"
#MC
#inFile = "/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root"
#inFile = "/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/0EA60289-18C4-E611-8A8F-008CFA110AB4.root"
#inFile = "file:step5_Mu2000_1.root"
#inFile = "file:MuMC_2000_20180227.root"
inFile = "/store/user/sthakur/TestLHEGeneration_Mu/MCGenerationStep5_Mu2000_2018_1_27/180128_024426/0000/step5_Mu2000_1.root"
process.source.fileNames = [inFile_]
process.maxEvents.input = cms.untracked.int32(1000)
#for multi CRAB
process.TFileService.fileName = cms.string("outFile_.root")

#-----------------------------
#CONFIG PARAMETERS 
#-----------------------------
procName='LOCALUSER'
#trigMenu = 'HLT2' #https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
trigMenu = 'HLT'
isFastsim = False
trigpath = ''

#-----------------------------
#START PROCESS CONFIGURATION 
#-----------------------------
process.setName_(procName)
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag  = cms.string('80X_mcRun2_asymptotic_2016_TrancheIV_v6')

#-----------------------------
#CONFIGURE THE EXTRA MODULE 
#-----------------------------
storeOutPath=False
producePDFweights=False
defineBasePreSelection(process,False, False)

#-----------------------------
# ADD THE ANALYSIS MODULE 
#-----------------------------
process.load('MiniTree.Selection.selection_cfi')
process.myMiniTreeProducer.MCTruth.isData = cms.bool(isData)
if isData:
    process.myMiniTreeProducer.MCTruth.sampleCode = cms.string("DATA")
else:
    #for multi CRAB
    process.myMiniTreeProducer.MCTruth.sampleCode = cms.string("sampCode_")
process.myMiniTreeProducer.MCTruth.producePDFweights = cms.bool(producePDFweights)
#process.myMiniTreeProducer.minEventQualityToStore = cms.int32(1)
process.myMiniTreeProducer.Trigger.source = cms.InputTag('TriggerResults::'+trigMenu)
process.myMiniTreeProducer.Trigger.myTrig = "HLT"
#process.myMiniTreeProducer.Trigger.myTrig = "HLT_Mu17_TrkIsoVVL"
process.myMiniTreeProducer.MCTruth.sampleChannel = cms.string('muon')
#https://github.com/cms-jet/JRDatabase/tree/master/textFiles/Spring16_25nsV10_MC
process.myMiniTreeProducer.Jets.resolutionsFile = cms.string('Spring16_25nsV10_MC_PtResolution_AK8PF.txt')
process.myMiniTreeProducer.Jets.scaleFactorsFile = cms.string('Spring16_25nsV10_MC_SF_AK8PF.txt')

#-----------------------------
#ANALYSIS SEQUENCE 
#-----------------------------
process.p  = cms.Path(process.allEventsFilter*process.basePreSel*process.myMiniTreeProducer)
process.schedule = cms.Schedule(process.p)
checkProcessSchedule(storeOutPath,True)

