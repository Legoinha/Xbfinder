### HiForest Configuration
# Input: miniAOD
# Type: data

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import Run3_pp_on_PbPb_2023
process = cms.Process('HiForest', Run3_pp_on_PbPb_2023)

process.options = cms.untracked.PSet(
        numberOfThreads = cms.untracked.uint32(4)
)
process.options.numberOfConcurrentLuminosityBlocks = 1
###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 132X, data")

# import subprocess, os
# version = subprocess.check_output(
#     ['git', '-C', os.path.expandvars('$CMSSW_BASE/src'), 'describe', '--tags'])
# if version == '':
#     version = 'no git info'
# process.HiForestInfo.HiForestVersion = cms.string(version)

###############################################################################

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        'root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v2/000/375/790/00000/56ad580f-b228-4f3c-b8e3-17f9d95c7654.root'
    ), 
)

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
    )

###############################################################################

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v7', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag

###############################################################################

# Define centrality binning
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")

###############################################################################

# root output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("HiForestMiniAOD.root"))

# # edm output for debugging purposes
# process.output = cms.OutputModule(
#     "PoolOutputModule",
#     fileName = cms.untracked.string('HiForestEDM.root'),
#     outputCommands = cms.untracked.vstring(
#         'keep *',
#         )
#     )

# process.output_path = cms.EndPath(process.output)

###############################################################################

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

#process.hiEvtAnalyzer.doCentrality = cms.bool(False)
#process.hiEvtAnalyzer.doHFfilters = cms.bool(False)

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data_2023_skimmed
process.hltobject.triggerNames = trigger_list_data_2023_skimmed

#process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
################################
# electrons, photons, muons
#process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
#process.ggHiNtuplizer.doMuons = cms.bool(False)
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
################################
# jet reco sequence
#process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_data_cff')
#process.load('HeavyIonsAnalysis.JetAnalysis.akPu4CaloJetSequence_pponPbPb_data_cff')
#process.akPu4CaloJetAnalyzer.doHiJetID = True
################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")
# muons (FTW)
process.load("HeavyIonsAnalysis.MuonAnalysis.unpackedMuons_cfi")
#process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
###############################################################################
# ZDC RecHit Producer
process.load("RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi")
process.load('HeavyIonsAnalysis.ZDCAnalysis.ZDCAnalyzersPbPb_cff')
###############################################################################

# main forest sequence
process.forest = cms.Path(
    process.HiForestInfo +
    process.centralityBin +
    process.hiEvtAnalyzer +
    process.hltanalysis +
    process.hltobject +
    process.l1object +
    process.trackSequencePbPb +
    #process.particleFlowAnalyser +
    #process.ggHiNtuplizer +
    #process.zdcdigi +
    #process.QWzdcreco +
    #process.zdcanalyzer +
    process.unpackedMuons #+
    #process.muonAnalyzer +
    #process.akPu4CaloJetAnalyzer
    )

#customisation

# Select the types of jets filled
addR3Jets = False
addR3FlowJets = False
addR4Jets = True
addR4FlowJets = True
addUnsubtractedR4Jets = True

# Choose which additional information is added to jet trees
doHIJetID = True             # Fill jet ID and composition information branches
doWTARecluster = True        # Add jet phi and eta for WTA axis

# this is only for non-reclustered jets
addCandidateTagging = False
'''
if addR3Jets or addR3FlowJets or addR4Jets or addR4FlowJets or addUnsubtractedR4Jets :
    process.load("HeavyIonsAnalysis.JetAnalysis.extraJets_cff")
    from HeavyIonsAnalysis.JetAnalysis.clusterJetsFromMiniAOD_cff import setupHeavyIonJets
    process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")

    if addR3Jets :
        process.jetsR3 = cms.Sequence()
        setupHeavyIonJets('akCs3PF', process.jetsR3, process, isMC = 0, radius = 0.30, JECTag = 'AK3PF', doFlow = False)
        process.akCs3PFpatJetCorrFactors.levels = ['L2Relative', 'L2L3Residual']
        process.akCs3PFJetAnalyzer = process.akCs4PFJetAnalyzer.clone(jetTag = "akCs3PFpatJets", jetName = 'akCs3PF', doHiJetID = doHIJetID, doWTARecluster = doWTARecluster)
        process.forest += process.extraJetsData * process.jetsR3 * process.akCs3PFJetAnalyzer

    if addR3FlowJets :
        process.jetsR3flow = cms.Sequence()
        setupHeavyIonJets('akCs3PFFlow', process.jetsR3flow, process, isMC = 0, radius = 0.30, JECTag = 'AK3PF', doFlow = True)
        process.akCs3PFFlowpatJetCorrFactors.levels = ['L2Relative', 'L2L3Residual']
        process.akFlowPuCs3PFJetAnalyzer = process.akCs4PFJetAnalyzer.clone(jetTag = "akCs3PFFlowpatJets", jetName = 'akCs3PFFlow', doHiJetID = doHIJetID, doWTARecluster = doWTARecluster)
        process.forest += process.extraFlowJetsData * process.jetsR3flow * process.akFlowPuCs3PFJetAnalyzer

    if addR4Jets :
        # Recluster using an alias "0" in order not to get mixed up with the default AK4 collections
        process.jetsR4 = cms.Sequence()
        setupHeavyIonJets('akCs0PF', process.jetsR4, process, isMC = 0, radius = 0.40, JECTag = 'AK4PF', doFlow = False)
        process.akCs0PFpatJetCorrFactors.levels = ['L2Relative', 'L2L3Residual']
        process.akCs4PFJetAnalyzer.jetTag = 'akCs0PFpatJets'
        process.akCs4PFJetAnalyzer.jetName = 'akCs0PF'
        process.akCs4PFJetAnalyzer.doHiJetID = doHIJetID
        process.akCs4PFJetAnalyzer.doWTARecluster = doWTARecluster
        process.forest += process.extraJetsData * process.jetsR4 * process.akCs4PFJetAnalyzer

    if addR4FlowJets :
        process.jetsR4flow = cms.Sequence()
        setupHeavyIonJets('akCs4PFFlow', process.jetsR4flow, process, isMC = 0, radius = 0.40, JECTag = 'AK4PF', doFlow = True)
        process.akCs4PFFlowpatJetCorrFactors.levels = ['L2Relative', 'L2L3Residual']
        process.akFlowPuCs4PFJetAnalyzer.jetTag = 'akCs4PFFlowpatJets'
        process.akFlowPuCs4PFJetAnalyzer.jetName = 'akCs4PFFlow'
        process.akFlowPuCs4PFJetAnalyzer.doHiJetID = doHIJetID
        process.akFlowPuCs4PFJetAnalyzer.doWTARecluster = doWTARecluster
        process.forest += process.extraFlowJetsData * process.jetsR4flow * process.akFlowPuCs4PFJetAnalyzer

    if addUnsubtractedR4Jets:
        process.load('HeavyIonsAnalysis.JetAnalysis.ak4PFJetSequence_ppref_data_cff')
        from HeavyIonsAnalysis.JetAnalysis.clusterJetsFromMiniAOD_cff import setupPprefJets
        process.unsubtractedJetR4 = cms.Sequence()
        setupPprefJets('ak04PF', process.unsubtractedJetR4, process, isMC = 0, radius = 0.40, JECTag = 'AK4PF')
        process.ak04PFpatJetCorrFactors.levels = ['L2Relative', 'L2L3Residual']
        process.ak4PFJetAnalyzer.jetTag = "ak04PFpatJets"
        process.ak4PFJetAnalyzer.jetName = "ak04PF"
        process.forest += process.unsubtractedJetR4 * process.ak4PFJetAnalyzer

if addCandidateTagging:
    process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")

    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
        btagDiscriminators = ['pfCombinedSecondaryVertexV2BJetTags', 'pfDeepCSVDiscriminatorsJetTags:BvsAll', 'pfDeepCSVDiscriminatorsJetTags:CvsB', 'pfDeepCSVDiscriminatorsJetTags:CvsL'], ## to add discriminators,
        btagPrefix = 'TEST',
    )

    process.updatedPatJets.addJetCorrFactors = False
    process.updatedPatJets.discriminatorSources = cms.VInputTag(
        cms.InputTag('pfDeepCSVJetTags:probb'),
        cms.InputTag('pfDeepCSVJetTags:probc'),
        cms.InputTag('pfDeepCSVJetTags:probudsg'),
        cms.InputTag('pfDeepCSVJetTags:probbb'),
    )

    process.akCs4PFJetAnalyzer.jetTag = "updatedPatJets"

    process.forest.insert(1,process.candidateBtagging*process.updatedPatJets)
'''
#########################
# Event Selection -> add the needed filters here
#########################

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')
process.pphfCoincFilter4Th2 = cms.Path(process.phfCoincFilter4Th2)
process.pphfCoincFilter1Th3 = cms.Path(process.phfCoincFilter1Th3)
process.pphfCoincFilter2Th3 = cms.Path(process.phfCoincFilter2Th3)
process.pphfCoincFilter3Th3 = cms.Path(process.phfCoincFilter3Th3)
process.pphfCoincFilter4Th3 = cms.Path(process.phfCoincFilter4Th3)
process.pphfCoincFilter5Th3 = cms.Path(process.phfCoincFilter5Th3)
process.pphfCoincFilter1Th4 = cms.Path(process.phfCoincFilter1Th4)
process.pphfCoincFilter2Th4 = cms.Path(process.phfCoincFilter2Th4)
process.pphfCoincFilter3Th4 = cms.Path(process.phfCoincFilter3Th4)
process.pphfCoincFilter4Th4 = cms.Path(process.phfCoincFilter4Th4)
process.pphfCoincFilter5Th4 = cms.Path(process.phfCoincFilter5Th4)
process.pphfCoincFilter1Th5 = cms.Path(process.phfCoincFilter1Th5)
process.pphfCoincFilter2Th5 = cms.Path(process.phfCoincFilter2Th5)
process.pphfCoincFilter3Th5 = cms.Path(process.phfCoincFilter3Th5)
process.pphfCoincFilter4Th5 = cms.Path(process.phfCoincFilter4Th5)
process.pphfCoincFilter5Th5 = cms.Path(process.phfCoincFilter5Th5)

process.pAna = cms.EndPath(process.skimanalysis)

#from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
#process.hltfilter = hltHighLevel.clone(
#    HLTPaths = [
#        #"HLT_HIZeroBias_v4",                                                     
#        "HLT_HIMinimumBias_v2",
#    ]
#)
#process.filterSequence = cms.Sequence(
#    process.hltfilter
#)
#
#process.superFilterPath = cms.Path(process.filterSequence)
#process.skimanalysis.superFilters = cms.vstring("superFilterPath")
#
#for path in process.paths:
#    getattr(process, path)._seq = process.filterSequence * getattr(process,path)._seq


#################### B finder #################
runOnMC = False ## !!
VtxLabel = "offlineSlimmedPrimaryVertices"
TrkLabel = "packedPFCandidates"
TrkChi2Label = "packedPFCandidateTrackChi2"
GenLabel = "prunedGenParticles"

from Bfinder.finderMaker.finderMaker_75X_cff import finderMaker_75X
finderMaker_75X(process, runOnMC, VtxLabel, TrkLabel, TrkChi2Label, GenLabel)
process.Bfinder.MVAMapLabel = cms.InputTag(TrkLabel,"MVAValues")
process.Bfinder.tkPtCut = cms.double(1) # before fit
process.Bfinder.tkEtaCut = cms.double(2.4) # before fit
process.Bfinder.jpsiPtCut = cms.double(0.0) # before fit
process.Bfinder.bPtCut = cms.vdouble(3.0, 5.0, 5.0, 5.0, 5.0, 3.0, 1.0) # before fit
process.Bfinder.Bchannel = cms.vint32(1, 0, 0, 1, 1, 1, 0)
process.Bfinder.VtxChiProbCut = cms.vdouble(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
process.Bfinder.svpvDistanceCut = cms.vdouble(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0)
process.Bfinder.doTkPreCut = cms.bool(True)
process.Bfinder.doMuPreCut = cms.bool(True)
process.Bfinder.MuonTriggerMatchingPath = cms.vstring("")
process.Bfinder.MuonTriggerMatchingFilter = cms.vstring("")
process.BfinderSequence.insert(0, process.unpackedMuons)
process.BfinderSequence.insert(0, process.unpackedTracksAndVertices)
# process.unpackedMuons.muonSelectors = cms.vstring() # uncomment for pp

process.Bfinder.centmin = cms.double(0)
process.Bfinder.centmax = cms.double(100)

process.p = cms.Path(process.BfinderSequence)




#######################################################################################################################
#######################################################################################################################
# Muon filtering before running Bfinder to significantly speed up the processing
MUONCUT = "isTrackerMuon && ((abs(eta) <= 1.0 && pt > 3.5) || (1.0 < abs(eta) <= 2.4 && pt > 1.2)) && innerTrack.hitPattern.trackerLayersWithMeasurement > 5 && innerTrack.hitPattern.pixelLayersWithMeasurement > 0"
process.muonSelector = cms.EDFilter("PATMuonRefSelector",
                                        src = cms.InputTag("slimmedMuons"),
                                        cut = cms.string(MUONCUT),
                                        filter = cms.bool(True)
)
process.atLeastTwoMuons = cms.EDFilter("MuonRefPatCount",
                                 src = cms.InputTag("slimmedMuons"),
                                  cut = cms.string(MUONCUT),
                                 minNumber = cms.uint32(2)
                                 )
# only process events containing at least one J/psi candidate
process.dimuonSelection = cms.EDProducer("CandViewShallowCloneCombiner",
                                    checkCharge = cms.bool(True),
                                    cut = cms.string("mass > 2.8 && mass < 3.4"),
                                    decay = cms.string("muonSelector@+ muonSelector@-")
                                    )
process.atLeastOneDimuon = cms.EDFilter("CandViewCountFilter",
                                        src = cms.InputTag("dimuonSelection"),
                                        minNumber = cms.uint32(1)
                                        )
process.p.replace(process.BfinderSequence, process.muonSelector * process.atLeastTwoMuons * process.dimuonSelection * process.atLeastOneDimuon * process.BfinderSequence)
#######################################################################################################################
#######################################################################################################################



###############################
import FWCore.ParameterSet.VarParsing as VarParsing
ivars = VarParsing.VarParsing('analysis')

ivars.maxEvents = 2000
ivars.outputFile='HiForestMINIAOD.root'
ivars.inputFiles='root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v2/000/375/790/00000/56ad580f-b228-4f3c-b8e3-17f9d95c7654.root'
ivars.parseArguments() # get and parse the command line arguments

import os
if os.path.exists(ivars.outputFile):
    os.remove(ivars.outputFile)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(ivars.inputFiles),
    # eventsToProcess = cms.untracked.VEventRange('1:236:29748033-1:236:29748033')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(ivars.maxEvents)
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(ivars.outputFile))