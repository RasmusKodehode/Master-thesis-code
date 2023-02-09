from AnaAlgorithm.AlgSequence import AlgSequence
from AnaAlgorithm.DualUseConfig import createService

def makeSequence (dataType) :

    algSeq = AlgSequence()

    # Set up the systematics loader/handler service:
    sysService = createService( 'CP::SystematicsSvc', 'SystematicsSvc', sequence = algSeq )
    sysService.sigmaRecommended = 1

    # Include, and then set up the muon analysis algorithm sequence:
    from MuonAnalysisAlgorithms.MuonAnalysisSequence import makeMuonAnalysisSequence
    muonSequenceMedium = makeMuonAnalysisSequence( dataType, deepCopyOutput = True, shallowViewOutput = False,
                                                   workingPoint = 'Medium.NonIso', postfix = 'medium' )
    muonSequenceMedium.configure( inputName = 'Muons',
                                  outputName = 'AnalysisMuons_%SYS%' )

    # Add the sequence to the job:
    algSeq += muonSequenceMedium

    # Include, and then set up the pileup analysis sequence:
    from AsgAnalysisAlgorithms.PileupAnalysisSequence import \
        makePileupAnalysisSequence
    lumicalcfiles = []
    prwfiles = []
    pileupSequence = makePileupAnalysisSequence( dataType,
                                                 userPileupConfigs=prwfiles,
                                                 userLumicalcFiles=lumicalcfiles,
                                               )
    pileupSequence.configure( inputName = 'EventInfo', outputName = 'EventInfo_%SYS%' )

    # Add the pileup sequence to the job:
    algSeq += pileupSequence


    return algSeq