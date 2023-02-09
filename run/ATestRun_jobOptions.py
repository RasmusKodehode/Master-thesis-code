#See: https://twiki.cern.ch/twiki/bin/viewauth/AtlasComputing/SoftwareTutorialxAODAnalysisInCMake for more details about anything here

testFile = "mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.AOD.e3601_e5984_s3126_s3136_r9364_r9315"

#override next line on command line with: --filesInput=XXX
jps.AthenaCommonFlags.FilesInput = [testFile] 

#Specify AccessMode (read mode) ... ClassAccess is good default for xAOD
jps.AthenaCommonFlags.AccessMode = "ClassAccess" 

jps.AthenaCommonFlags.HistOutputs = ["ANALYSIS:MyxAODAnalysis.outputs.root"]
svcMgr.THistSvc.MaxFileSize=-1 #speeds up jobs that output lots of histograms

# Create the algorithm's configuration.
#from AnaAlgorithm.DualUseConfig import createAlgorithm
#alg = createAlgorithm ( 'MyxAODAnalysis', 'AnalysisAlg' )
alg = CfgMgr.MyxAODAnalysis("AnalysisAlg")

# later on we'll add some configuration options for our algorithm that go here
from AnaAlgorithm.DualUseConfig import addPrivateTool

# add the GRL tool to the algorithm
addPrivateTool( alg, 'grlTool', 'GoodRunsListSelectionTool' )

# configure the properties of the GRL tool
fullGRLFilePath = os.getenv ("ALRB_TutorialData") + "/data16_13TeV.periodAllYear_DetStatus-v89-pro21-01_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml"
alg.grlTool.GoodRunsListVec = [ fullGRLFilePath ]
alg.grlTool.PassThrough = 0 # if true (default) will ignore result of GRL and will just pass all events

inputFilePath = os.getenv( 'ALRB_TutorialData' ) + '/mc21_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_PHYS.e8357_s3802_r13508_p5057/'
ROOT.SH.ScanDir().filePattern( 'DAOD_PHYS.28625583._000007.pool.root.1' ).scan( sh, inputFilePath )
sh.printContent()

dataType = "mc" 

# Add our algorithm to the main alg sequence
athAlgSeq += alg

# limit the number of events (for testing purposes)
theApp.EvtMax = -1

# optional include for reducing printout from athena
include("AthAnalysisBaseComps/SuppressLogging.py")