################################################################
## In principle all you have to setup is defined in this file ##
################################################################
from configManager import configMgr
from ROOT import kBlack,kWhite,kGray,kRed,kPink,kMagenta,kViolet,kBlue,kAzure,kCyan,kTeal,kGreen,kSpring,kYellow,kOrange
from configWriter import fitConfig,Measurement,Channel,Sample
from systematic import Systematic
from math import sqrt
import os

# Setup for ATLAS plotting
from ROOT import gROOT
gROOT.LoadMacro("./macros/AtlasStyle.C")
import ROOT
ROOT.SetAtlasStyle()

##########################

# Set run mode
# runMode = "exclusion"                     # run with -l
# runMode = "discovery"                   # run with -z
runMode = "model_independent_discovery" # run with -l

# Set observed and expected number of events in counting experiment
ndata = 1 # Number of events observed in data
nbkg = 1.8 # Number of predicted bkg events
nsig = 1 # Number of predicted signal events
nbkgErr = 0.54 # (Absolute) Statistical error on bkg estimate
nsigErr = 1 # (Absolute) Statistical error on signal estimate
lumiError = 0.028 # Relative luminosity uncertainty
up  = 1.3;
dwn = 0.7;
ucorb_up      = up;
ucorb_dwn     = dwn;
ucors_up      = 1.0;
ucors_dwn     = 1.0;
corb_up       = 1.0;
corb_dwn      = 1.0;
cors_up       = 1.0;
cors_dwn      = 1.0;
normFactorMin = 0.;
normFactorMax = 100.;

##########################

# Set uncorrelated systematics for bkg and signal (1 +- relative uncertainties)
ucb = Systematic("ucb", configMgr.weights, ucorb_up,ucorb_dwn, "user","userOverallSys")
ucs = Systematic("ucs", configMgr.weights, ucors_up,ucors_dwn, "user","userOverallSys")

# correlated systematic between background and signal (1 +- relative uncertainties)
corb = Systematic("cor",configMgr.weights, [corb_up],[corb_dwn], "user","userHistoSys")
cors = Systematic("cor",configMgr.weights, [cors_up],[cors_dwn], "user","userHistoSys")

##########################

# Setting the parameters of the hypothesis test
#configMgr.nTOYs=5000
configMgr.calculatorType=2 # 2=asymptotic calculator, 0=frequentist calculator
configMgr.testStatType=3   # 3=one-sided profile likelihood test statistic (LHC default)
configMgr.nPoints=20       # number of values scanned of signal-strength for upper-limit determination of signal strength.

##########################

# Give the analysis a name
configMgr.analysisName = "MyUserAnalysis"
configMgr.outputFileName = "results/%s_Output.root"%configMgr.analysisName

# Define cuts
configMgr.cutsDict["UserRegion"] = "1."

# Define weights
configMgr.weights = "1."

# Define samples
bkgSample = Sample("Bkg",kGreen-9)
bkgSample.setStatConfig(True)
bkgSample.buildHisto([nbkg],"UserRegion","cuts")
bkgSample.buildStatErrors([nbkgErr],"UserRegion","cuts") ###
if(runMode=="exclusion"):
	bkgSample.addSystematic(corb)
bkgSample.addSystematic(ucb)

sigSample = Sample("Sig",kPink)
sigSample.setNormFactor("mu_Sig",1.,normFactorMin,normFactorMax)
sigSample.setStatConfig(False)
sigSample.setNormByTheory(False)
sigSample.buildHisto([nsig],"UserRegion","cuts")
sigSample.buildStatErrors([nsigErr],"UserRegion","cuts") ###
sigSample.addSystematic(cors) ###
sigSample.addSystematic(ucs) ###

dataSample = Sample("Data",kBlack)
dataSample.setData()
dataSample.buildHisto([ndata],"UserRegion","cuts")

# Define top-level
ana = configMgr.addTopLevelXML("SPlusB")
if(runMode=="exclusion"):
	ana.addSamples([bkgSample,sigSample,dataSample])
else:
	ana.addSamples([bkgSample,dataSample])
if(runMode=="exclusion"):
	ana.setSignalSample(sigSample)

# Define measurement
meas = ana.addMeasurement(name="NormalMeasurement",lumi=1.0,lumiErr=lumiError)
meas.addPOI("mu_Sig")
meas.addParamSetting("Lumi","const",1.0)

# Add the channel
chan = ana.addChannel("cuts",["UserRegion"],1,0.,1.)
if(runMode=="discovery" or runMode=="model_independent_discovery"):
	chan.addDiscoverySamples(["Sig"],[1.],[normFactorMin],[normFactorMax],[ROOT.kPink])
ana.setSignalChannels([chan])

# These lines are needed for the user analysis to run
# Make sure file is re-made when executing HistFactory
if configMgr.executeHistFactory:
    if os.path.isfile("data/%s.root"%configMgr.analysisName):
        os.remove("data/%s.root"%configMgr.analysisName)
