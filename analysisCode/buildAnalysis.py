import ROOT
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--doMCLooping', action='store_true', dest='doMCLooping', default=False, help='go!')
parser.add_option('--doRhalphabet', action='store_true', dest='doRhalphabet', default=False, help='go!')
parser.add_option('--doData', action='store_true', dest='doData', default=False, help='go!')
parser.add_option('--doPlots', action='store_true', dest='doPlots', default=False, help='go!')
parser.add_option("--rholo", dest="rholo", default = 0.,help="mass of LSP", metavar="MLSP")
parser.add_option("--rhohi", dest="rhohi", default = 4.,help="mass of LSP", metavar="MLSP")
parser.add_option("--DDTcut", dest="DDTcut", default = 0.45,help="mass of LSP", metavar="MLSP")

# parser.add_option("--rholo", dest="rholo", default = 0.,help="mass of LSP", metavar="MLSP")
(options, args) = parser.parse_args()

from MCContainer import *
from rhalphabet import *
from plotHelpers import makeCanvas, makeCanvasDataMC, makeCanvasShapeComparison,makeCanvas2D,makeCanvasDataMC_wpred

#
import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadTopMargin(0.06);
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPadRightMargin(0.10);
ROOT.gStyle.SetPalette(1);
ROOT.gStyle.SetPaintTextFormat("1.1f");
ROOT.gStyle.SetOptFit(0000);

###############################################################################################################
# M A I N 
###############################################################################################################
def main(): 

	idir = "/Users/ntran/Documents/Research/CMS/WZpToQQ/bkgEst/sklim-v2";

	####################################################################################
	# do mc looping - a class that holds histograms
	if options.doMCLooping: 
		
		bkgContainers = [];
		bkgNames = ["QCD.root","W.root"];
		bkgLabels = ["QCD","W(qq)"];
		for i in range(len(bkgNames)):
			bkgContainers.append( MCContainer( idir+"/"+bkgNames[i], 0.44, bkgLabels[i], 10 ) );

		sigContainers = [];
		sigNames = [];
		sigNames.append("ZPrimeToQQ_50GeV_v4_mc.root")
		sigNames.append("ZPrimeToQQ_100GeV_v4_mc.root")
		sigNames.append("ZPrimeToQQ_150GeV_v4_mc.root")
		sigNames.append("ZPrimeToQQ_200GeV_v4_mc.root")
		sigNames.append("ZPrimeToQQ_250GeV_v4_mc.root")
		sigNames.append("ZPrimeToQQ_300GeV_v4_mc.root")	
		sigLabels = ["Z\'(50 GeV)","Z\'(100 GeV)","Z\'(150 GeV)","Z\'(200 GeV)","Z\'(250 GeV)","Z\'(300 GeV)"]
		for i in range(len(bkgNames)):
			sigContainers.append( MCContainer( idir+"/"+sigNames[i], 0.44, sigLabels[i], 1 ) );

	####################################################################################
	# do background estimation
	if options.doRhalphabet: 
		# isData = True;
		# theRhalphabet = rhalphabet(idir+"/"+"JetHT.root",1,"rhalphabet",1, True);
		# theRhalphabet.GetPredictedDistributions( idir+"/"+"JetHT.root", 1, 5, isData);
		
		# there is a flag to do a closure test as well
		isData = False;
		theRhalphabet = rhalphabet(idir+"/"+"QCD.root",1,"rhalphabetClosure",1, False);
		theRhalphabet.GetPredictedDistributions( idir+"/"+"QCD.root", 0.44, 10, isData );

	####################################################################################
	# do the loop on data
	theData = None;
	if options.doData:
		isData = True;
		theData = MCContainer( idir+"/"+"JetHT.root", 1, "data", 5, isData);

	####################################################################################
	# do some plotting
	if options.doPlots: 
		BuildPlots(bkgContainers,sigContainers,theRhalphabet,theData);

	# c = ROOT.TCanvas('c','c',1000,800);
	# bkgContainers[0].h_jetpt.Draw();
	# c.SaveAs("plots/test.pdf");

def BuildPlots(bkgContainers,sigContainers,theRhalphabet,theData):

	print "making plots...";

	# makeCanvasDataMC_wpred( theData.h_jetmsd_passcut,
	# 						theRhalphabet.grpred_jetmsd, 
	# 						[bkgContainers[0].h_jetmsd_passcut],
	# 						['qcd'],
	# 						'jetmsd_pred',
	# 						'plots/results/',
	# 						False);

	makeCanvasDataMC_wpred( bkgContainers[0].h_jetmsd_passcut,
							theRhalphabet.grpred_jetmsd, 
							[bkgContainers[0].h_jetmsd_passcut],
							['qcd'],
							'jetmsd_pred',
							'plots/results/',
							False);

	makeCanvasDataMC_wpred( bkgContainers[0].h_rhoDDT_passcut,
							theRhalphabet.grpred_rhoDDT, 
							[bkgContainers[0].h_rhoDDT_passcut],
							['qcd'],
							'rhoDDT_pred',
							'plots/results/',
							False);



#----------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	main();
#----------------------------------------------------------------------------------------------------------------
