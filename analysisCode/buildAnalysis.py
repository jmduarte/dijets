import ROOT
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array

from MCContainer import *
from ralphabet import *

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
# parser.add_option("--rholo", dest="rholo", default = 0.,help="mass of LSP", metavar="MLSP")
(options, args) = parser.parse_args()

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
if __name__ == '__main__':

	idir = "/Users/ntran/Documents/Research/CMS/WZpToQQ/bkgEst/sklim-v2";
	doMCLooping = False

	####################################################################################
	# do mc looping - a class that holds histograms
	if doMCLooping: 
		
		bkgContainers = [];
		bkgNames = ["QCD.root","W.root"];
		bkgLabels = ["QCD","W(qq)"];
		for i in range(len(bkgNames)):
			bkgContainers.append( MCContainer( idir+"/"+bkgNames[i], 0.44, bkgLabels[i], 100 ) );

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
	theRhalphabet = rhalphabet(idir+"/"+"JetHT.root",1,"rhalphabet",100);
	# there is a flag to do a closure test as well

	####################################################################################
	# do the loop on data
	theData = MCContainer( idir+"/"+"JetHT.root", 1, "data", 100 );

	####################################################################################
	# do some plotting
	# BuildPlots(bkgContainers,sigContainers,theRhalphabet,theData);

	# c = ROOT.TCanvas('c','c',1000,800);
	# bkgContainers[0].h_jetpt.Draw();
	# c.SaveAs("plots/test.pdf");
