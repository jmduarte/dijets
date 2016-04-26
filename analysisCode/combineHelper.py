import ROOT
from ROOT import TFile, TTree, TChain, gPad, gDirectory, TVirtualFitter
import math
import sys
import time
import array
from optparse import OptionParser
from operator import add

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--doMCLooping', action='store_true', dest='doMCLooping', default=False, help='go!')
parser.add_option('--doRhalphabet', action='store_true', dest='doRhalphabet', default=False, help='go!')
parser.add_option('--doData', action='store_true', dest='doData', default=False, help='go!')
parser.add_option('--doPlots', action='store_true', dest='doPlots', default=False, help='go!')
parser.add_option('--doCards', action='store_true', dest='doCards', default=False, help='go!')

parser.add_option("--lumi", dest="lumi", default = 0.44,help="mass of LSP", metavar="MLSP")
parser.add_option("--rholo", dest="rholo", default = 0.,help="mass of LSP", metavar="MLSP")
parser.add_option("--rhohi", dest="rhohi", default = 4.,help="mass of LSP", metavar="MLSP")
parser.add_option("--DDTcut", dest="DDTcut", default = 0.45,help="mass of LSP", metavar="MLSP")

# parser.add_option("--rholo", dest="rholo", default = 0.,help="mass of LSP", metavar="MLSP")
(options, args) = parser.parse_args()

from MCContainer import *
from rhalphabet import *

def buildDatacards(bkgContainers,sigContainers,theRhalphabet,theData):

	# Create a card for each signal mass point
	for i in range(len(sigContainers)):

		tag = sigContainers[i]._tag;
		basename = "combine_"+tag+".root";
		fname = "plots/datacards/"+basename;
		fout = ROOT.TFile(fname,"RECREATE");

		theSignalName = "sig";
		theQCDName = "qcd";
		theWincName = "Winc";

		theSignalShape       = sigContainers[i].h_jetmsd_passcut;
		theQCDShape          = bkgContainers[0].h_jetmsd_passcut; # this would eventually become the rhalphabet piece
		theWincShape         = bkgContainers[1].h_jetmsd_passcut;
		theDataObs           = theQCDShape.Clone("data_obs"); # and this would eventually become the data :)
		theDataObs.Add(theWincShape);

		fout.cd();
		theSignalShape.SetName(theSignalName);
		theQCDShape.SetName(theQCDName);
		theWincShape.SetName(theWincName);
		theDataObs.SetName("data_obs")
		theSignalShape.Write();
		theQCDShape.Write();
		theWincShape.Write();
		theDataObs.Write();

		ofile = open("plots/datacards/combine_"+tag+".dat",'w');

		#write the damn thing
		allLines = [];

		line = "#card for signal = %s \n" % (tag);
		allLines.append(line);

		line = "imax 1 #number of channels \n";
		allLines.append(line);
		line = "jmax 2 #number of backgrounds \n";
		allLines.append(line);
		allLines.append("kmax * nuissance \n");
		allLines.append("------------ \n");

		line = "shapes * * %s $PROCESS $PROCESS_$SYSTEMATIC \n" % (basename);
		allLines.append(line);
		allLines.append("------------ \n");


		allLines.append("bin ch0 \n");
		line = "observation %s \n" % (str(round(theDataObs.Integral(),4)));
		allLines.append(line);
		allLines.append("------------ \n");

		line = "bin ch0 ch0 ch0 \n"
		allLines.append(line);
		line = "process 0 1 2 \n";
		allLines.append(line);
		line = "process %s %s %s \n" % ( theSignalName, theQCDName, theWincName );
		allLines.append(line);
		line = "rate %s %s %s \n" % ( str(round(theSignalShape.Integral(),4)), str(round(theQCDShape.Integral(),4)), str(round(theWincShape.Integral(),4)) )
		allLines.append(line);		
		allLines.append("------------ \n");

		line = "lumi_13TEV lnN 1.027 - - \n"
		allLines.append(line);		
		line = "bkgd_QCDNorm lnN - 1.3 - \n"
		allLines.append(line);		
		line = "bkgd_WincNorm lnN - - 1.2 \n"
		allLines.append(line);		

		for l in allLines:
			ofile.write(l);
		ofile.close();

		fout.Close();
















