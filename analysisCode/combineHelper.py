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
parser.add_option("--rhohi", dest="rhohi", default = 6.,help="mass of LSP", metavar="MLSP")
parser.add_option("--DDTcut", dest="DDTcut", default = 0.38,help="mass of LSP", metavar="MLSP")
parser.add_option('--qcdClosure', action='store_true', dest='qcdClosure', default=False, help='go!')

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
		theZincName = "Zinc";

		theSignalShape       = sigContainers[i].h_jetmsd_passcut;
		theWincShape         = bkgContainers[1].h_jetmsd_passcut;
		theZincShape         = bkgContainers[2].h_jetmsd_passcut;

		##theQCDShape          = bkgContainers[0].h_jetmsd_passcut; # this would eventually become the rhalphabet piece		
		theQCDShape = theRhalphabet.hpred_jetmsd

		theDataObs = None;
		if options.qcdClosure:
			theDataObs = theQCDShape.Clone("data_obs"); # and this would eventually become the data :)
			theDataObs.Add(theWincShape);
			theDataObs.Add(theZincShape);
		else:
			theDataObs = theData.h_jetmsd_passcut;

		fout.cd();
		theSignalShape.SetName(theSignalName);
		theQCDShape.SetName(theQCDName);
		theWincShape.SetName(theWincName);
		theZincShape.SetName(theZincName);
		theDataObs.SetName("data_obs")
		theSignalShape.Write();
		theQCDShape.Write();
		theWincShape.Write();
		theZincShape.Write();
		theDataObs.Write();

		###############################################################
		# write the card
		ofile = open("plots/datacards/combine_"+tag+".dat",'w');

		#write the damn thing
		allLines = [];

		line = "#card for signal = %s \n" % (tag);
		allLines.append(line);

		line = "imax 1 #number of channels \n";
		allLines.append(line);
		line = "jmax * #number of backgrounds \n";
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

		line = "bin ch0 ch0 ch0 ch0 \n"
		allLines.append(line);
		line = "process 0 1 2 3 \n";
		allLines.append(line);
		line = "process %s %s %s %s \n" % ( theSignalName, theQCDName, theWincName, theZincName );
		allLines.append(line);
		line = "rate %s %s %s %s \n" % ( str(round(theSignalShape.Integral(),4)), str(round(theQCDShape.Integral(),4)), str(round(theWincShape.Integral(),4)), str(round(theZincShape.Integral(),4)) )
		allLines.append(line);		
		allLines.append("------------ \n");

		line = "lumi_13TEV lnN 1.027 - - - \n"
		allLines.append(line);		
		line = "bkgd_QCDMClosure lnN - 1.15 - - \n"
		allLines.append(line);		
		line = "bkgd_WincNorm lnN - - 1.2 - \n"
		allLines.append(line);				
		line = "bkgd_ZincNorm lnN - - - 1.2 \n"
		allLines.append(line);		

		WriteBinByBinUncertainties(theQCDShape," - 1 - - ",allLines,True);

		for l in allLines:
			ofile.write(l);
		ofile.close();

		###############################################################
		########### make a card with W as signal
		if i == 0:

			ofileW = open("plots/datacards/combine_WZsignal.dat",'w');

			#write the damn thing
			allLines = [];

			line = "#card for signal = %s \n" % (tag);
			allLines.append(line);

			line = "imax 1 #number of channels \n";
			allLines.append(line);
			line = "jmax * #number of backgrounds \n";
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
			line = "process -1 0 1 \n";
			allLines.append(line);
			line = "process %s %s %s \n" % ( theWincName, theZincName, theQCDName );
			allLines.append(line);
			line = "rate %s %s %s \n" % ( str(round(theWincShape.Integral(),4)), str(round(theZincShape.Integral(),4)), str(round(theQCDShape.Integral(),4)) )
			allLines.append(line);		
			allLines.append("------------ \n");

			line = "lumi_13TEV lnN 1.027 1.027 - \n"
			allLines.append(line);		
			line = "bkgd_QCDNorm lnN - - 1.15 \n"
			allLines.append(line);		
			line = "bkgd_WincNorm lnN 1.2 - - \n"
			allLines.append(line);		
			line = "bkgd_ZincNorm lnN - 1.2 - \n"
			allLines.append(line);

			WriteBinByBinUncertainties(theQCDShape," - - 1 ",allLines,False);					

			for l in allLines:
				ofileW.write(l);
			ofileW.close();	

		fout.Close();


def WriteBinByBinUncertainties(theH, tag, allLines, write=False):

	h_unc_up = [];
	h_unc_dn = [];
	text_unc = [];
	# text_unc_dn = [];
	for i in range(theH.GetNbinsX()):

		samplename = theH.GetName();
		uncname = "BinByBin"+str(i);
		sysnameUp = samplename + "_" + uncname + "Up";
		sysnameDn = samplename + "_" + uncname + "Down";

		text_unc.append( uncname + " shapeN2 " + tag + " \n" );

		h_unc_up.append(theH.Clone());
		h_unc_dn.append(theH.Clone());
		curbincontent = theH.GetBinContent(i+1);
		curbinerror = theH.GetBinError(i+1);
		
		h_unc_up[i].SetBinContent(i+1, curbincontent+curbinerror);
		h_unc_dn[i].SetBinContent(i+1, curbincontent-curbinerror);
		h_unc_up[i].SetName(sysnameUp);
		h_unc_dn[i].SetName(sysnameDn);

	if write:
		for i in range(len(text_unc)):
			# print text_unc_up[i], text_unc_dn[i], h_unc_up[i].GetName(), h_unc_dn[i].GetName();
			h_unc_up[i].Write();
			h_unc_dn[i].Write();

	for i in range(len(text_unc)):		
		allLines.append( text_unc[i] );








