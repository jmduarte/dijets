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
parser.add_option("--jetNum", dest="jetNum", default = 0,help="mass of LSP", metavar="MLSP")
parser.add_option("--ZPrimeMass", dest="ZPrimeMass", default = 50., type=float,help="mass of LSP", metavar="MLSP")

# parser.add_option("--rholo", dest="rholo", default = 0.,help="mass of LSP", metavar="MLSP")
(options, args) = parser.parse_args()

from MCContainer import *
from rhalphabet import *

def buildDatacards(bkgContainers,sigContainer,theRhalphabet,theData,jetNum):

	# Create a card for each signal mass point
	# for i in range(len(sigContainers)):
	includeTop = False;

	tag = sigContainer._tag;
	basename = "combine_"+tag+"_"+str(jetNum)+".root";
	fname = "plots"+str(jetNum)+"/datacards/"+basename;
	fout = ROOT.TFile(fname,"RECREATE");

	theSignalName = "sig";
	theQCDName = "qcd";
	theWincName = "Winc";
	theZincName = "Zinc";
	theTopName  = "top";

	theSignalShape       = sigContainer.h_peakshape; ### note this new name!!!
	theWincShape         = bkgContainers[1].h_peakshape;
	theZincShape         = bkgContainers[2].h_peakshape;
	theTopShape          = bkgContainers[3].h_jetmsd_passcut;

	# theQCDShape          = bkgContainers[0].h_jetmsd_passcut; # this would eventually become the rhalphabet piece		
	theQCDShape = theRhalphabet.hpred_jetmsd
	checkTheQCDShapeForZeroBins( theQCDShape );

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
	if includeTop: theTopShape.SetName(theTopName);
	theDataObs.SetName("data_obs")
	theSignalShape.Write();
	theQCDShape.Write();
	theWincShape.Write();
	theZincShape.Write();
	if includeTop: theTopShape.Write();
	theDataObs.Write();

	# peak shifts
	theSignalShape_shiftUp       = sigContainer.h_peakshape_shiftUp;
	theWincShape_shiftUp         = bkgContainers[1].h_peakshape_shiftUp;
	theZincShape_shiftUp         = bkgContainers[2].h_peakshape_shiftUp;		
	theSignalShape_shiftDn       = sigContainer.h_peakshape_shiftDn;
	theWincShape_shiftDn         = bkgContainers[1].h_peakshape_shiftDn;
	theZincShape_shiftDn         = bkgContainers[2].h_peakshape_shiftDn;		
	theSignalShape_smearUp       = sigContainer.h_peakshape_smearUp;
	theWincShape_smearUp         = bkgContainers[1].h_peakshape_smearUp;
	theZincShape_smearUp         = bkgContainers[2].h_peakshape_smearUp;		
	theSignalShape_smearDn       = sigContainer.h_peakshape_smearDn;
	theWincShape_smearDn         = bkgContainers[1].h_peakshape_smearDn;
	theZincShape_smearDn         = bkgContainers[2].h_peakshape_smearDn;		

	theSignalShape_shiftUp.SetName(theSignalName+"_shiftUp")
	theSignalShape_shiftDn.SetName(theSignalName+"_shiftDown")
	theSignalShape_smearUp.SetName(theSignalName+"_smearUp")
	theSignalShape_smearDn.SetName(theSignalName+"_smearDown")
	theSignalShape_shiftUp.Write();
	theSignalShape_shiftDn.Write();
	theSignalShape_smearUp.Write();
	theSignalShape_smearDn.Write();

	theWincShape_shiftUp.SetName(theWincName+"_shiftUp")
	theWincShape_shiftDn.SetName(theWincName+"_shiftDown")
	theWincShape_smearUp.SetName(theWincName+"_smearUp")
	theWincShape_smearDn.SetName(theWincName+"_smearDown")
	theWincShape_shiftUp.Write();
	theWincShape_shiftDn.Write();
	theWincShape_smearUp.Write();
	theWincShape_smearDn.Write();

	theZincShape_shiftUp.SetName(theZincName+"_shiftUp")
	theZincShape_shiftDn.SetName(theZincName+"_shiftDown")
	theZincShape_smearUp.SetName(theZincName+"_smearUp")
	theZincShape_smearDn.SetName(theZincName+"_smearDown")		
	theZincShape_shiftUp.Write();
	theZincShape_shiftDn.Write();
	theZincShape_smearUp.Write();
	theZincShape_smearDn.Write();		

	# qcd closure
	qcdHist1ClosureUp = theQCDShape.Clone();
	qcdHist1ClosureDn = theQCDShape.Clone();
	qcdHist1ClosureUp.SetName("qcd_qcdClosure1Up")
	qcdHist1ClosureDn.SetName("qcd_qcdClosure1Down")
	qcdHist2ClosureUp = theQCDShape.Clone();
	qcdHist2ClosureDn = theQCDShape.Clone();
	qcdHist2ClosureUp.SetName("qcd_qcdClosure2Up")
	qcdHist2ClosureDn.SetName("qcd_qcdClosure2Down")
	qcdHist3ClosureUp = theQCDShape.Clone();
	qcdHist3ClosureDn = theQCDShape.Clone();
	qcdHist3ClosureUp.SetName("qcd_qcdClosure3Up")
	qcdHist3ClosureDn.SetName("qcd_qcdClosure3Down")	
	
	bound1 = 10; 
	bound2 = 30;
	if theQCDShape.GetNbinsX() == 30: 
		bound1 = 5; 
		bound2 = 15;

	for i in range(theQCDShape.GetNbinsX()):
		
		closurefactor = 0.01 + 0.005*float(i);

		if i < bound1:
			qcdHist1ClosureUp.SetBinContent(i+1, (1+closurefactor)*theQCDShape.GetBinContent(i+1));
			qcdHist1ClosureDn.SetBinContent(i+1, (1-closurefactor)*theQCDShape.GetBinContent(i+1));

		if i >= bound1 and i < bound2:
			qcdHist2ClosureUp.SetBinContent(i+1, (1+0.05)*theQCDShape.GetBinContent(i+1));
			qcdHist2ClosureDn.SetBinContent(i+1, (1-0.05)*theQCDShape.GetBinContent(i+1));

		if i >= bound2:
			qcdHist3ClosureUp.SetBinContent(i+1, 1.15*theQCDShape.GetBinContent(i+1));
			qcdHist3ClosureDn.SetBinContent(i+1, 0.85*theQCDShape.GetBinContent(i+1));
	
	lBFile = ROOT.TFile("qcdhists.root")
	lBQCD  = lBFile.Get("h_jetmsd_passcutQCD")
	lHist  = lBFile.Get("qcd")
	lBQCD.Divide(lHist)
	lUp3   = theQCDShape.Clone("qcd_mcclosureUp")
	lDown3 = theQCDShape.Clone("qcd_mcclosureDown")
	for i0 in range(1,lUp3.GetNbinsX()):
		if lBQCD.GetBinContent(i0) > 0:
			lUp3  .SetBinContent(i0,lUp3  .GetBinContent(i0)*lBQCD.GetBinContent(i0))
			lDown3.SetBinContent(i0,lDown3.GetBinContent(i0)*1./(lBQCD.GetBinContent(i0)))

	fout.cd();

	qcdHist1ClosureUp.Write();
	qcdHist1ClosureDn.Write();
	qcdHist2ClosureUp.Write();
	qcdHist2ClosureDn.Write();
	qcdHist3ClosureUp.Write();
	qcdHist3ClosureDn.Write();
	lUp3.Write();
	lDown3.Write();

	###############################################################
	# write the card
	ofile = open("plots"+str(jetNum)+"/datacards/combine_"+tag+"_"+str(jetNum)+".dat",'w');

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

	if includeTop: 
		line = "bin ch0 ch0 ch0 ch0 ch0 \n"
		allLines.append(line);
		line = "process 0 1 2 3 4 \n";
		allLines.append(line);
		line = "process %s %s %s %s %s \n" % ( theSignalName, theQCDName, theWincName, theZincName, theTopName );
		allLines.append(line);
		line = "rate %s %s %s %s %s \n" % ( str(round(theSignalShape.Integral(),4)), str(round(theQCDShape.Integral(),4)), str(round(theWincShape.Integral(),4)), str(round(theZincShape.Integral(),4)), str(round(theTopShape.Integral(),4)) )
		allLines.append(line);		
		allLines.append("------------ \n");
		
		line = "lumi_13TEV lnN 1.027 - - - 1.027 \n"
		allLines.append(line);		
		# line = "bkgd_QCDMClosure lnN - 1.15 - - - \n"
		# allLines.append(line);		
		line = "bkgd_VincNorm lnN 1.12 - 1.12 1.12 - \n"
		allLines.append(line);				
		line = "bkgd_TopNorm lnN - - - - 1.3 \n"
		allLines.append(line);				

		line = "shift shapeN2 1 - 1 1 - \n"
		allLines.append(line);		
		line = "smear shapeN2 1 - 1 1 - \n"
		allLines.append(line);		
		line = "WtagSF lnN 1.20 - 1.20 1.20 - \n"
		allLines.append(line);			

		WriteBinByBinUncertainties(theQCDShape," - 1 - - - ",allLines,"BinByBin",True);
		WriteBinByBinUncertainties(theSignalShape," 1 - - - - ",allLines,"BinByBinSig",True);		

		line = "qcdClosure1 shapeN2 - 1 - - - \n";
		allLines.append(line);		
		line = "qcdClosure2 shapeN2 - 1 - - - \n";
		allLines.append(line);		
		line = "qcdClosure3 shapeN2 - 1 - - - \n";
		allLines.append(line);		
		line = "mcclosure shapeN2 - 1 - - - \n";
		allLines.append(line);		

	else:
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
		# line = "bkgd_QCDMClosure lnN - 1.15 - - \n"
		# allLines.append(line);		
		line = "bkgd_VincNorm lnN 1.12 - 1.12 1.12 \n"
		allLines.append(line);				

		line = "shift shapeN2 1 - 1 1 \n"
		allLines.append(line);		
		line = "smear shapeN2 1 - 1 1 \n"
		allLines.append(line);		
		line = "WtagSF lnN 1.2 - 1.2 1.2\n"
		allLines.append(line);			

		WriteBinByBinUncertainties(theQCDShape," - 1 - - ",allLines,"BinByBin",True);
		WriteBinByBinUncertainties(theSignalShape," 1 - - - ",allLines,"BinByBinSig",True);		
		
		line = "qcdClosure1 shapeN2 - 1 - - \n";
		allLines.append(line);		
		line = "qcdClosure2 shapeN2 - 1 - - \n";
		allLines.append(line);		
		line = "qcdClosure3 shapeN2 - 1 - - \n";
		allLines.append(line);		
		line = "mcclosure shapeN2 - 1 - - \n";
		allLines.append(line);		

	for l in allLines:
		ofile.write(l);
	ofile.close();

		# ###############################################################
		# ########### make a card with W as signal
		# if i == (len(sigContainers)-1):

		# 	ofileW = open("plots"+str(jetNum)+"/datacards/combine_WZsignal"+"_"+str(jetNum)+".dat",'w');

		# 	#write the damn thing
		# 	allLines = [];

		# 	line = "#card for signal = %s \n" % (tag);
		# 	allLines.append(line);

		# 	line = "imax 1 #number of channels \n";
		# 	allLines.append(line);
		# 	line = "jmax * #number of backgrounds \n";
		# 	allLines.append(line);
		# 	allLines.append("kmax * nuissance \n");
		# 	allLines.append("------------ \n");

		# 	line = "shapes * * %s $PROCESS $PROCESS_$SYSTEMATIC \n" % (basename);
		# 	allLines.append(line);
		# 	allLines.append("------------ \n");


		# 	allLines.append("bin ch0 \n");
		# 	line = "observation %s \n" % (str(round(theDataObs.Integral(),4)));
		# 	allLines.append(line);
		# 	allLines.append("------------ \n");

		# 	line = "bin ch0 ch0 ch0 \n"
		# 	allLines.append(line);
		# 	line = "process -1 0 1 \n";
		# 	allLines.append(line);
		# 	line = "process %s %s %s \n" % ( theWincName, theZincName, theQCDName );
		# 	allLines.append(line);
		# 	line = "rate %s %s %s \n" % ( str(round(theWincShape.Integral(),4)), str(round(theZincShape.Integral(),4)), str(round(theQCDShape.Integral(),4)) )
		# 	allLines.append(line);		
		# 	allLines.append("------------ \n");

		# 	line = "lumi_13TEV lnN 1.027 1.027 - \n"
		# 	allLines.append(line);		
		# 	line = "bkgd_QCDNorm lnN - - 1.15 \n"
		# 	allLines.append(line);		
		# 	line = "bkgd_WincNorm lnN 1.2 - - \n"
		# 	allLines.append(line);		
		# 	line = "bkgd_ZincNorm lnN - 1.2 - \n"
		# 	allLines.append(line);

		# 	line = "shift shapeN2 1 1 - \n"
		# 	allLines.append(line);			
		# 	line = "smear shapeN2 1 1 - \n"
		# 	allLines.append(line);			
		# 	line = "WtagSF lnN 1.2 1.2 - \n"
		# 	allLines.append(line);			

		# 	WriteBinByBinUncertainties(theQCDShape," - - 1 ",allLines,"BinByBin",False);					

		# 	for l in allLines:
		# 		ofileW.write(l);
		# 	ofileW.close();	

	fout.Close();


def WriteBinByBinUncertainties(theH, tag, allLines, sysname, write=False):

	h_unc_up = [];
	h_unc_dn = [];
	text_unc = [];
	# text_unc_dn = [];
	for i in range(theH.GetNbinsX()):

		samplename = theH.GetName();
		uncname = sysname+str(i);
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

def checkTheQCDShapeForZeroBins( h ):

	for i in range(h.GetNbinsX()):
		if h.GetBinContent(i+1) <= 0.:
			h.SetBinContent(i+1,0.01);
			h.SetBinError(i+1,1.87);	
			print "bin #",i+1," got modified!"		








