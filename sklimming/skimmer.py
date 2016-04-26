#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time

import ROOT

# ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
# ROOT.setTDRStyle();
# ROOT.gStyle.SetPadTopMargin(0.06);
# ROOT.gStyle.SetPadLeftMargin(0.16);
# ROOT.gStyle.SetPadRightMargin(0.10);
# ROOT.gStyle.SetPalette(1);
# ROOT.gStyle.SetPaintTextFormat("1.1f");

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--train', action='store_true', dest='train', default=False, help='no X11 windows')

(options, args) = parser.parse_args()

############################################################

# observableTraining takes in 2 root files, list of observables, spectator observables ... launches a CONDOR job
# TMVAhelper.py is used by observableTraining
# analysis.py defines the list of trainings, the signal and background process

########################################################################################################################
########################################################################################################################

def getSampleXS(name):

	QCDFudgeFactor = 1./7.;

	xs = -999;
	if 'ttbar' in name and '0_600' in name: xs = 435.841349313;
	elif 'ttbar' in name and '600_1200' in name: xs = 32.6826950495;
	elif 'ttbar' in name and '1200_1900' in name: xs = 2.0582766;
	elif 'ttbar' in name and '1900_2700' in name: xs = 0.168654823762;
	elif 'ttbar' in name and '2700_3600' in name: xs = 0.01588917;
	elif 'ttbar' in name and '3600_4600' in name: xs = 0.0015207095742;
	elif 'ttbar' in name and '4600_100000' in name: xs = 0.000141367514851;

	elif 'Wjets' in name and '0_400' in name: xs = 6363.608615;
	elif 'Wjets' in name and '400_900' in name: xs = 432.88054467;
	elif 'Wjets' in name and '900_1600' in name: xs = 23.3849837436;
	elif 'Wjets' in name and '1600_2500' in name: xs = 1.65095612398;
	elif 'Wjets' in name and '2500_3500' in name: xs = 0.129479800515;
	elif 'Wjets' in name and '3500_4600' in name: xs = 0.0121681361371;
	elif 'Wjets' in name and '4600_5700' in name: xs = 0.001086189415;
	elif 'Wjets' in name and '5700_10000' in name: xs = 0.000106971454914;

	# elif 'QCD' in name and '0_400' in name: xs = 1.18e6;
	# elif 'QCD' in name and '400_800' in name: xs = 6.09e4;
	# elif 'QCD' in name and '800_1300' in name: xs = 2550.0;
	# elif 'QCD' in name and '1300_2000' in name: xs = 196.0;
	# elif 'QCD' in name and '2000_2900' in name: xs = 14.7;
	# elif 'QCD' in name and '2900_3900' in name: xs = 1.06;
	# elif 'QCD' in name and '3900_4900' in name: xs = 8.71e-2;
	# elif 'QCD' in name and '4900_5900' in name: xs = 8.51e-3;
	# elif 'QCD' in name and '5900_6900' in name: xs = 8.47e-4;
	# elif 'QCD' in name and '6900_100000' in name: xs = 8.46e-5;
	
	# elif 'QCD' in name and '0_300' in name: xs = 10289941.6637;
	elif 'QCD' in name and '300_600' in name: xs = 266000.0;
	elif 'QCD' in name and '600_1000' in name: xs = 12546.4206;
	elif 'QCD' in name and '1000_1600' in name: xs = 972.061327;
	elif 'QCD' in name and '1600_2400' in name: xs = 68.2020051;
	elif 'QCD' in name and '2400_3300' in name: xs = 4.84015207;
	elif 'QCD' in name and '3300_4300' in name: xs = 0.429700105;
	elif 'QCD' in name and '4300_5300' in name: xs = 0.0386151612871;
	elif 'QCD' in name and '5300_6300' in name: xs = 0.0038;
	elif 'QCD' in name and '6300_100000' in name: xs = 0.000413720491984;

	# elif 'znunu' in name and '0_400' in name: xs = 136.6;
	# elif 'znunu' in name and '400_900' in name: xs = 9.37;
	# elif 'znunu' in name and '900_1600' in name: xs = 0.59;
	# elif 'znunu' in name and '1600_2500' in name: xs = 0.044;
	# elif 'znunu' in name and '2500_3500' in name: xs = 3.53e-3;
	# elif 'znunu' in name and '3500_100000' in name: xs = 3.50e-4;
	elif 'znunu' in name and '0_600' in name: xs = 539.8;
	elif 'znunu' in name and '600_1200' in name: xs = 26.72;
	elif 'znunu' in name and '1200_2000' in name: xs = 1.926;
	elif 'znunu' in name and '2000_2900' in name: xs = 0.151;
	elif 'znunu' in name and '2900_3900' in name: xs = 1.51e-2;
	elif 'znunu' in name and '3900_5000' in name: xs = 1.56e-3;
	elif 'znunu' in name and '5000_100000' in name: xs = 1.53e-4;

	elif 'Gj' in name: xs = -999;

	else: 
		raise NameError("too many keys!")
		return;		

	return xs;

def getLHEWeight(files):
	
	print "processing: ",files[0]
	theXS = getSampleXS(files[0])
	totalEvents = 0;

	for f in files:
		f1 = ROOT.TFile(f,'read');
		thekey = None;
		if len(ROOT.gDirectory.GetListOfKeys()) > 1: 
			raise NameError("too many keys!")
			return;
		for key in ROOT.gDirectory.GetListOfKeys(): thekey = key;
		tree = thekey.ReadObj();
		totalEvents += tree.GetEntriesFast();

	print totalEvents, theXS;
	if theXS == -999: return 1;
	else: return (theXS/totalEvents)

def sklimAdd(fn,odir):

	basename = os.path.basename( fn );

	f1 = ROOT.TFile(fn,'read');
	tree = f1.Get("Events");

	ofile = ROOT.TFile(odir+'/'+basename,'RECREATE');
	ofile.cd();
	otree = tree.CloneTree(0);
	otree.SetName("otree");
	otree.SetBranchStatus("*res_*",0);
	# otree.SetBranchStatus("bst8_PUPPIjet0_pt",1);

	# lheWeight = array( 'f', [ 0. ] );  
	# MHTOvHT = array( 'f', [ 0. ] );  
	# b_lheWeight = otree.Branch("lheWeight",lheWeight,"lheWeight/F"); #xs*lumi/nev, take lumi as 1fb-1
	# b_MHTOvHT = otree.Branch("MHTOvHT",MHTOvHT,"MHTOvHT/F"); #xs*lumi/nev, take lumi as 1fb-1

	nent = tree.GetEntriesFast();
	for i in range(nent):

		if(i % (1 * nent/100) == 0):
			sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
			sys.stdout.flush();

		tree.GetEntry(i);
		# print tree.HT, tree.mT2, tree.alphaT, tree.dRazor, tree.mRazor, tree.sumJetMass

		if tree.bst8_PUPPIjet0_pt > 200:
			# throw out NaN values...
			# print tree.HT, tree.mT2, tree.alphaT, tree.dRazor, tree.mRazor, tree.sumJetMass

			# curalphaT = tree.alphaT;
			# curdRazor = tree.dRazor;
			# if math.isnan(curalphaT) or math.isnan(curdRazor): continue;
			# # print tree.HT, tree.mT2, tree.alphaT, tree.dRazor, tree.mRazor, tree.sumJetMass
			# if int(tree.HT) % 2 == modval: continue;
			# # print int(tree.HT)
			# lheWeight[0] = float(weight);
			# MHTOvHT[0] = tree.MHT/math.sqrt(tree.HT);
			otree.Fill();   

	print "\n"
	#otree.Print();
	otree.AutoSave();
	ofile.cd();
	otree.Write();
	ofile.Close();

def getFilesRecursively(dir,searchstring,additionalstring = None):
	
	# thesearchstring = "_"+searchstring+"_";
	thesearchstring = searchstring;
	
	theadditionalstring = None;
	if not additionalstring == None: 
		theadditionalstring = additionalstring;

	cfiles = [];
	for root, dirs, files in os.walk(dir):
		for file in files:
			# print file	
			if thesearchstring in file:
				if theadditionalstring == None or theadditionalstring in file:
					cfiles.append(os.path.join(root, file))
	return cfiles;

if __name__ == '__main__':

	DataDir = 'zprimebits_notrigger/'
	# DataDir = "/Users/ntran/Documents/Research/Ext/DissectingJetsPlusMET/sampleProcessing/DissectingJetsPlusMET/localData/Backgrounds/Backgrounds_13TEV/TTBAR/";
	# DataDir = "/Users/ntran/Documents/Research/Ext/DissectingJetsPlusMET/sampleProcessing/DissectingJetsPlusMET/andrewBkg/";
	# DataDir = "/Users/ntran/Documents/Research/Ext/DissectingJetsPlusMET/sampleProcessing/DissectingJetsPlusMET/rawData-v3/";
	OutDir = 'sklim-v0'

	tags = [];
	# tags.append( ['QCD'] );
	# tags.append( ['W.ro'] );
	# tags.append( ['ZPrimeToQQ_50GeV_v4_mc'] );
	# tags.append( ['ZPrimeToQQ_100GeV_v4_mc'] );
	# tags.append( ['ZPrimeToQQ_150GeV_v4_mc'] );
	# tags.append( ['ZPrimeToQQ_200GeV_v4_mc'] );
	# tags.append( ['ZPrimeToQQ_250GeV_v4_mc'] );
	# tags.append( ['ZPrimeToQQ_300GeV_v4_mc'] );
	tags.append( ['JetHT'] );


	# make a tmp dir
	#####
	postfix = '';
	for i in range(len(tags)):
		
		filesToConvert = getFilesRecursively(DataDir,tags[i][0]);
		print "files To Convert = ",filesToConvert
		# curweight = getLHEWeight( filesToConvert );
		

		for f in filesToConvert:
			print f;
			sklimAdd(f,OutDir);
		## hadd stuff
	# 	oname = OutDir + '/ProcJPM_'+tags[i][0]+"_"+tags[i][1]+"-"+postfix+".root";
	# 	# oname = OutDir + '/ProcJPM_'+tags[i][0]+"_"+tags[i][1]+".root";
	# 	haddCmmd = 'hadd -f '+oname+" ";
	# 	for f in filesToConvert:
	# 		ofile = OutDir + "/" + os.path.basename( f );
	# 		haddCmmd += ofile + " ";
	# 	print haddCmmd;
	# 	os.system(haddCmmd);

	# for files in os.listdir(OutDir):
	# 	if 'slim' in files: os.system('rm '+OutDir+'/'+files)

	# cmmd = 'hadd -f %s/ProcJPM_ttbar-%s.root %s/*ttbar*-%s.root' % (OutDir,postfix,OutDir,postfix)
	# os.system(cmmd)
	# cmmd = 'hadd -f %s/ProcJPM_Wjets-%s.root %s/*Wjets*-%s.root' % (OutDir,postfix,OutDir,postfix)
	# os.system(cmmd)
	# cmmd = 'hadd -f %s/ProcJPM_QCD-%s.root %s/*QCD*-%s.root' % (OutDir,postfix,OutDir,postfix)
	# os.system(cmmd)
	# cmmd = 'hadd -f %s/ProcJPM_znunu-%s.root %s/*znunu*-%s.root' % (OutDir,postfix,OutDir,postfix)
	# os.system(cmmd)



