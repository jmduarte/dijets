import ROOT
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array
import os

sys.path.append('../fitting')
from hist import hist

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--doMCLooping', action='store_true', dest='doMCLooping', default=True, help='go!')
parser.add_option('--doRhalphabet', action='store_true', dest='doRhalphabet', default=False, help='go!')
parser.add_option('--doData', action='store_true', dest='doData', default=True, help='go!')
parser.add_option('--doPlots', action='store_true', dest='doPlots', default=True, help='go!')
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
from plotHelpers import makeCanvas, makeCanvasDataMC, makeCanvasShapeComparison,makeCanvas2D,makeCanvasDataMC_wpred,makeCanvasDataMC_MONEY,makeCanvasComparison,makeROCFromHisto,plotROCs
from combineHelper import buildDatacards
#
import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadTopMargin(0.06);
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPadRightMargin(0.10);
ROOT.gStyle.SetPalette(1);
ROOT.gStyle.SetPaintTextFormat("1.1f");
ROOT.gStyle.SetOptFit(0000);

"""
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 50. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 60. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 75. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 90. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 100. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 110. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 125. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 135. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 150. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 165. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 180. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 200. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 250. &
python buildAnalysis.py -b --doMCLooping --doRhalphabet --doCards --doPlots --doData --lumi 2.7 --DDTcut 0.38 --rholo 0.0 --rhohi 6. --jetNum 0 --ZPrimeMass 300. &
"""


###############################################################################################################
# M A I N 
###############################################################################################################

sigmass = [50.,75.,100.,125.,150.,200.,250.,300]
NMassBins = [60,60,60,60,60,60,30,30];
# sigmass = [50.,75.,100.,125.,150.,200.]
# NMassBins = 60;
# sigmass = [250.,300.]
# NMassBins = 30;

def main(): 
	
	#ROOT::Math::MinimizerOptions::SetDefaultMinimizer(fitter)
	# ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2");

	# output directories
	if not os.path.exists('plots'+str(options.jetNum)+'/results'): os.makedirs('plots'+str(options.jetNum)+'/results')
	if not os.path.exists('plots'+str(options.jetNum)+'/yields'): os.makedirs('plots'+str(options.jetNum)+'/yields')
	if not os.path.exists('plots'+str(options.jetNum)+'/shapes'): os.makedirs('plots'+str(options.jetNum)+'/shapes')
	if not os.path.exists('plots'+str(options.jetNum)+'/rhalphabet'): os.makedirs('plots'+str(options.jetNum)+'/rhalphabet')
	if not os.path.exists('plots'+str(options.jetNum)+'/datacards'): os.makedirs('plots'+str(options.jetNum)+'/datacards')

	idir = "../sklimming/sklim-v0-Jun16";
	# idir = "/tmp/cmantill/"

	####################################################################################
	# do mc looping - a class that holds histograms
	bkgContainers = None;
	sigContainers = None;

	### signal corrections and other directors
	## Mass shift [GeV]     : -0.590 +/- 0.872
	## Mass resolution SF:  1.094 +/- 0.123	
	# -1.6 GeV +/- 1.2 GeV
	sig_mass_shift = 0.99;
	sig_mass_shift_unc = 0.015;
	# sig_res_shift = 0.95;
	sig_res_shift = 1.094;
	sig_res_shift_unc = 0.123;
	qcdSF = 100;
	
	if options.doMCLooping: 

		sigContainers = [];
		sigNames = [];

		sigNames.append("VectorDiJet1Jet_M50_mc.root")
		sigNames.append("VectorDiJet1Jet_M75_mc.root")
		sigNames.append("VectorDiJet1Jet_M100_mc.root")		
		sigNames.append("VectorDiJet1Jet_M125_mc.root")		
		sigNames.append("VectorDiJet1Jet_M150_mc.root")		
		sigNames.append("VectorDiJet1Jet_M200_mc.root")	
		sigNames.append("VectorDiJet1Jet_M250_mc.root")		
		sigNames.append("VectorDiJet1Jet_M300_mc.root")	
		sigLabels = ["Z\'(50 GeV)","Z\'(75 GeV)","Z\'(100 GeV)","Z\'(125 GeV)","Z\'(150 GeV)","Z\'(200 GeV)","Z\'(250 GeV)","Z\'(300 GeV)"]
		sigTags = ["Zprime50","Zprime75","Zprime100","Zprime125","Zprime150","Zprime200","Zprime250","Zprime300"];

		for i in range(0,len(sigNames)):
			sigContainers.append( MCContainer( idir+"/"+sigNames[i], float(options.lumi), sigLabels[i], sigTags[i], 1, False, options.jetNum,NMassBins[i] ) );
			# k-factor is 1.2
			sigContainers[i].morphSignal("h_peakshape",sigmass[i],
										               sig_mass_shift,sig_mass_shift_unc,
										               sig_res_shift,sig_res_shift_unc);

			# hsig = [];
			# hsig.append( getattr( sigContainers[i], "h_peakshape" ) );
			# hsig.append( getattr( sigContainers[i], "h_peakshape_shiftUp" ) );
			# hsig.append( getattr( sigContainers[i], "h_peakshape_smearUp" ) );
			# hsig.append( getattr( sigContainers[i], "h_peakshape_shiftDn" ) );
			# hsig.append( getattr( sigContainers[i], "h_peakshape_smearDn" ) );

			# makeCanvasShapeComparison(hsig,["cen","shiftup","smearup","shiftdn","smeardn"],"mcsignalshapes_"+sigTags[i],"plots"+str(options.jetNum)+"/shapes/");

			# dummyaxis = sigContainers[i].h_peakshape_matched.GetXaxis();

		###### creating interpolated sig containers
		signalMorphersHists  = [];
		signalMorphersMasses = [];
		for i in range(6):
			signalMorphersHists.append( sigContainers[i].h_peakshape_matched );
			signalMorphersMasses.append( sigmass[i] );

		morphedHistContainer = hist(signalMorphersMasses,signalMorphersHists);

		ctmp = ROOT.TCanvas("ctmp","ctmp",1000,800);
		htmp = morphedHistContainer.morph(165.);
		htmp.SetLineColor(6);
		htmp.Draw();
		morphedHistContainer.morph(180.).Draw("sames");
		ctmp.SaveAs("mtmp.pdf");


		
		interpolatedMasses = [60.0,90.0,110.0,135.0,165.0,180.0];
		additionalSigContainers = [];
		sigLabels = ["Z\'(60 GeV)","Z\'(90 GeV)","Z\'(110 GeV)","Z\'(135 GeV)","Z\'(165 GeV)","Z\'(180 GeV)"];
		sigTags = ["Zprime60","Zprime90","Zprime110","Zprime135","Zprime165","Zprime180"];
		interpolatedMasses_nbins = [60,60,60,60,60,60];
		isMorphed=True;
		for i in range(len(interpolatedMasses)):
			additionalSigContainers.append( MCContainer( 'notapplicable', float(options.lumi), sigLabels[i], sigTags[i], 1, False, options.jetNum,interpolatedMasses_nbins[i],isMorphed ) );
			

			print "interpolating, ",sigTags[i],interpolatedMasses[i],sigLabels[i]
			additionalSigContainers[i].morphSignal( "h_peakshape",interpolatedMasses[i],
										               sig_mass_shift,sig_mass_shift_unc,
										               sig_res_shift,sig_res_shift_unc,
										               morphedHistContainer.morph(interpolatedMasses[i]));

			# hsig = [];
			# hsig.append( getattr( additionalSigContainers[i], "h_peakshape" ) );
			# hsig.append( getattr( additionalSigContainers[i], "h_peakshape_shiftUp" ) );
			# hsig.append( getattr( additionalSigContainers[i], "h_peakshape_smearUp" ) );
			# hsig.append( getattr( additionalSigContainers[i], "h_peakshape_shiftDn" ) );
			# hsig.append( getattr( additionalSigContainers[i], "h_peakshape_smearDn" ) );

			# makeCanvasShapeComparison(hsig,["cen","shiftup","smearup","shiftdn","smeardn"],"mcsignalshapes_"+sigTags[i],"plots"+str(options.jetNum)+"/shapes/");

		for i in range(len(additionalSigContainers)):
			for j in range(len(sigmass)-1):
				if interpolatedMasses[i] > sigmass[j] and interpolatedMasses[i] < sigmass[j+1]:
					sigContainers.insert(j+1,additionalSigContainers[i]);
					sigmass.insert(j+1,interpolatedMasses[i]);
					NMassBins.insert(j+1,interpolatedMasses_nbins[i]);
					print "newsigmass = ", sigmass
					break;

		for s in sigContainers:
			print s._tag,s._name;


		bkgContainers = [];
		bkgNames = ["QCD.root","W.root","DY.root","TTT.root"];
		bkgLabels = ["QCD","W(qq)","Z+jets","top"];
		bkgTags = ["QCD","Winc","Zinc","top"];
		bkgmass = [0.0,80.4,91.2,80.4];
		bkgsf = [1.,0.95,0.95,0.95]; # put in the W tag SF! 
		# bkgNames = ["QCD.root","W.root","DY.root"];
		# bkgLabels = ["QCD","W(qq)","Z+jets"];
		# bkgTags = ["QCD","Winc","Zinc"];
		# bkgmass = [0.0,80.4,91.2];
		# bkgsf = [1.,0.95,0.95]; # put in the W tag SF! 
		for i in range(0,len(bkgNames)):
			tmpsf = qcdSF;
			if i > 0: tmpsf = 1;
			bkgContainers.append( MCContainer( idir+"/"+bkgNames[i], float(options.lumi)*bkgsf[i], bkgLabels[i], bkgTags[i], tmpsf, False, options.jetNum,NMassBins[sigmass.index(options.ZPrimeMass)] ) );
			# bkgContainers.append( MCContainer( idir+"/"+bkgNames[i], float(options.lumi), bkgLabels[i], bkgTags[i], tmpsf, False, options.jetNum, 60 ) );
			if i == 1 or i == 2: 
				bkgContainers[i].morphSignal("h_peakshape",bkgmass[i],
											               sig_mass_shift,sig_mass_shift_unc,
											               sig_res_shift,sig_res_shift_unc);		

	####################################################################################
	# do background estimation
	theRhalphabet = None;
	if options.doRhalphabet: 
		print "Now doing the rhalphabet!"
		if not options.qcdClosure:		
			isData = True;
			extractTFs = True;
			# for i in range(len(sigmass)):
			theRhalphabet = rhalphabet(idir+"/"+"JetHT350.root",1,"rhalphabet",1, extractTFs, options.jetNum,int(options.ZPrimeMass),isData,NMassBins[sigmass.index(options.ZPrimeMass)]) ;
			theRhalphabet.GetPredictedDistributions( idir+"/"+"JetHT.root", 1, 1, isData);
			# extractTFs = False;

			# theRhalphabet.append( rhalphabet(idir+"/"+"JetHT.root",1,"rhalphabet",1, extractTFs, options.jetNum, 85. ); #add a prediction without the W
			# theRhalphabet[len(sigmass)].GetPredictedDistributions( idir+"/"+"JetHT.root", 1, 5, isData);

		# there is a flag to do a closure test as well
		if options.qcdClosure:
			isData = False;
			extractTFs = False;
			theRhalphabet = rhalphabet(idir+"/"+"QCD.root",options.lumi,"rhalphabetClosure",1, extractTFs, options.jetNum,int(options.ZPrimeMass),isData,NMassBins[sigmass.index(options.ZPrimeMass)]) ;
			theRhalphabet.GetPredictedDistributions( idir+"/"+"QCD.root", options.lumi, 100, isData );

	####################################################################################
	# do the loop on data
	theData = None;
	if options.doData:
		print "Now doing the data!...."
		isData = True;
		theData = MCContainer( idir+"/"+"JetHT.root", 1, "data" ,"data" , 1, isData, options.jetNum, NMassBins[sigmass.index(options.ZPrimeMass)], False, True );

	####################################################################################
	# do some plotting
	if options.doCards: 
		buildDatacards(bkgContainers,sigContainers[sigmass.index(options.ZPrimeMass)],theRhalphabet,theData,options.jetNum);

	####################################################################################
	# do some plotting
	if options.doPlots: 
		BuildPlots(bkgContainers,sigContainers,theRhalphabet,theData);


def BuildPlots(bkgContainers,sigContainers,theRhalphabet,theData):

	print "making plots...";

	if options.doMCLooping and options.doRhalphabet and options.doData: 

		# sigmass = [50.,75.,100.,125.,150.,200.,250.,300.]
		# sigmass = [50.,75.,100.,125.,150.,200.]
		# sigmass = [300.]
		amBlind = False;
		# for i in range(len(sigmass)):

		makeCanvasDataMC_wpred( theData.h_jetmsd_passcut,
								theRhalphabet.grpred_jetmsd, 
								[bkgContainers[0].h_jetmsd_passcut],
								['qcd'],
								'jetmsd_pred'+str(options.jetNum)+"_"+str(int(options.ZPrimeMass)),
								'plots'+str(options.jetNum)+'/results/',
								amBlind);

		makeCanvasDataMC_wpred( theData.h_rhoDDT_passcut,
								theRhalphabet.grpred_rhoDDT, 
								[bkgContainers[0].h_rhoDDT_passcut],
								['qcd'],
								'rhoDDT_pred'+str(options.jetNum)+"_"+str(int(options.ZPrimeMass)),
								'plots'+str(options.jetNum)+'/results/',
								amBlind);

		makeCanvasDataMC_MONEY( theData.h_jetmsd_passcut,
								theRhalphabet.grpred_jetmsd, 
								# [bkgContainers[1].h_jetmsd_passcut,bkgContainers[2].h_jetmsd_passcut,bkgContainers[3].h_jetmsd_passcut,sigContainers[sigmass.index(options.ZPrimeMass)].h_peakshape],
								[bkgContainers[1].h_jetmsd_passcut,bkgContainers[2].h_jetmsd_passcut,sigContainers[sigmass.index(options.ZPrimeMass)].h_peakshape],
								# ['W(qq)','Z(qq)','top','Z\' signal'],
								['W(qq)','Z(qq)','Z\' signal'],
								'jetmsd_final'+str(options.jetNum)+"_"+str(int(options.ZPrimeMass)),
								'plots'+str(options.jetNum)+'/results/',
								amBlind);		

	if options.doMCLooping and options.doRhalphabet and not options.doData:

		for i in range(len(sigContainers)):

			makeCanvasDataMC_wpred( bkgContainers[0].h_jetmsd_passcut,
									theRhalphabet.grpred_jetmsd, 
									[bkgContainers[0].h_jetmsd_passcut],
									['qcd'],
									'jetmsd_pred'+str(options.jetNum),
									'plots'+str(options.jetNum)+'/results/',
									False);

			histholder = ROOT.TFile("qcdhists.root","RECREATE");
			histholder.cd();
			bkgContainers[0].h_jetmsd_passcut.Write();
			theRhalphabet.grpred_jetmsd.Write();
			theRhalphabet.hpred_jetmsd.Write();
			histholder.Close();
	

			# makeCanvasDataMC_wpred( bkgContainers[0].h_rhoDDT_passcut,
			# 						theRhalphabet.grpred_rhoDDT, 
			# 						[bkgContainers[0].h_rhoDDT_passcut],
			# 						['qcd'],
			# 						'rhoDDT_pred'+str(options.jetNum),
			# 						'plots'+str(options.jetNum)+'/results/',
			# 						False);

			# makeCanvasDataMC_MONEY( bkgContainers[0].h_jetmsd_passcut,
			# 						theRhalphabet[i].grpred_jetmsd, 
			# 						[bkgContainers[1].h_jetmsd_passcut,bkgContainers[2].h_jetmsd_passcut,sigContainers[0].h_jetmsd_passcut,sigContainers[2].h_jetmsd_passcut],
			# 						['W(qq)','Z(qq)','Z\' (50 GeV)','Z\' (150 GeV)'],
			# 						'jetmsd_final'+str(options.jetNum),
			# 						'plots'+str(options.jetNum)+'/results/',
			# 						False);

	# if options.doMCLooping and options.doData:

	# 	names = [];
	# 	names.append( "h_jetpt" );
	# 	names.append( "h_jeteta" );
	# 	names.append( "h_jett21" );
	# 	names.append( "h_jett21DDT" );
	# 	names.append( "h_jetmsd" );

	# 	for n in names: 
	
	# 		harray = [];
	# 		hlabels = [];
	# 		for b in bkgContainers: 
	# 			harray.append( getattr( b, n ) );
	# 			hlabels.append( b._name );
	# 		# for s in sigContainers: 
	# 		# 	harray.append( getattr( s, n ) );
	# 		# 	hlabels.append( s._name );
	# 		hd = getattr( theData, n );
	# 		makeCanvasDataMC(hd,harray,hlabels,"mc_"+n,"plots/yields/");
	# 		makeCanvasDataMC(hd,harray,hlabels,"mc_"+n,"plots/shapes/");
			
	if options.doMCLooping:

		names = [];
		names.append( "h_jetmsd" );
		names.append( "h_jetmsd_passcut" );
		names.append( "h_jett21_masswindow" );
		names.append( "h_jett21DDT_masswindow" );		
		names.append( "h_peakshape" );		

		harrays = [];
		for n in names: 
	
			harray = [];
			hlabels = [];
			# for b in bkgContainers: 
			# 	harray.append( getattr( b, n ) );
			# 	hlabels.append( b._name );
			for i,s in enumerate(sigContainers):
				if i > 3: 
					curh = getattr( s, n );
					curh.GetXaxis().SetRangeUser(70,330);
					curh.GetYaxis().SetTitle( "fraction of events" );

					harray.append( curh );
					hlabels.append( s._name );

			makeCanvasComparison(harray,hlabels,"mc_"+n,"plots"+str(options.jetNum)+"/yields/");
			makeCanvasShapeComparison(harray,hlabels,"mc_"+n,"plots"+str(options.jetNum)+"/shapes/");

			harrays.append(harray);
			
		# #roc comp
		# gr_t21    = makeROCFromHisto([harrays[2][1],harrays[2][0]], False);
		# gr_t21ddt = makeROCFromHisto([harrays[3][1],harrays[3][0]], False);
		# plotROCs([gr_t21,gr_t21ddt],["t21","t21DDT"],"plots"+str(options.jetNum)+'/shapes/','roccomp')


#----------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	main();
#----------------------------------------------------------------------------------------------------------------
