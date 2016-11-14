import ROOT
import sys
import math
import array
# from ROOT import RooRealVar,RooDataSet,RooArgSet,RooArgList,RooGenericPdf

from plotHelpers import makeCanvas, makeCanvasDataMC, makeCanvasShapeComparison,makeCanvas2D
from Fits import LinearFit,QuadraticFit

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--doMCLooping', action='store_true', dest='doMCLooping', default=False, help='go!')
parser.add_option('--doRhalphabet', action='store_true', dest='doRhalphabet', default=False, help='go!')
parser.add_option('--doData', action='store_true', dest='doData', default=False, help='go!')
parser.add_option('--doPlots', action='store_true', dest='doPlots', default=False, help='go!')
parser.add_option('--doCards', action='store_true', dest='doCards', default=False, help='go!')

parser.add_option('--extractTF', action='store_true', dest='extractTF', default=False, help='go!')
parser.add_option("--lumi", dest="lumi", default = 0.44,help="mass of LSP", metavar="MLSP")
parser.add_option("--rholo", dest="rholo", default = 0.,help="mass of LSP", metavar="MLSP")
parser.add_option("--rhohi", dest="rhohi", default = 6.,help="mass of LSP", metavar="MLSP")
parser.add_option("--DDTcut", dest="DDTcut", default = 0.38,help="mass of LSP", metavar="MLSP")
parser.add_option('--qcdClosure', action='store_true', dest='qcdClosure', default=False, help='go!')
parser.add_option("--jetNum", dest="jetNum", default = 0,help="mass of LSP", metavar="MLSP")
parser.add_option("--ZPrimeMass", dest="ZPrimeMass", default = 50., type=float,help="mass of LSP", metavar="MLSP")

(options, args) = parser.parse_args()

class rhalphabet:

	################################################################################################
	# init
	def __init__( self , filename, lumi, name, scaleFactor = 1, extractTFs = True,jetNum=0,zprimemass=85.,isData = False,nmassbins=60):

		## fixed parameters 
		# self._ptbins = 5;
		# self._ptlo = 500;
		# self._pthi = 1000;
		# self._ptbins = 9;
		# self._ptlo = 350;
		# self._pthi = 1250;
		self._ptbins = 7;
		self._ptlo = 350;
		self._pthi = 1050;
		self._nrhobins = 50;
		self._nmassbins = nmassbins;
		self._ZPrimeMass = zprimemass;
		self._isData = isData;

		self._filename = filename;
		self._tf = ROOT.TFile( self._filename );
		self._tt = self._tf.Get("otree");
		self._scaleFactor = scaleFactor;
		self._name = name;
		self._lumi = float(lumi);
		self._jetNum = str(jetNum);

		tfopt = "READ";
		if extractTFs: tfopt = "RECREATE"
		self._storageFile = ROOT.TFile("plots"+str(self._jetNum)+"/rhalphabet/storedHistos.root",tfopt);
		self.hys = None;
		self.TF_pafa = None;

		if extractTFs: 
			self.FillTFMap();
			self._storageFile.cd();
			self.hys = makeCanvas2D( self.TF_pafa,"map","plots"+str(self._jetNum)+"/rhalphabet" );
			for h in self.hys: h.Write();
			self.TF_pafa.Write();
		else:
			self.hys = [];
			for i in range(self._ptbins): self.hys.append( self._storageFile.Get("hys"+str(i)) );
			self.TF_pafa = self._storageFile.Get("TF_pafa");			

		self.do2DFit();

		# else:
		# 	# get stored information

		self._storageFile.Close();

	################################################################################################
	# fill map
	def FillTFMap( self ):

		DDTCUT = float(options.DDTcut);

		self.TF_fail = ROOT.TH2F("TF_fail",";rho;pT",self._nrhobins,float(options.rholo),float(options.rhohi),self._ptbins,self._ptlo,self._pthi)
		self.TF_pafa = ROOT.TH2F("TF_pafa",";rho;pT",self._nrhobins,float(options.rholo),float(options.rhohi),self._ptbins,self._ptlo,self._pthi)	
		self.TF_fail.Sumw2();
		self.TF_pafa.Sumw2();

		t = self._tt;

		# looping
		nent = t.GetEntries();
		for i in range(t.GetEntries()):
			
			# preamble
			t.GetEntry(i);
			if(i % (1 * nent/100) == 0):
				sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
				sys.stdout.flush();
			if i % int(self._scaleFactor) != 0: continue;

			# cutting
			jpt = getattr(t,"bst8_PUPPIjet"+self._jetNum+"_pt");
			jmsd = getattr(t,"bst8_PUPPIjet"+self._jetNum+"_msd");		
			weight = self._scaleFactor*self._lumi*getattr(self._tt,"scale1fb")*getattr(self._tt,"kfactor")*getattr(self._tt,"puWeight");
			if self._isData: weight = 1;

			if jmsd > 0:

				jt21 = getattr(t,"bst8_PUPPIjet"+self._jetNum+"_tau21");
				rhP = math.log(jmsd*jmsd/jpt);			
				jt21P = jt21 + 0.063*rhP;

				if jt21P < DDTCUT: self.TF_pafa.Fill( rhP, jpt, weight );
				else: self.TF_fail.Fill( rhP, jpt, weight );

		print "\n";

		### IF it's data, need to subract the electroweak component! 
		if self._isData: 
			print "subtracting EWK contributions..."
			self._fileEWK = ROOT.TFile("/Users/ntran/Documents/Research/CMS/WZpToQQ/dijetsGH/dijets/sklimming/sklim-v0-Jun16/EWK.root");
			self._treeEWK = self._fileEWK.Get("otree");
			self.TF_failEWK = ROOT.TH2F("TF_failEWK",";rho;pT",self._nrhobins,float(options.rholo),float(options.rhohi),self._ptbins,self._ptlo,self._pthi)
			self.TF_pafaEWK = ROOT.TH2F("TF_pafaEWK",";rho;pT",self._nrhobins,float(options.rholo),float(options.rhohi),self._ptbins,self._ptlo,self._pthi)	
			self.TF_failEWK.Sumw2();
			self.TF_pafaEWK.Sumw2();
			print "self._treeEWK.GetEntries() = ",self._treeEWK.GetEntries()
			nent = self._treeEWK.GetEntries();
			for i in range( self._treeEWK.GetEntries() ):
				self._treeEWK.GetEntry(i);
				# preamble
				if(i % (1 * nent/100) == 0):
					sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
					sys.stdout.flush();
				if i % int(self._scaleFactor) != 0: continue;

				# cutting
				jpt = getattr(self._treeEWK,"bst8_PUPPIjet"+self._jetNum+"_pt");
				jmsd = getattr(self._treeEWK,"bst8_PUPPIjet"+self._jetNum+"_msd");		
				weight = self._lumi*getattr(self._treeEWK,"scale1fb")*getattr(self._treeEWK,"kfactor")*getattr(self._treeEWK,"puWeight");

				if jmsd > 0:

					jt21 = getattr(self._treeEWK,"bst8_PUPPIjet"+self._jetNum+"_tau21");
					rhP = math.log(jmsd*jmsd/jpt);			
					jt21P = jt21 + 0.063*rhP;

					if jt21P < DDTCUT: self.TF_pafaEWK.Fill( rhP, jpt, weight );
					else: self.TF_failEWK.Fill( rhP, jpt, weight );				

			print self.TF_fail.Integral(),self.TF_pafa.Integral(),self.TF_failEWK.Integral(),self.TF_pafaEWK.Integral();
			self.TF_fail.Add(self.TF_failEWK,-1);			
			print self.TF_fail.Integral(),self.TF_pafa.Integral(),self.TF_failEWK.Integral(),self.TF_pafaEWK.Integral();
			self.TF_pafa.Add(self.TF_pafaEWK,-1);			
			print self.TF_fail.Integral(),self.TF_pafa.Integral(),self.TF_failEWK.Integral(),self.TF_pafaEWK.Integral();

		#####!!!! this uncertinaty may be wrong, binomial instead of poisson
		
		if not self._isData:
			# make the MC uncertainties like data
			for i in range(self.TF_fail.GetNbinsX()):
				for j in range(self.TF_fail.GetNbinsY()):
					self.TF_fail.SetBinError(i+1,j+1, math.sqrt( self.TF_fail.GetBinContent(i+1,j+1) ) )
					self.TF_pafa.SetBinError(i+1,j+1, math.sqrt( self.TF_pafa.GetBinContent(i+1,j+1) ) )

		self.TF_pafa.Divide( self.TF_fail );

	################################################################################################
	# fit 2D space
	def do2DFit( self ):
		THEFITRHOMAX = 3.8;
		THEFITRHOMIN = 0.2;

		# turn the histograms into graphs, that way exclude the various mass windows
		grs = [];
		for i in range(len(self.hys)): 
			curpt = self.TF_pafa.GetYaxis().GetBinCenter(i+1);
			# xlo = math.log(70.*70./curpt);
			xlo = math.log(999.*999./curpt);
			xhi = math.log(100.*100./curpt);
			grs.append( self.convertHistToGraph( self.hys[i], xlo, xhi ) );
			# grs.append( self.convertHistToGraph( self.hys[i] ) );

		# do the fits

		##########################################################################################
		######## THE RHO FITS
		print "----------- the rho fits..."
		##########################################################################################
		self.theRhoFits = [];
		ctr = 0
		for g in grs:
			# self.theRhoFits.append( QuadraticFit([0.01,6e-5,0.01],float(options.rholo),float(options.rhohi),"linFit1", "EMRNS") );
			# g.Fit(self.theRhoFits[ctr].fit,"RQ")
			# fitter = ROOT.TVirtualFitter.GetFitter()
			# self.theRhoFits[ctr].Converter(fitter)

			# self.theRhoFits.append( ROOT.TF1("rhofit"+str(ctr),'([0]+[1]*x+[2]*x*x)',float(options.rholo),float(options.rhohi)) ); 	
			cur_rhofitmin = THEFITRHOMIN+float(ctr)*0.1;
			# cur_rhofitmin = THEFITRHOMIN;
			cur_rhofitmax = THEFITRHOMAX+float(ctr)*0.2;
			print "rho limits = ",cur_rhofitmin,cur_rhofitmax
			self.theRhoFits.append( ROOT.TF1("rhofit"+str(ctr),'([0]+[1]*x+[2]*x*x)',cur_rhofitmin,cur_rhofitmax ) ); 	
			
			self.theRhoFits[ctr].SetParameter(0,0.02);
			self.theRhoFits[ctr].SetParameter(1,0.01);					
			self.theRhoFits[ctr].SetParameter(2,0.);						

			# self.theRhoFits[ctr].SetParLimits(0,0.0,0.1);
			# self.theRhoFits[ctr].SetParLimits(1,-0.1,0.1);					
			# self.theRhoFits[ctr].SetParLimits(2,-0.1,0.1);					

			g.Fit(self.theRhoFits[ctr],"R");
			par0 = self.theRhoFits[ctr].GetParameter(0);
			par1 = self.theRhoFits[ctr].GetParameter(1);
			par2 = self.theRhoFits[ctr].GetParameter(2);

			print "Fit statistics = "; 
			print "      NDOF = ", self.theRhoFits[ctr].GetNDF();
			print "      ChiS = ", self.theRhoFits[ctr].GetChisquare();
			print "      Prob = ", self.theRhoFits[ctr].GetProb();

			ctr+=1;		

		for i in range(len(grs)):
			curcan = ROOT.TCanvas("curcan"+str(i),"curcan"+str(i),1200,800);
			curhr2 = curcan.DrawFrame(float(options.rholo),0,float(options.rhohi),0.1);
			curhr2.GetYaxis().SetTitle("npass/nfail in a pT bin");
			curhr2.GetXaxis().SetTitle("#rho^{DDT}");	
			curcan.SetGrid();
			grs[i].Draw("pe");
			
			# self.theRhoFits[i].fit.SetLineColor(2);	
			# self.theRhoFits[i].ErrUp.SetLineColor(2);	
			# self.theRhoFits[i].ErrDn.SetLineColor(2);	
			# self.theRhoFits[i].ErrUp.SetLineStyle(2);	
			# self.theRhoFits[i].ErrDn.SetLineStyle(2);	
			# self.theRhoFits[i].fit.Draw("sames");	
			# self.theRhoFits[i].ErrUp.Draw("sames");	
			# self.theRhoFits[i].ErrDn.Draw("sames");	

			self.theRhoFits[i].Draw("sames");

			curcan.SaveAs("plots"+str(self._jetNum)+"/rhalphabet/map_rhodependence_bin"+str(i)+"_"+self._jetNum+"_"+str(self._ZPrimeMass)+".pdf")

		##############################
		hpt_par0 = ROOT.TH1F("hpt_par0",";pT;par0",self._ptbins,self._ptlo,self._pthi);
		hpt_par1 = ROOT.TH1F("hpt_par1",";pT;par1",self._ptbins,self._ptlo,self._pthi);
		hpt_par2 = ROOT.TH1F("hpt_par2",";pT;par2",self._ptbins,self._ptlo,self._pthi);
		
		binno = 1;
		for f in self.theRhoFits:

			hpt_par0.SetBinContent( binno, f.GetParameter(0) );
			hpt_par0.SetBinError( binno, f.GetParError(0) );
			hpt_par1.SetBinContent( binno, f.GetParameter(1) );
			hpt_par1.SetBinError( binno, f.GetParError(1) );
			hpt_par2.SetBinContent( binno, f.GetParameter(2) );
			hpt_par2.SetBinError( binno, f.GetParError(2) );

			binno += 1;

		hpt_pars = [hpt_par0,hpt_par1,hpt_par2];

		##########################################################################################
		######## THE pT FITS
		print "--------the Pt fits..."
		##########################################################################################		
		self.thePtFits = [];
		ctr = 0;
		for h in hpt_pars:
			# if ctr == 0: self.thePtFits.append( QuadraticFit([0.01,6e-5],self._ptlo,self._pthi,"linFit2", "EMRNS") );
			# else: self.thePtFits.append( QuadraticFit([0.01,6e-5,0.01],self._ptlo,self._pthi,"linFit2", "EMRNS") );
			
			# # self.thePtFits[ctr].SetParLimits(0,0.005,0.05);
			# curfitresult = h.Fit(self.thePtFits[ctr].fit,"RS")
			# fitter = ROOT.TVirtualFitter.GetFitter()
			# self.thePtFits[ctr].Converter(fitter)
			self.thePtFits.append( ROOT.TF1("ptfit"+str(ctr),'([0]+[1]*x+[2]*x*x+[3]*x*x*x)',self._ptlo,self._pthi ) ); 	

			self.thePtFits[ctr].SetParameter(0,0.05);
			self.thePtFits[ctr].SetParameter(1,0.005);					
			self.thePtFits[ctr].SetParameter(2,0.);						
			self.thePtFits[ctr].SetParameter(3,0.);						
			# self.thePtFits[ctr].SetParLimits(0,0.006,0.1);
			# self.thePtFits[ctr].SetParLimits(1,-0.1,0.1);					
			# self.thePtFits[ctr].SetParLimits(2,-0.1,0.1);					

			h.Fit(self.thePtFits[ctr],"R");
			# par0 = self.thePtFits[ctr].GetParameter(0);
			# par1 = self.thePtFits[ctr].GetParameter(1);
			# par2 = self.thePtFits[ctr].GetParameter(2);

			# print "Fit statistics = "; 
			# print "      NDOF = ", self.thePtFits[ctr].fit.GetNDF();
			# print "      ChiS = ", self.thePtFits[ctr].fit.GetChisquare();
			# print "      Prob = ", self.thePtFits[ctr].fit.GetProb();

			ctr+=1

		for i in range(len(hpt_pars)):
			curcan2 = ROOT.TCanvas("curcan2"+str(i),"curcan2"+str(i),1200,800);
			# curhr2 = curcan.DrawFrame(ptlo[0],0,pthi[nptbins-1],0.2);
			# curhr2.GetYaxis().SetTitle("npass/nfail in a #rho^{DDT} bin");
			# curhr2.GetXaxis().SetTitle("pT");	
			curcan2.SetGrid();
			hpt_pars[i].Draw("pe");
			# self.thePtFits[i].fit.SetLineColor(2);	
			# self.thePtFits[i].ErrUp.SetLineColor(2);	
			# self.thePtFits[i].ErrDn.SetLineColor(2);	
			# self.thePtFits[i].ErrUp.SetLineStyle(2);	
			# self.thePtFits[i].ErrDn.SetLineStyle(2);	
			# self.thePtFits[i].fit.Draw("sames");	
			# self.thePtFits[i].ErrUp.Draw("sames");	
			# self.thePtFits[i].ErrDn.Draw("sames");	
			self.thePtFits[i].Draw("sames");
			curcan2.SaveAs("plots"+str(self._jetNum)+"/rhalphabet/map_ptdependence_par"+str(i)+"_"+self._jetNum+"_"+str(self._ZPrimeMass)+".pdf")	

		# self.effPlane = ROOT.TF2("TransferPlane", "(([0]+ [1]*y + [2]*y*y + [9]*y*y*y) + ([3]+ [4]*y + [5]*y*y + [10]*y*y*y)*x + ([6]+ [7]*y + [8]*y*y + [11]*y*y*y)*x*x + ([12]+ [13]*y + [14]*y*y)*x*x*x)",THEFITRHOMIN,6.,self._ptlo,self._pthi);
		self.effPlane = ROOT.TF2("TransferPlane", "(([0]+ [1]*y + [2]*y*y + [9]*y*y*y) + ([3]+ [4]*y + [5]*y*y + [10]*y*y*y)*x + ([6]+ [7]*y + [8]*y*y + [11]*y*y*y)*x*x)",THEFITRHOMIN,6.,self._ptlo,self._pthi);

		# par0_2d = self.thePtFits[0].GetParameter(0);
		# par1_2d = self.thePtFits[0].GetParameter(1);
		# par2_2d = self.thePtFits[0].GetParameter(2);
		
		# par3_2d = self.thePtFits[1].GetParameter(0);
		# par4_2d = self.thePtFits[1].GetParameter(1);
		# par5_2d = self.thePtFits[1].GetParameter(2);
		
		# par6_2d = self.thePtFits[2].GetParameter(0);
		# par7_2d = self.thePtFits[2].GetParameter(1);
		# par8_2d = self.thePtFits[2].GetParameter(2);

		# par9_2d  = self.thePtFits[0].GetParameter(3);
		# par10_2d = self.thePtFits[1].GetParameter(3);
		# par11_2d = self.thePtFits[2].GetParameter(3);

		# self.effPlane.SetParameter(0,par0_2d);
		# # self.effPlane.SetParError(0,self.thePtFits[0].GetParError(0));
		# self.effPlane.SetParameter(1,self.thePtFits[0].GetParameter(1));
		# # self.effPlane.SetParError(1,self.thePtFits[0].GetParError(1));
		# self.effPlane.SetParameter(2,self.thePtFits[0].GetParameter(2));
		# # self.effPlane.SetParError(2,self.thePtFits[0].GetParError(2)); 

		# self.effPlane.SetParameter(3,self.thePtFits[1].GetParameter(0));
		# # self.effPlane.SetParError(3,self.thePtFits[1].GetParError(0));
		# self.effPlane.SetParameter(4,self.thePtFits[1].GetParameter(1));
		# # self.effPlane.SetParError(4,self.thePtFits[1].GetParError(1));
		# self.effPlane.SetParameter(5,self.thePtFits[1].GetParameter(2));
		# # self.effPlane.SetParError(5,self.thePtFits[1].GetParError(2));

		# self.effPlane.SetParameter(6,self.thePtFits[2].GetParameter(0));
		# # self.effPlane.SetParError(6,self.thePtFits[2].GetParError(0));
		# self.effPlane.SetParameter(7,self.thePtFits[2].GetParameter(1));
		# # self.effPlane.SetParError(7,self.thePtFits[2].GetParError(1));
		# self.effPlane.SetParameter(8,self.thePtFits[2].GetParameter(2));
		# # self.effPlane.SetParError(8,self.thePtFits[2].GetParError(2));

		# self.effPlane.SetParameter(9,self.thePtFits[0].GetParameter(3));
		# self.effPlane.SetParameter(10,self.thePtFits[1].GetParameter(3));
		# self.effPlane.SetParameter(11,self.thePtFits[2].GetParameter(3));

		# self.effPlane.SetParameter(12,0.0);
		# self.effPlane.SetParameter(13,0.0);
		# self.effPlane.SetParameter(14,0.0);
		# self.effPlane.SetParError(9,self.thePtFits[0].GetParError(3));
		# self.effPlane.SetParError(10,self.thePtFits[1].GetParError(3));
		# self.effPlane.SetParError(11,self.thePtFits[2].GetParError(3));

		# self.effPlane.SetParameter(11,0);
		# self.effPlane.FixParameter(6,0);
		# self.effPlane.FixParameter(7,0);
		# self.effPlane.FixParameter(8,0);
		# self.effPlane.FixParameter(9,0);
		# self.effPlane.FixParameter(10,0);
		# self.effPlane.FixParameter(11,0);
		# self.effPlane.FixParameter(12,0);
		# self.effPlane.FixParameter(13,0);
		# self.effPlane.FixParameter(14,0);

		self.effPlane_prefit = self.effPlane.Clone();

		txta = ROOT.TLatex(0.05,0.90,"CMS");
		txta.SetNDC();
		txtb = ROOT.TLatex(0.12,0.90,"Preliminary");
		txtb.SetNDC(); txtb.SetTextFont(52);
		txta.SetTextSize(0.040);
		txtb.SetTextSize(0.040);
		cplane = ROOT.TCanvas("cplane","cplane",1000,800);
		cplane.SetBottomMargin(0.15);
		cplane.SetTopMargin(0.08);
		self.effPlane.SetTitle(";#rho^{DDT};pT;N_{pass}/N_{fail}")
		self.effPlane.GetXaxis().SetTitleOffset(1.0)
		self.effPlane.GetYaxis().SetTitleOffset(1.0)
		self.effPlane.GetZaxis().SetTitleOffset(1.0)
		self.effPlane.GetXaxis().SetLabelSize(0.03);
		self.effPlane.GetYaxis().SetLabelSize(0.03);
		self.effPlane.GetZaxis().SetLabelSize(0.03);
		self.effPlane.Draw("surf1");
		txta.Draw();
		txtb.Draw();
		# ROOT.gPad.SetPhi(0.6);
		cplane.SaveAs("plots"+str(self._jetNum)+"/rhalphabet/map_parameterized"+self._jetNum+"_"+str(self._ZPrimeMass)+".pdf")			

		# need to update this to the TGraph2D with the mass windows missing?
		self.TF_pafa_gr2D = self.convertHistToGraph2D(self.TF_pafa,91,89,self._ZPrimeMass);
		self.theFitResult0 = self.TF_pafa_gr2D.Fit(self.effPlane,"RS");
		self.theFitResult = self.TF_pafa_gr2D.Fit(self.effPlane,"RS");

		# get baseline parameters:
		self.base_parameters = [];
		for i in range(self.effPlane.GetNpar()):
			self.base_parameters.append(self.effPlane.GetParameter(i));

		lCovMatrix      = self.theFitResult.GetCovarianceMatrix();
		lCovMatrixEigen = ROOT.TMatrixDSymEigen(lCovMatrix);
		lEigVecs = lCovMatrixEigen.GetEigenVectors();
		lEigVals = lCovMatrixEigen.GetEigenValues();
		for v in lEigVals: print "eigen values = ", v
		print lEigVecs.GetNrows(),lEigVecs.GetNcols();

		
		self.eigen_parameters_up = [];
		self.eigen_parameters_dn = [];
		for i in range(lEigVecs.GetNcols()):
			cureigen_up = []
			cureigen_dn = []
			for j in range(lEigVecs.GetNrows()):
				
				# if self.effPlane.GetParameter(j) != 0: print round(lEigVecs(j,i),10),
				# else: print "0.0",
				cureigen = lEigVals(j);
				if cureigen < 0: cureigen = 0.;
				cureigen_up.append( self.effPlane.GetParameter(j) + math.sqrt(cureigen)*lEigVecs(j,i)*self.effPlane.GetParameter(j) );
				cureigen_dn.append( self.effPlane.GetParameter(j) - math.sqrt(cureigen)*lEigVecs(j,i)*self.effPlane.GetParameter(j) );
				# cureigen_up.append( self.effPlane.GetParameter(j) + math.sqrt(lEigVals(j))*lEigVecs(j,i) );
				# cureigen_dn.append( self.effPlane.GetParameter(j) - math.sqrt(lEigVals(j))*lEigVecs(j,i) );

			# print "\n"
			self.eigen_parameters_up.append(cureigen_up);
			self.eigen_parameters_dn.append(cureigen_dn);
		
		# print self.base_parameters
		# print len(self.eigen_parameters_up),len(self.eigen_parameters_dn)
		# print self.eigen_parameters_up;

		# TMatrixD  lEigVecs(3,3);    
		# lEigVecs = TMatrixDSymEigen(lCovMatrix).GetEigenVectors();
		# TVectorD  lEigVals(3);      
		# lEigVals = TMatrixDSymEigen(lCovMatrix).GetEigenValues();


		print "Fit statistics = "; 
		print "      NDOF = ", self.effPlane.GetNDF();
		print "      ChiS = ", self.effPlane.GetChisquare();
		print "      Prob = ", self.effPlane.GetProb();
		print "      GOF  = ", 2.*self.theFitResult.MinFcnValue();
		# curpoint = array.array('d', []);
		# curpoint.append( 2 );
		# curpoint.append( 1400 );
		# curerr = array.array('d', []);
		# curerr.append(0.)
		# r.GetConfidenceIntervals(1,1,1,curpoint,curerr, 0.683, False);
		# print curerr

		cplane2 = ROOT.TCanvas("cplane2","cplane2",1000,800);
		self.effPlane.SetTitle(";#rho^{DDT};p_{T} (GeV);N_{pass}/N_{fail}")
		cplane2.SetBottomMargin(0.15);
		cplane2.SetTopMargin(0.08);		
		# self.TF_pafa_gr2D.Draw("pe");
		self.effPlane.GetXaxis().SetTitleOffset(1.2)
		self.effPlane.GetYaxis().SetTitleOffset(1.5)
		self.effPlane.GetZaxis().SetTitleOffset(1.0)
		self.effPlane.GetXaxis().SetLabelSize(0.03);
		self.effPlane.GetYaxis().SetLabelSize(0.03);
		self.effPlane.GetZaxis().SetLabelSize(0.03);		
		self.effPlane.Draw("surf1");
		txta.Draw();
		txtb.Draw();
		# ROOT.gPad.SetPhi(0.6);
		cplane2.SaveAs("plots"+str(self._jetNum)+"/rhalphabet/map_parameterized_postfit"+self._jetNum+"_"+str(self._ZPrimeMass)+".pdf")	

		#make a diff of the parameterized plane and the original
		self.TF_diff = ROOT.TH2F("TF_diff",";rho;pT",20,float(options.rholo),float(options.rhohi),self._ptbins,self._ptlo,self._pthi)	
		for i in range(self.TF_pafa.GetNbinsX()):
			for j in range(self.TF_pafa.GetNbinsY()):
				paramval = self.effPlane.Eval( self.TF_pafa.GetXaxis().GetBinCenter(i+1) , self.TF_pafa.GetYaxis().GetBinCenter(j+1) )
				tableval = self.TF_pafa.GetBinContent(i+1,j+1);
				if tableval != 0: self.TF_diff.SetBinContent( i+1,j+1, (paramval-tableval)/tableval );
				else: self.TF_diff.SetBinContent( i+1,j+1, 0 );
		
		cplane3 = ROOT.TCanvas("cplane3","cplane3",1000,800);
		# ROOT.gPad.GetRightMargin(0.18);
		self.TF_diff.SetTitle(";#rho^{DDT};pT;N_{pass}/N_{fail}")
		self.TF_diff.Draw("colz");
		cplane3.SaveAs("plots"+str(self._jetNum)+"/rhalphabet/map_diff"+self._jetNum+"_"+str(self._ZPrimeMass)+".pdf")
		print self.TF_pafa	

		# #make a diff of the parameterized plane and the original
		# chivals = [];
		# chibyhand = 0.;
		# ndofbyhand = 0;
		# for i in range(self.TF_pafa_gr2D.GetN()):
		# 	if self.TF_pafa_gr2D.GetX()[i] < 0.2 or self.TF_pafa_gr2D.GetX()[i] > 3.8: continue;
		# 	curval,curerr = self.predictedVal(self.TF_pafa_gr2D.GetX()[i],self.TF_pafa_gr2D.GetY()[i]);
		# 	# tmperr = math.sqrt((0.15*curval)*(0.15*curval)+(curerr/curval)*(curerr/curval));
		# 	# curerr = curval*tmperr;

		# 	toterror1 = curerr*curerr + self.TF_pafa_gr2D.GetEZ()[i]*self.TF_pafa_gr2D.GetEZ()[i];
		# 	toterror2 = self.TF_pafa_gr2D.GetEZ()[i]*self.TF_pafa_gr2D.GetEZ()[i];
		# 	# print self.TF_pafa_gr2D.GetX()[i],self.TF_pafa_gr2D.GetY()[i],round(self.TF_pafa_gr2D.GetZ()[i],6),round(self.TF_pafa_gr2D.GetEZ()[i],6), round(curval,6), round(curerr,6), "-----", round( (self.TF_pafa_gr2D.GetZ()[i]-curval)*(self.TF_pafa_gr2D.GetZ()[i]-curval)/toterror2 ,6) ;
		# 	chivals.append( (self.TF_pafa_gr2D.GetZ()[i]-curval)/math.sqrt(toterror1) );
		# 	# chibyhand += ((self.TF_pafa_gr2D.GetZ()[i]-curval)/math.sqrt(toterror2))*((self.TF_pafa_gr2D.GetZ()[i]-curval)/math.sqrt(toterror2));
		# 	chibyhand += (self.TF_pafa_gr2D.GetZ()[i]-curval)*(self.TF_pafa_gr2D.GetZ()[i]-curval)/curval;
		# 	ndofbyhand += 1;

		# h_pulls = ROOT.TH1F("fitpulls",";point no. ;pull",ndofbyhand,0.5,ndofbyhand+0.5);
		# for i in range(ndofbyhand): h_pulls.SetBinContent(i+1,chivals[i]);
		# print "chi2 by hand = ", chibyhand, ndofbyhand;
		# cplane4 = ROOT.TCanvas("cplane4","cplane4",1000,800);
		# h_pulls.Draw("hist");
		# cplane4.SaveAs("plots"+str(self._jetNum)+"/rhalphabet/map_pulls"+self._jetNum+"_"+str(self._ZPrimeMass)+".pdf")

	################################################################################################
	# predicted distribution from the given input file
	def GetPredictedDistributions(self, filename, lumi=1, scaleFactor=1, isData=False):
		
		DDTCUT = float(options.DDTcut);

		self.pred_fn = filename;
		self.pred_tf = ROOT.TFile( self.pred_fn );
		self.pred_tt = self.pred_tf.Get("otree");
		self.pred_scaleFactor = scaleFactor;

		jetmassbins = [30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,
					   95,100,105,110,115,120,125,130,135,140,145,150,
					   157,164,171,178,185,192,199,206,213,220,
					   230,240,250,260,270,280,290,300,315,330]

		nbinsmass = self._nmassbins; 
		print "[GetPredictedMassDistributions] = ",nbinsmass

		# Define histogram
		self.hpred_jetmsd         = ROOT.TH1F("hpred_jetmsd"+self._name,"; soft drop mass (GeV);", nbinsmass, 30, 330);
		self.hpred_jetmsd_errup   = ROOT.TH1F("hpred_jetmsd_errup"+self._name,"; soft drop mass (GeV);", nbinsmass, 30, 330);
		self.hpred_jetmsd_failcut = ROOT.TH1F("hpred_jetmsd_failcut"+self._name,"; soft drop mass (GeV);", nbinsmass, 30, 330);
		self.hpred_rhoDDT         = ROOT.TH1F("hpred_rhoDDT"+self._name,"; #rho^{DDT};", self._nrhobins, float(options.rholo),float(options.rhohi));
		self.hpred_rhoDDT_errup    = ROOT.TH1F("hpred_rhoDDT_errup"+self._name,"; #rho^{DDT};", self._nrhobins, float(options.rholo),float(options.rhohi));
		self.hpred_rhoDDT_failcut = ROOT.TH1F("hpred_rhoDDT_failcut"+self._name,"; #rho^{DDT};", self._nrhobins, float(options.rholo),float(options.rhohi));
		
		self.hpred_jetmsd_altup = [];
		self.hpred_jetmsd_altdn = [];
		for i in range(self.effPlane.GetNpar()): self.hpred_jetmsd_altup.append( ROOT.TH1F("hpred_jetmsd_up"+self._name+str(i),"; soft drop mass (GeV);", nbinsmass, 30, 330) );
		for i in range(self.effPlane.GetNpar()): self.hpred_jetmsd_altdn.append( ROOT.TH1F("hpred_jetmsd_dn"+self._name+str(i),"; soft drop mass (GeV);", nbinsmass, 30, 330) );

		self.hpred_jetmsd.Sumw2();
		self.hpred_jetmsd_errup.Sumw2();
		self.hpred_jetmsd_failcut.Sumw2();
		self.hpred_rhoDDT.Sumw2();
		self.hpred_rhoDDT_errup.Sumw2();
		self.hpred_rhoDDT_failcut.Sumw2();

		# looping
		nent = self.pred_tt.GetEntries();
		print "size = ",nent
		for i in range(self.pred_tt.GetEntries()):

			# if i > 20000: break;

			# preamble
			self.pred_tt.GetEntry(i);
			if(i % (1 * nent/100) == 0):
				sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
				sys.stdout.flush();
			
			if i % self.pred_scaleFactor != 0: continue;		

			jpt = getattr(self.pred_tt,"bst8_PUPPIjet"+self._jetNum+"_pt");
			jeta = getattr(self.pred_tt,"bst8_PUPPIjet"+self._jetNum+"_eta");
			jmsd = getattr(self.pred_tt,"bst8_PUPPIjet"+self._jetNum+"_msd");		
			if jmsd == 0.: jmsd = 0.01;
			# weight = float(self.pred_scaleFactor)*float(lumi)*getattr(self.pred_tt,"scale1fb");
			weight = self.pred_scaleFactor*self._lumi*getattr(self.pred_tt,"scale1fb")*getattr(self.pred_tt,"kfactor")*getattr(self.pred_tt,"puWeight");
			if isData: weight = 1;

			if jpt < 500: continue;

			jt21 = getattr(self.pred_tt,"bst8_PUPPIjet"+self._jetNum+"_tau21");
			rhP = math.log(jmsd*jmsd/jpt);                  
			jt21P = jt21 + 0.063*rhP;

			if rhP > 0.2 and jt21P > DDTCUT:
				self.hpred_jetmsd_failcut.Fill( jmsd, weight );
				self.hpred_rhoDDT_failcut.Fill( rhP, weight );
 				curval,curerr,altvals_up,altvals_dn = self.predictedVal(rhP,jpt);
 				# curerr = 0.; # this will ensure keeping only the statistical errors from fit region
				self.hpred_jetmsd.Fill( jmsd, weight*curval );
				self.hpred_rhoDDT.Fill( rhP, weight*curval );
				self.hpred_jetmsd_errup.Fill( jmsd, weight*(curval+curerr) );
				self.hpred_rhoDDT_errup.Fill( rhP, weight*(curval+curerr) );

				for i,v in enumerate(altvals_up): self.hpred_jetmsd_altup[i].Fill( jmsd, weight*v );
				for i,v in enumerate(altvals_dn): self.hpred_jetmsd_altdn[i].Fill( jmsd, weight*v );

		print "\n";
		print "predicted integral = ", self.hpred_jetmsd.Integral(), self.hpred_jetmsd.GetNbinsX();
		hpred_jetmsdClone = self.hpred_jetmsd.Clone("hpred_jetmsdClone");

		if self._isData: 
			hpred_jetmsd_WZToSubtract = ROOT.TH1F("hpred_jetmsd_WZToSubtract"+self._name,"; soft drop mass (GeV);", nbinsmass, 30, 330);
			hpred_rhoDDT_WZToSubtract = ROOT.TH1F("hpred_rhoDDT_WZToSubtract"+self._name,"; #rho^{DDT};", self._nrhobins, float(options.rholo),float(options.rhohi)); 
			print "subtracting EWK contributions from fail template..."
			self._fileEWK = ROOT.TFile("/Users/ntran/Documents/Research/CMS/WZpToQQ/dijetsGH/dijets/sklimming/sklim-v0-Jun16/EWK.root");
			self._treeEWK = self._fileEWK.Get("otree");
			nent = self._treeEWK.GetEntries();
			for i in range( self._treeEWK.GetEntries() ):		

				# preamble
				self._treeEWK.GetEntry(i);
				if(i % (1 * nent/100) == 0):
					sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
					sys.stdout.flush();
				
				jpt = getattr(self._treeEWK,"bst8_PUPPIjet"+self._jetNum+"_pt");
				jeta = getattr(self._treeEWK,"bst8_PUPPIjet"+self._jetNum+"_eta");
				jmsd = getattr(self._treeEWK,"bst8_PUPPIjet"+self._jetNum+"_msd");		
				if jmsd == 0.: jmsd = 0.01;
				weight = self._lumi*getattr(self._treeEWK,"scale1fb")*getattr(self._treeEWK,"kfactor")*getattr(self._treeEWK,"puWeight");
				# weight = weight*-1.;
				if jpt < 500: continue;

				jt21 = getattr(self._treeEWK,"bst8_PUPPIjet"+self._jetNum+"_tau21");
				rhP = math.log(jmsd*jmsd/jpt);                  
				jt21P = jt21 + 0.063*rhP;

				if rhP > 0.2 and jt21P > DDTCUT:
	 				curval,curerr,altvals_up,altvals_dn = self.predictedVal(rhP,jpt);
	 				# curerr = 0.; # this will ensure keeping only the statistical errors from fit region
					hpred_jetmsd_WZToSubtract.Fill( jmsd, weight*curval );
					hpred_rhoDDT_WZToSubtract.Fill( rhP, weight*curval );		

			for i in range(self.hpred_jetmsd.GetNbinsX()): self.hpred_jetmsd.SetBinContent(i+1, self.hpred_jetmsd.GetBinContent(i+1)-hpred_jetmsd_WZToSubtract.GetBinContent(i+1)  );
			for i in range(self.hpred_rhoDDT.GetNbinsX()): self.hpred_rhoDDT.SetBinContent(i+1, self.hpred_rhoDDT.GetBinContent(i+1)-hpred_rhoDDT_WZToSubtract.GetBinContent(i+1)  );
			for h in self.hpred_jetmsd_altup:
				for i in range(h.GetNbinsX()): h.SetBinContent(i+1, h.GetBinContent(i+1)-hpred_jetmsd_WZToSubtract.GetBinContent(i+1)  );	
			for h in self.hpred_jetmsd_altdn: 
				for i in range(h.GetNbinsX()): h.SetBinContent(i+1, h.GetBinContent(i+1)-hpred_jetmsd_WZToSubtract.GetBinContent(i+1)  );	
				
		print "\n";
		print "predicted integral after EWK subtraction = ", self.hpred_jetmsd.Integral(), self.hpred_jetmsd.GetNbinsX();

		# adjust histogram errors
		self.grpred_jetmsd = self.AdjustHistogramErrors(self.hpred_jetmsd,self.hpred_jetmsd_errup);
		self.grpred_rhoDDT = self.AdjustHistogramErrors(self.hpred_rhoDDT,self.hpred_rhoDDT_errup);

		ceigtest = ROOT.TCanvas("ceigtest","ceigtest",1000,800);
		hpred_jetmsd_altup_cloned = [];
		hpred_jetmsd_altdn_cloned = [];
		for h in self.hpred_jetmsd_altup: hpred_jetmsd_altup_cloned.append(h.Clone());
		for h in self.hpred_jetmsd_altdn: hpred_jetmsd_altdn_cloned.append(h.Clone());

		for i,h in enumerate(hpred_jetmsd_altup_cloned): 
			# h = h0.Clone();
			h.Scale(self.hpred_jetmsd.Integral()/h.Integral())
			h.Divide(self.hpred_jetmsd)
			if i == 0: 
				h.SetMaximum(1.2);
				h.SetMinimum(0.8);
				h.GetYaxis().SetTitle("Shape variations for fit eigenvectors");
			h.Draw('histsames');
		for i,h in enumerate(hpred_jetmsd_altdn_cloned): 
			# h = h0.Clone();
			h.SetLineColor(2);
			h.Scale(self.hpred_jetmsd.Integral()/h.Integral())
			h.Divide(self.hpred_jetmsd)
			h.Draw('histsames');
		# ROOT.gPad.SetLogy();			
		ceigtest.SaveAs("plots"+str(self._jetNum)+"/rhalphabet/eigenvectors"+self._jetNum+"_"+str(self._ZPrimeMass)+".pdf")

		ceigtest2 = ROOT.TCanvas("ceigtest2","ceigtest2",1000,800);
		hpred_jetmsd_altup_cloned = [];
		hpred_jetmsd_altdn_cloned = [];
		for h in self.hpred_jetmsd_altup: hpred_jetmsd_altup_cloned.append(h.Clone());
		for h in self.hpred_jetmsd_altdn: hpred_jetmsd_altdn_cloned.append(h.Clone());		
		for i,h in enumerate(hpred_jetmsd_altup_cloned): 
			h.Divide(self.hpred_jetmsd)
			if i == 0: 
				h.SetMaximum(1.2);
				h.SetMinimum(0.8);
				h.GetYaxis().SetTitle("Shape variations for fit eigenvectors");
			h.Draw('histsames');
		for i,h in enumerate(hpred_jetmsd_altdn_cloned): 
			h.SetLineColor(2);
			h.Divide(self.hpred_jetmsd)
			h.Draw('histsames');
		ceigtest2.SaveAs("plots"+str(self._jetNum)+"/rhalphabet/eigenvectors_raw_"+self._jetNum+"_"+str(self._ZPrimeMass)+".pdf")

		efname = "plots"+str(self._jetNum)+"/rhalphabet/eigenshapes"+self._jetNum+"_"+str(self._ZPrimeMass)+".root";
		efile = ROOT.TFile(efname,"RECREATE")
		for h in self.hpred_jetmsd_altup: h.Write();
		for h in self.hpred_jetmsd_altdn: h.Write();
		self.hpred_jetmsd.Write();
		hpred_jetmsdClone.Write();
		self.grpred_jetmsd.Write();
		efile.Close();



	def predictedVal(self,xval,yval):

		curpoint = array.array('d', []);
		curpoint.append( xval );
		curpoint.append( yval );
		curerr = array.array('d', []);
		curerr.append(0.)
		self.theFitResult.GetConfidenceIntervals(1,1,1,curpoint,curerr,0.683,False);	
		
		# print self.base_parameters
		# print self.eigen_parameters_up[0]
		# print self.eigen_parameters_up[10]

		altval_up = []
		altval_dn = []
		curval = self.effPlane.Eval( xval, yval );

		for i in range(len(self.base_parameters)): self.effPlane.SetParameter(i,self.base_parameters[i]);
		# print self.effPlane.Eval( xval, yval ), xval, yval

		for i in range(len(self.eigen_parameters_up)):
		
			for j in range(len(self.eigen_parameters_up[i])): self.effPlane.SetParameter(j,self.eigen_parameters_up[i][j]); 
			altval_up.append( self.effPlane.Eval( xval, yval ) );
			for j in range(len(self.eigen_parameters_dn[i])): self.effPlane.SetParameter(j,self.eigen_parameters_dn[i][j]); 
			altval_dn.append( self.effPlane.Eval( xval, yval ) );
		
		# print curval,altval_up[0],altval_up[1],altval_up[10]
		# print curval,altval_dn[0],altval_dn[1],altval_dn[10]
		# print altval_up
		# print altval_up
		# curval2 = self.effPlane_prefit.Eval( xval, yval );
		# print curval,curval2,curerr[0]

		# storageFile2 = ROOT.TFile("plots/rhalphabet/storedHistos.root","READ");
		# tfpafa = storageFile2.Get("TF_pafa");
		# # print self.TF_pafa
		# # print xval,yval,self.TF_pafa.GetTitle();#,self.TF_pafa.GetXaxis().FindBin(xval), self.TF_pafa.GetYaxis().FindBin(yval)
		# curval = tfpafa.GetBinContent( tfpafa.GetXaxis().FindBin(xval), tfpafa.GetYaxis().FindBin(yval) );

		return curval,curerr[0],altval_up,altval_dn;

	################################################################################################
	# helpers
	def convertHistToGraph( self, h, xlo=9999, xhi=-9999, zprimemass=85 ):

		x = array.array('d', [])
		y = array.array('d', [])
		xe = array.array('d', [])
		ye = array.array('d', [])
		ctr = 0;
		for i in range(h.GetNbinsX()):
			curx = h.GetBinCenter(i+1);
			if curx > xlo and curx < xhi: continue;
			
			x.append( h.GetBinCenter(i+1) );
			y.append( h.GetBinContent(i+1) );
			xe.append( h.GetBinWidth(i+1)/2 );
			ye.append( h.GetBinError(i+1) );
			ctr+=1;

		gr = ROOT.TGraphErrors(ctr,x,y,xe,ye);

		return gr;

	def convertHistToGraph2D( self, h, mlo=9999, mhi=-9999, zprimemass=85 ):

		zmlo  = zprimemass*0.90;
		zmhi  = zprimemass*1.10;

		x = array.array('d', [])
		y = array.array('d', [])
		z = array.array('d', [])
		xe = array.array('d', [])
		ye = array.array('d', [])
		ze = array.array('d', [])
		
		ctr = 0;
		for i in range(h.GetNbinsX()):
			for j in range(h.GetNbinsY()):
	
				curx = h.GetXaxis().GetBinCenter(i+1);
				cury = h.GetYaxis().GetBinCenter(j+1);

				xlo = math.log(mlo*mlo/cury);
				xhi = math.log(mhi*mhi/cury);
				xmlo = math.log(zmlo*zmlo/cury);
				xmhi = math.log(zmhi*zmhi/cury);

				# print "windows in rho: ", xlo, xhi, xmlo, xmhi;

				if curx > xlo and curx < xhi: continue;
				if curx > xmlo and curx < xmhi: continue;

				# f(0)   = 3.0;
				# f(100) = 3.2
				# f(300) = 3.6
				# f(400) = 3.8;
				# f(500) = 4.0;
				# f(600) = 4.2;
				# f(x) = 3.0 + 0.2*pT/100
				# if curx > 3.0+ 0.10*cury/100: continue; ## for use with 250 only...
				if curx > 3.0+ 0.2*cury/100: continue; ## default
				# if curx > 3.0+ 0.15*cury/100: continue; ## use for 125, 180
				# if curx > 3.8: continue;
				# if curx > 3.8: continue;
				# if curx > 3.4 + 0.1*cury/100: continue;
				# if curx < 0.2: continue;
				if curx < -0.0 + 0.1*cury/100: continue;
				# if curx > 3.4 + 0.1*cury/100: continue;
				
				x.append( h.GetXaxis().GetBinCenter(i+1) );
				y.append( h.GetYaxis().GetBinCenter(j+1) );
				z.append( h.GetBinContent(i+1,j+1) );
				xe.append( h.GetXaxis().GetBinWidth(i+1)/2 );
				ye.append( h.GetYaxis().GetBinWidth(j+1)/2 );
				ze.append( h.GetBinError(i+1,j+1) );
				ctr+=1;

		print "totalPointsInGraph = ", ctr, len(x);
		gr = ROOT.TGraph2DErrors(ctr,x,y,z,xe,ye,ze);
		print "gr npx = ", gr.GetN();
		return gr;		

	def AdjustHistogramErrors(self,h,herr):

		ax_vals = array.array('d', []);
		axerr_vals = array.array('d', []);
		ay_vals = array.array('d', []);
		ayerr_vals = array.array('d', []);

		for i in range(h.GetNbinsX()):
			
			curerror = h.GetBinError(i+1);
			errup = herr.GetBinContent(i+1) - h.GetBinContent(i+1);
			mccloseerr = 0.00*h.GetBinContent(i+1);
			totalerrorup = math.sqrt(curerror*curerror + errup*errup + mccloseerr*mccloseerr);
			h.SetBinError( i+1, totalerrorup );

			# print "errors = ", curerror, errup

			ax_vals.append( h.GetBinCenter(i+1) )
			axerr_vals.append( h.GetBinWidth(i+1)/2 )
			ay_vals.append( h.GetBinContent(i+1) )
			ayerr_vals.append( totalerrorup )

		gr_pred = ROOT.TGraphErrors(len(ax_vals),ax_vals,ay_vals,axerr_vals,ayerr_vals);
		return gr_pred

