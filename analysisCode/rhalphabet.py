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

parser.add_option("--lumi", dest="lumi", default = 0.44,help="mass of LSP", metavar="MLSP")
parser.add_option("--rholo", dest="rholo", default = 0.,help="mass of LSP", metavar="MLSP")
parser.add_option("--rhohi", dest="rhohi", default = 6.,help="mass of LSP", metavar="MLSP")
parser.add_option("--DDTcut", dest="DDTcut", default = 0.38,help="mass of LSP", metavar="MLSP")
parser.add_option('--qcdClosure', action='store_true', dest='qcdClosure', default=False, help='go!')

(options, args) = parser.parse_args()

class rhalphabet:

	################################################################################################
	# init
	def __init__( self , filename, lumi, name, scaleFactor = 1, extractTFs = True):

		## fixed parameters 
		# self._ptbins = 5;
		# self._ptlo = 500;
		# self._pthi = 1000;
		self._ptbins = 5;
		self._ptlo = 350;
		self._pthi = 850;

		self._filename = filename;
		self._tf = ROOT.TFile( self._filename );
		self._tt = self._tf.Get("otree");
		self._sf = scaleFactor;
		self._name = name;

		tfopt = "READ";
		if extractTFs: tfopt = "RECREATE"
		self._storageFile = ROOT.TFile("plots/rhalphabet/storedHistos.root",tfopt);
		self.hys = None;
		self.TF_pafa = None;

		if extractTFs: 
			self.FillTFMap();
			self._storageFile.cd();
			self.hys = makeCanvas2D( self.TF_pafa,"map","plots/rhalphabet" );
			for h in self.hys: h.Write();
			self.TF_pafa.Write();
		else:
			self.hys = [];
			for i in range(self._ptbins): self.hys.append( self._storageFile.Get("hys"+str(i)) );
			self.TF_pafa = self._storageFile.Get("TF_pafa");			

		self.do2DFit();
		# self.do2DRooFit();

		# else:
		# 	# get stored information

		self._storageFile.Close();

	################################################################################################
	# fill map
	def FillTFMap( self ):

		DDTCUT = float(options.DDTcut);

		self.TF_fail = ROOT.TH2F("TF_fail",";rho;pT",30,float(options.rholo),float(options.rhohi),self._ptbins,self._ptlo,self._pthi)
		self.TF_pafa = ROOT.TH2F("TF_pafa",";rho;pT",30,float(options.rholo),float(options.rhohi),self._ptbins,self._ptlo,self._pthi)	
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
			if i % int(self._sf) != 0: continue;

			# cutting
			jpt = getattr(t,"bst8_PUPPIjet0_pt");
			jmsd = getattr(t,"bst8_PUPPIjet0_msd");		
			weight = getattr(t,"scale1fb");

			if jmsd > 0:

				jt21 = getattr(t,"bst8_PUPPIjet0_tau21");
				rhP = math.log(jmsd*jmsd/jpt);			
				jt21P = jt21 + 0.063*rhP;

				if jt21P < DDTCUT: self.TF_pafa.Fill( rhP, jpt, weight );
				else: self.TF_fail.Fill( rhP, jpt, weight );

		print "\n";

		#####!!!! this uncertinaty may be wrong, binomial instead of poisson
		self.TF_pafa.Divide( self.TF_fail );

	################################################################################################
	# fit 2D space
	def do2DFit( self ):
		THEFITRHOMAX = 3.8;

		# turn the histograms into graphs, that way exclude the various mass windows
		grs = [];
		for i in range(len(self.hys)): 
			curpt = self.TF_pafa.GetYaxis().GetBinCenter(i+1);
			xlo = math.log(70.*70./curpt);
			xhi = math.log(90.*90./curpt);
			grs.append( self.convertHistToGraph( self.hys[i], xlo, xhi ) );
			# grs.append( self.convertHistToGraph( self.hys[i] ) );

		# do the fits
		self.theRhoFits = [];
		ctr = 0
		for g in grs:
			# self.theRhoFits.append( QuadraticFit([0.01,6e-5,0.01],float(options.rholo),float(options.rhohi),"linFit1", "EMRNS") );
			# g.Fit(self.theRhoFits[ctr].fit,"RQ")
			# fitter = ROOT.TVirtualFitter.GetFitter()
			# self.theRhoFits[ctr].Converter(fitter)

			# self.theRhoFits.append( ROOT.TF1("rhofit"+str(ctr),'([0]+[1]*x+[2]*x*x)',float(options.rholo),float(options.rhohi)) ); 	
			self.theRhoFits.append( ROOT.TF1("rhofit"+str(ctr),'([0]+[1]*x+[2]*x*x)',float(options.rholo),THEFITRHOMAX ) ); 	
			
			self.theRhoFits[ctr].SetParameter(0,0.05);
			self.theRhoFits[ctr].SetParameter(1,0.005);					
			self.theRhoFits[ctr].SetParameter(2,0.);						

			self.theRhoFits[ctr].SetParLimits(0,0.,0.1);
			self.theRhoFits[ctr].SetParLimits(1,-0.1,0.1);					
			self.theRhoFits[ctr].SetParLimits(2,-0.1,0.1);					

			g.Fit(self.theRhoFits[ctr],"R");
			par0 = self.theRhoFits[ctr].GetParameter(0);
			par1 = self.theRhoFits[ctr].GetParameter(1);
			par2 = self.theRhoFits[ctr].GetParameter(2);

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

			curcan.SaveAs("plots/rhalphabet/map_rhodependence_bin"+str(i)+".pdf")

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
		
		self.thePtFits = [];
		ctr = 0;
		for h in hpt_pars:
			self.thePtFits.append( QuadraticFit([0.01,6e-5,0.01],self._ptlo,self._pthi,"linFit2", "EMRNS") );
			curfitresult = h.Fit(self.thePtFits[ctr].fit,"RS")
			fitter = ROOT.TVirtualFitter.GetFitter()
			self.thePtFits[ctr].Converter(fitter)

			ctr+=1

		for i in range(len(hpt_pars)):
			curcan2 = ROOT.TCanvas("curcan2"+str(i),"curcan2"+str(i),1200,800);
			# curhr2 = curcan.DrawFrame(ptlo[0],0,pthi[nptbins-1],0.2);
			# curhr2.GetYaxis().SetTitle("npass/nfail in a #rho^{DDT} bin");
			# curhr2.GetXaxis().SetTitle("pT");	
			curcan2.SetGrid();
			hpt_pars[i].Draw("pe");
			self.thePtFits[i].fit.SetLineColor(2);	
			self.thePtFits[i].ErrUp.SetLineColor(2);	
			self.thePtFits[i].ErrDn.SetLineColor(2);	
			self.thePtFits[i].ErrUp.SetLineStyle(2);	
			self.thePtFits[i].ErrDn.SetLineStyle(2);	
			self.thePtFits[i].fit.Draw("sames");	
			self.thePtFits[i].ErrUp.Draw("sames");	
			self.thePtFits[i].ErrDn.Draw("sames");	
			curcan2.SaveAs("plots/rhalphabet/map_ptdependence_par"+str(i)+".pdf")	

		self.effPlane = ROOT.TF2("TransferPlane", "(([0]+ [1]*y + [2]*y*y) + ([3]+ [4]*y + [5]*y*y)*x + ([6]+ [7]*y + [8]*y*y)*x*x)",float(options.rholo),THEFITRHOMAX,self._ptlo,self._pthi);

		self.effPlane.SetParameter(0,self.thePtFits[0].fit.GetParameter(0));
		self.effPlane.SetParError(0,self.thePtFits[0].fit.GetParError(0));
		self.effPlane.SetParameter(1,self.thePtFits[0].fit.GetParameter(1));
		self.effPlane.SetParError(1,self.thePtFits[0].fit.GetParError(1));
		self.effPlane.SetParameter(2,self.thePtFits[0].fit.GetParameter(2));
		self.effPlane.SetParError(2,self.thePtFits[0].fit.GetParError(2));

		self.effPlane.SetParameter(3,self.thePtFits[1].fit.GetParameter(0));
		self.effPlane.SetParError(3,self.thePtFits[1].fit.GetParError(0));
		self.effPlane.SetParameter(4,self.thePtFits[1].fit.GetParameter(1));
		self.effPlane.SetParError(4,self.thePtFits[1].fit.GetParError(1));
		self.effPlane.SetParameter(5,self.thePtFits[1].fit.GetParameter(2));
		self.effPlane.SetParError(5,self.thePtFits[1].fit.GetParError(2));

		self.effPlane.SetParameter(6,self.thePtFits[2].fit.GetParameter(0));
		self.effPlane.SetParError(6,self.thePtFits[2].fit.GetParError(0));
		self.effPlane.SetParameter(7,self.thePtFits[2].fit.GetParameter(1));
		self.effPlane.SetParError(7,self.thePtFits[2].fit.GetParError(1));
		self.effPlane.SetParameter(8,self.thePtFits[2].fit.GetParameter(2));
		self.effPlane.SetParError(8,self.thePtFits[2].fit.GetParError(2));

		self.effPlane_prefit = self.effPlane.Clone();

		cplane = ROOT.TCanvas("cplane","cplane",1000,800);
		self.effPlane.SetTitle(";#rho^{DDT};pT;N_{pass}/N_{fail}")
		self.effPlane.Draw("surf1");
		# ROOT.gPad.SetPhi(0.6);
		cplane.SaveAs("plots/rhalphabet/map_parameterized.pdf")			

		# need to update this to the TGraph2D with the mass windows missing?
		self.TF_pafa_gr2D = self.convertHistToGraph2D(self.TF_pafa,70,90);
		self.theFitResult = self.TF_pafa_gr2D.Fit(self.effPlane,"RS");
		# curpoint = array.array('d', []);
		# curpoint.append( 2 );
		# curpoint.append( 1400 );
		# curerr = array.array('d', []);
		# curerr.append(0.)
		# r.GetConfidenceIntervals(1,1,1,curpoint,curerr, 0.683, False);
		# print curerr

		cplane2 = ROOT.TCanvas("cplane2","cplane2",1000,800);
		self.effPlane.SetTitle(";#rho^{DDT};pT;N_{pass}/N_{fail}")
		self.effPlane.Draw("surf1");
		# ROOT.gPad.SetPhi(0.6);
		cplane2.SaveAs("plots/rhalphabet/map_parameterized_postfit.pdf")	

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
		cplane3.SaveAs("plots/rhalphabet/map_diff.pdf")
		print self.TF_pafa	

	################################################################################################
	# predicted distribution from the given input file
	def GetPredictedDistributions(self, filename, lumi=1, scaleFactor=1, isData=False):
		
		DDTCUT = float(options.DDTcut);

		self.pred_fn = filename;
		self.pred_tf = ROOT.TFile( self.pred_fn );
		self.pred_tt = self.pred_tf.Get("otree");
		self.pred_scaleFactor = scaleFactor;

		# Define histogram
		self.hpred_jetmsd         = ROOT.TH1F("hpred_jetmsd"+self._name,"; soft drop mass (GeV);", 60, 30, 330);
		self.hpred_jetmsd_errup   = ROOT.TH1F("hpred_jetmsd_errup"+self._name,"; soft drop mass (GeV);", 60, 30, 330);
		self.hpred_jetmsd_failcut = ROOT.TH1F("hpred_jetmsd_failcut"+self._name,"; soft drop mass (GeV);", 60, 30, 330);
		self.hpred_rhoDDT         = ROOT.TH1F("hpred_rhoDDT"+self._name,"; soft drop mass (GeV);", 20, float(options.rholo),float(options.rhohi));
		self.hpred_rhoDDT_errup    = ROOT.TH1F("hpred_rhoDDT_errup"+self._name,"; soft drop mass (GeV);", 20, float(options.rholo),float(options.rhohi));
		self.hpred_rhoDDT_failcut = ROOT.TH1F("hpred_rhoDDT_failcut"+self._name,"; soft drop mass (GeV);", 20, float(options.rholo),float(options.rhohi));
		
		self.hpred_jetmsd.Sumw2();
		self.hpred_jetmsd_errup.Sumw2();
		self.hpred_jetmsd_failcut.Sumw2();
		self.hpred_rhoDDT.Sumw2();
		self.hpred_rhoDDT_errup.Sumw2();
		self.hpred_rhoDDT_failcut.Sumw2();

		# looping
		nent = self.pred_tt.GetEntries();
		for i in range(self.pred_tt.GetEntries()):

			# preamble
			self.pred_tt.GetEntry(i);
			if(i % (1 * nent/100) == 0):
				sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
				sys.stdout.flush();
			
			if i % self.pred_scaleFactor != 0: continue;		

			jpt = getattr(self.pred_tt,"bst8_PUPPIjet0_pt");
			jeta = getattr(self.pred_tt,"bst8_PUPPIjet0_eta");
			jmsd = getattr(self.pred_tt,"bst8_PUPPIjet0_msd");		
			if jmsd == 0.: jmsd = 0.01;
			weight = float(self.pred_scaleFactor)*float(lumi)*getattr(self.pred_tt,"scale1fb");
			if isData: weight = 1;

			if jpt < 500: continue;

			jt21 = getattr(self.pred_tt,"bst8_PUPPIjet0_tau21");
			rhP = math.log(jmsd*jmsd/jpt);                  
			jt21P = jt21 + 0.063*rhP;

			if rhP > 0 and jt21P > DDTCUT:
				self.hpred_jetmsd_failcut.Fill( jmsd, weight );
				self.hpred_rhoDDT_failcut.Fill( rhP, weight );
 				curval,curerr = self.predictedVal(rhP,jpt);
				self.hpred_jetmsd.Fill( jmsd, weight*curval );
				self.hpred_rhoDDT.Fill( rhP, weight*curval );
				self.hpred_jetmsd_errup.Fill( jmsd, weight*(curval+curerr) );
				self.hpred_rhoDDT_errup.Fill( rhP, weight*(curval+curerr) );

		print "\n"

		# adjust histogram errors
		self.grpred_jetmsd = self.AdjustHistogramErrors(self.hpred_jetmsd,self.hpred_jetmsd_errup);
		self.grpred_rhoDDT = self.AdjustHistogramErrors(self.hpred_rhoDDT,self.hpred_rhoDDT_errup);

	def predictedVal(self,xval,yval):

		curpoint = array.array('d', []);
		curpoint.append( xval );
		curpoint.append( yval );
		curerr = array.array('d', []);
		curerr.append(0.)
		self.theFitResult.GetConfidenceIntervals(1,1,1,curpoint,curerr, 0.683, False);	
		
		curval = self.effPlane.Eval( xval, yval );
		curval2 = self.effPlane_prefit.Eval( xval, yval );
		# print curval,curval2,curerr[0]

		# storageFile2 = ROOT.TFile("plots/rhalphabet/storedHistos.root","READ");
		# tfpafa = storageFile2.Get("TF_pafa");
		# # print self.TF_pafa
		# # print xval,yval,self.TF_pafa.GetTitle();#,self.TF_pafa.GetXaxis().FindBin(xval), self.TF_pafa.GetYaxis().FindBin(yval)
		# curval = tfpafa.GetBinContent( tfpafa.GetXaxis().FindBin(xval), tfpafa.GetYaxis().FindBin(yval) );

		return curval, curerr[0];

	################################################################################################
	# helpers
	def convertHistToGraph( self, h, xlo=9999, xhi=-9999 ):

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

	def convertHistToGraph2D( self, h, mlo=9999, mhi=-9999 ):

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
				if curx > xlo and curx < xhi: continue;
				
				x.append( h.GetXaxis().GetBinCenter(i+1) );
				y.append( h.GetYaxis().GetBinCenter(j+1) );
				z.append( h.GetBinContent(i+1,j+1) );
				xe.append( h.GetXaxis().GetBinWidth(i+1)/2 );
				ye.append( h.GetYaxis().GetBinWidth(j+1)/2 );
				ze.append( h.GetBinError(i+1,j+1) );
				ctr+=1;

		gr = ROOT.TGraph2DErrors(ctr,x,y,z,xe,ye,ze);

		return gr;		

	def AdjustHistogramErrors(self,h,herr):

		ax_vals = array.array('d', []);
		axerr_vals = array.array('d', []);
		ay_vals = array.array('d', []);
		ayerr_vals = array.array('d', []);

		for i in range(h.GetNbinsX()):
			
			curerror = h.GetBinError(i+1);
			errup = herr.GetBinContent(i+1) - h.GetBinContent(i+1);
			totalerrorup = math.sqrt(curerror*curerror + errup*errup);
			h.SetBinError( i+1, totalerrorup );

			ax_vals.append( h.GetBinCenter(i+1) )
			axerr_vals.append( h.GetBinWidth(i+1)/2 )
			ay_vals.append( h.GetBinContent(i+1) )
			ayerr_vals.append( totalerrorup )

		gr_pred = ROOT.TGraphErrors(len(ax_vals),ax_vals,ay_vals,axerr_vals,ayerr_vals);
		return gr_pred

