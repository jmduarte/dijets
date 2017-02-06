import ROOT
from ROOT import TFile, TTree, TChain, gPad, gDirectory
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option("--lumi", dest="lumi", default = 4,help="mass of LSP", metavar="MLSP")
parser.add_option("--binning",dest="binning",default="RA2bBins",help="Select binning to be used: Classic, SMJ, extSMJ", metavar="binning")

(options, args) = parser.parse_args()

#
import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadTopMargin(0.06);
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPadRightMargin(0.05);
ROOT.gStyle.SetPalette(1);
ROOT.gStyle.SetPaintTextFormat("1.1f");
ROOT.gStyle.SetOptFit(0000);

#########################################################################################################

def makeCanvas(hists,normalize=False,legnames=None):

	color = [1,2,4,6,7,8,3,5]

	c = ROOT.TCanvas("c"+hists[0].GetName(),"c"+hists[0].GetName(),900,800);

	max = -999;

	txta = ROOT.TLatex(0.20,0.85,"CMS");
	txta.SetNDC(); txta.SetTextSize( 0.065 );
	txtb = ROOT.TLatex(0.32,0.85,"Simulation");
	txtb.SetNDC(); txtb.SetTextFont(52); txtb.SetTextSize( 0.065 );
	txtc = ROOT.TLatex(0.20,0.79,"SM QCD Multijet");
	txtc.SetNDC(); txtc.SetTextFont(42); txtb.SetTextSize( 0.055 );

	leg = ROOT.TLegend(0.5,0.20,0.8,0.45);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	leg.SetTextFont(42);
	leg.SetTextSize(0.045);
	if legnames != None: 
		for i in range(len(legnames)):
			leg.AddEntry(hists[i],legnames[i],'pe') 	

	for i in range(len(hists)):
		hists[i].SetLineColor(color[i])
		hists[i].SetMarkerColor(color[i])
		hists[i].SetMarkerStyle(24)
		if hists[i].GetMaximum() > max: 
			max = hists[i].GetMaximum();
			hists[0].SetMaximum(max*1.25);
		if normalize and hists[i].Integral() > 0: hists[i].Scale(1./hists[i].Integral())
		option = 'pe';
		if i > 0: option+='sames';
		hists[i].Draw(option);

	txta.Draw();
	txtb.Draw();
	txtc.Draw();
	if legnames != None: leg.Draw();
	c.SaveAs("plots_qcd/"+hists[0].GetName()+".pdf")
	ROOT.gPad.SetLogy();
	c.SaveAs("plots_qcd/"+hists[0].GetName()+"_log.pdf")

if __name__ == '__main__':


	f = ROOT.TFile("../sklimming/sklim-v0-Jun16/QCD.root");
	t = f.Get("otree");

	h_pt = ROOT.TH1F("h_pt",";pt;N",50,0,3000);
	h_msd = ROOT.TH1F("h_msd",";mass (soft drop);N",50,0,300);
	h_msd_postCut = ROOT.TH1F("h_msd_postCut",";mass (soft drop);N",50,0,300);
	h_msd_antiCut = ROOT.TH1F("h_msd_antiCut",";mass (soft drop);N",50,0,300);

	h_t21 = ROOT.TH1F("h_t21",";#tau_{21};N",50,0,1);
	h_t21P = ROOT.TH1F("h_t21P",";#tau_{21}';N",50,0,1.2);

	n_ptbins = 5; 
	ptlo = 500; 
	ptwidth = 150;
	pthi = 1250;
	t21PCut = 0.45;

	hn_rho_V_t21 = [];
	hn_rhP_V_t21 = [];
	hn_rhP_V_t21P = [];
	hn_logm_V_t21P = [];
	hn_msd = [];
	hn_t21 = [];
	hn_t21P = [];
	hn_rho = [];
	hn_rhP = [];

	hn_msd_postCut = [];
	hn_rho_postCut = [];
	hn_rhP_postCut = [];
	hn_msd_antiCut = [];
	hn_rho_antiCut = [];	
	hn_rhP_antiCut = [];	

	hn_rhP_passfail = [];	

	for i in range(n_ptbins):
		hn_rho_V_t21.append( ROOT.TH2F("hn_rho_V_t21"+str(i),";#rho;t21",20,-15,0,20,0,1) );
		hn_rhP_V_t21.append( ROOT.TH2F("hn_rhP_V_t21"+str(i),";#rho';t21",20,-2,7,20,0,1) );
		hn_rhP_V_t21P.append( ROOT.TH2F("hn_rhP_V_t21P"+str(i),";#rho';t21'",20,-2,7,20,0,1.2) );
		hn_logm_V_t21P.append( ROOT.TH2F("hn_logm_V_t21P"+str(i),";log(m^{2});t21'",20,0,10,20,0,1.2) );
		
		hn_msd.append( ROOT.TH1F("hn_msd"+str(i),";mass;N",25,0,500) );
		hn_t21.append( ROOT.TH1F("hn_t21"+str(i),";#tau_{2}/#tau_{1};N",25,0,1) );
		hn_t21P.append( ROOT.TH1F("hn_t21P"+str(i),";#tau_{21}';N",25,-0.2,1.2) );		
		hn_rho.append( ROOT.TH1F("hn_rho"+str(i),";#rho;N",25,-15,0) );
		hn_rhP.append( ROOT.TH1F("hn_rhP"+str(i),";#rho';N",10,0.3,4) );

		hn_msd_postCut.append( ROOT.TH1F("hn_msd_postCut"+str(i),";mass;N",25,0,500) );
		hn_rho_postCut.append( ROOT.TH1F("hn_rho_postCut"+str(i),";#rho;N",25,-15,0) );
		hn_rhP_postCut.append( ROOT.TH1F("hn_rhP_postCut"+str(i),";#rho';N",25,-2,7) );
		hn_msd_antiCut.append( ROOT.TH1F("hn_msd_antiCut"+str(i),";mass;N",25,0,500) );
		hn_rho_antiCut.append( ROOT.TH1F("hn_rho_antiCut"+str(i),";#rho;N",25,-15,0) );
		hn_rhP_antiCut.append( ROOT.TH1F("hn_rhP_antiCut"+str(i),";#rho';N",10,0.3,4) );

		hn_rhP_passfail.append( ROOT.TH1F("hn_rhP_passfail"+str(i),";#rho';N",10,0.3,4) );

	nent = t.GetEntries();
	for i in range(t.GetEntries()):


		if(i % (1 * nent/100) == 0):
			sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
			sys.stdout.flush();

		if i % 100 != 0: continue;

		t.GetEntry(i);

		jpt = getattr(t,"bst8_PUPPIjet0_pt");
		# jpt = jp4.Pt();
		jmsd = getattr(t,"bst8_PUPPIjet0_msd");		
		weight = getattr(t,"scale1fb");
		h_pt.Fill( jpt, weight );

		if jpt > ptlo and jmsd > 0:

			jt21 = getattr(t,"bst8_PUPPIjet0_tau21");
			h_msd.Fill( jmsd, weight );
			h_t21.Fill( jt21, weight );

			rho = math.log(jmsd*jmsd/jpt/jpt);
			rhP = math.log(jmsd*jmsd/jpt);
			jt21P = jt21 + 0.063*rhP;
		 	h_t21P.Fill( jt21P, weight );

			for j in range(n_ptbins):
				if jpt > (ptlo + j*float(ptwidth)) and jpt < (ptlo + (j+1)*float(ptwidth)):
					hn_rho_V_t21[j].Fill(rho,jt21,weight);
					hn_rhP_V_t21[j].Fill(rhP,jt21,weight);
					hn_rhP_V_t21P[j].Fill(rhP,jt21P,weight);
					hn_logm_V_t21P[j].Fill(math.log(jmsd*jmsd),jt21P,weight);

					hn_msd[j].Fill(jmsd,weight);
					hn_t21[j].Fill(jt21,weight);
					hn_t21P[j].Fill(jt21P,weight);
					hn_rho[j].Fill(rho,weight);
					hn_rhP[j].Fill(rhP,weight);

					if jt21P < t21PCut and rhP > 0.: 
						hn_msd_postCut[j].Fill(jmsd,weight);
						hn_rho_postCut[j].Fill(rho,weight);
						hn_rhP_postCut[j].Fill(rhP,weight);
						hn_rhP_passfail[j].Fill(rhP,weight);

					if jt21P > t21PCut and rhP > 0.: 
						hn_msd_antiCut[j].Fill(jmsd,weight);
						hn_rho_antiCut[j].Fill(rho,weight);
						hn_rhP_antiCut[j].Fill(rhP,weight);

			if jt21P < t21PCut and rhP > 0.3: h_msd_postCut.Fill( jmsd, weight );
			if jt21P > t21PCut and rhP > 0.3: h_msd_antiCut.Fill( jmsd, weight );

		if i > 1e11: break;			      


	makeCanvas([h_pt]);
	makeCanvas([h_msd]);
	makeCanvas([h_msd_postCut]);
	makeCanvas([h_msd_antiCut]);

	makeCanvas([h_t21]);
	makeCanvas([h_t21P]);

	hn_rho_V_t21_prof = [];
	hn_rhP_V_t21_prof = [];	
	hn_rhP_V_t21P_prof = [];	
	hn_logm_V_t21P_prof = [];	
	for i in range(len(hn_rho_V_t21)): hn_rho_V_t21_prof.append( hn_rho_V_t21[i].ProfileX() );
	for i in range(len(hn_rhP_V_t21)): hn_rhP_V_t21_prof.append( hn_rhP_V_t21[i].ProfileX() );
	for i in range(len(hn_rhP_V_t21P)): hn_rhP_V_t21P_prof.append( hn_rhP_V_t21P[i].ProfileX() );
	for i in range(len(hn_logm_V_t21P)): hn_logm_V_t21P_prof.append( hn_logm_V_t21P[i].ProfileX() );
	hn_rho_V_t21_prof[0].GetYaxis().SetRangeUser(0.0,1.0);
	hn_rhP_V_t21_prof[0].GetYaxis().SetRangeUser(0.0,1.0);
	hn_rhP_V_t21P_prof[0].GetYaxis().SetRangeUser(0.2,0.8);
	hn_logm_V_t21P_prof[0].GetYaxis().SetRangeUser(0.0,1.0);
	hn_rho_V_t21_prof[0].SetTitle(";#rho;< #tau_{21 }>");
	hn_rhP_V_t21_prof[0].SetTitle(";#rho^{DDT};< #tau_{21} >");
	hn_rhP_V_t21P_prof[0].SetTitle(";#rho^{DDT};< #tau_{21}^{DDT} >");

	# fitting
	for i in range(n_ptbins):
		mvs_fit = ROOT.TF1("P1","pol1",1,4);
		hn_rhP_V_t21_prof[i].Fit(mvs_fit,"RQ")
		print "fit parameters = ", mvs_fit.GetParameter(0), mvs_fit.GetParameter(1);
		intercept = mvs_fit.GetParameter(0);
		slope = mvs_fit.GetParameter(1);  

	normalize = False;

	lnames = [];
	lnames.append( "p_{T} = 500-650 GeV" );
	lnames.append( "p_{T} = 650-800 GeV" );
	lnames.append( "p_{T} = 800-950 GeV" );
	lnames.append( "p_{T} = 950-1100 GeV" );
	lnames.append( "p_{T} = 1100-1250 GeV" );

	makeCanvas( hn_rho_V_t21_prof, normalize, lnames );
	makeCanvas( hn_rhP_V_t21_prof, normalize, lnames );
	makeCanvas( hn_rhP_V_t21P_prof, normalize, lnames );
	makeCanvas( hn_logm_V_t21P_prof, normalize, lnames );
	makeCanvas( hn_msd,normalize );
	makeCanvas( hn_t21,normalize );
	makeCanvas( hn_t21P,normalize );
	makeCanvas( hn_rho,normalize );
	makeCanvas( hn_rhP,normalize );	

	makeCanvas( hn_msd_postCut,normalize );	
	makeCanvas( hn_rho_postCut,normalize );	
	makeCanvas( hn_rhP_postCut,normalize );	
	makeCanvas( hn_msd_antiCut,normalize );	
	makeCanvas( hn_rho_antiCut,normalize );	
	makeCanvas( hn_rhP_antiCut,normalize );

	for i in range(len(hn_rhP_passfail)): hn_rhP_passfail[i].Divide(hn_rhP_antiCut[i])
	makeCanvas ( hn_rhP_passfail );

	###################################################################################################
	pTbins = array.array('d', []);
	for i in range(n_ptbins): pTbins.append( ptlo + i*ptwidth + (ptwidth/2.) );

	PFRatioPerBin = [];
	for j in range(hn_rhP_passfail[0].GetNbinsX()): 
		PFRatioPerBin.append( array.array('d', []) );

	for i in range(len(hn_rhP_passfail)):
		for j in range(hn_rhP_passfail[0].GetNbinsX()):
			PFRatioPerBin[j].append( hn_rhP_passfail[i].GetBinContent(j+1) );

	tgs = [];
	for a in PFRatioPerBin: tgs.append( ROOT.TGraph( len(pTbins), pTbins, a ) ); 
	
	print tgs

	c2 = ROOT.TCanvas("c2","c2",1200,800);
	hr2 = c2.DrawFrame(500,0,1300,0.2);
	hr2.GetYaxis().SetTitle("npass/nfail in a given bin");
	hr2.GetXaxis().SetTitle("pT");
	c2.SetGrid(); 

	colors = [1,2,4,6,7,1,2,4,6,7]
	styles = [1,1,1,1,1,2,2,2,2,2]
	gctr = 0;
	for g in tgs: 
		curc = gctr % 3 + 1
		g.SetLineColor(colors[gctr]);
		g.SetLineWidth(2);
		g.SetLineStyle(styles[gctr])
		g.SetMarkerStyle(20);
		g.SetMarkerColor(colors[gctr]);
		g.Draw("lp");
		gctr += 1;

	c2.SaveAs("plots_qcd/testg.pdf");




