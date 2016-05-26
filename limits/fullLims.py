import ROOT
from optparse import OptionParser
from operator import add
import math
import sys
import time
from array import array
import os
import glob

########################################################################################
########################################################################################
def main():
	
	sigXS = [11.,10.,10.,10.,10.,10.]; # in pb
	masses = [50,100,150,200,250,300];

	idir = "datacards";

	#--------------------------------
	results0 = [];
	results1 = [];
	resultsA = [];	
	for i in range(len(masses)): results0.append( getAsymLimits(idir+'/higgsCombineZprime'+str(masses[i])+'_0.Asymptotic.mH120.root','Zp'+str(masses[i])) );
	for i in range(len(masses)): results1.append( getAsymLimits(idir+'/higgsCombineZprime'+str(masses[i])+'_1.Asymptotic.mH120.root','Zp'+str(masses[i])) );
	for i in range(len(masses)): resultsA.append( getAsymLimits(idir+'/higgsCombineZprime'+str(masses[i])+'_all.Asymptotic.mH120.root','Zp'+str(masses[i])) );



	names   = [];
	l0_obs   = [];
	l0_m2sig = [];
	l0_m1sig = [];
	l0_exp   = [];
	l0_p1sig = [];
	l0_p2sig = [];

	l1_obs   = [];
	l1_m2sig = [];
	l1_m1sig = [];
	l1_exp   = [];
	l1_p1sig = [];
	l1_p2sig = [];

	lA_obs   = [];
	lA_m2sig = [];
	lA_m1sig = [];
	lA_exp   = [];
	lA_p1sig = [];
	lA_p2sig = [];

	ctr = 0
	for r in results0:
		names.append( "Zp"+str(masses[ctr]) );
		l0_m2sig.append(r[1]);
		l0_m1sig.append(r[2]);
		l0_exp.append(r[3]);
		l0_p1sig.append(r[4]);
		l0_p2sig.append(r[5]);
		l0_obs.append(r[0]);
		ctr+=1;

	for r in results1:
		l1_m2sig.append(r[1]);
		l1_m1sig.append(r[2]);
		l1_exp.append(r[3]);
		l1_p1sig.append(r[4]);
		l1_p2sig.append(r[5]);
		l1_obs.append(r[0]);

	for r in resultsA:
		lA_m2sig.append(r[1]);
		lA_m1sig.append(r[2]);
		lA_exp.append(r[3]);
		lA_p1sig.append(r[4]);
		lA_p2sig.append(r[5]);
		lA_obs.append(r[0]);		

	gr0_mu_exp    = makeAGraph( masses, l0_exp, 2, 2 );
	gr0_mu_obs    = makeAGraph( masses, l0_obs, 2, 1 );
	gr0_mu_1sigma = makeAFillGraph( masses, l0_m1sig, l0_p1sig, 0, 3, 1001 );
	gr0_mu_2sigma = makeAFillGraph( masses, l0_m2sig, l0_p2sig, 0, 5, 1001 );

	gr1_mu_exp    = makeAGraph( masses, l1_exp, 4, 2 );
	gr1_mu_obs    = makeAGraph( masses, l1_obs, 4, 1 );
	gr1_mu_1sigma = makeAFillGraph( masses, l1_m1sig, l1_p1sig, 0, 3, 1001 );
	gr1_mu_2sigma = makeAFillGraph( masses, l1_m2sig, l1_p2sig, 0, 5, 1001 );

	grA_mu_exp    = makeAGraph( masses, lA_exp, 1, 2 );
	grA_mu_obs    = makeAGraph( masses, lA_obs, 1, 1 );
	grA_mu_1sigma = makeAFillGraph( masses, lA_m1sig, lA_p1sig, 0, 3, 1001 );
	grA_mu_2sigma = makeAFillGraph( masses, lA_m2sig, lA_p2sig, 0, 5, 1001 );

	# convert to g_B limits! 
	l_obs_gB = [];
	l_exp_gB = [];
	l_m2sig_gB = [];
	l_m1sig_gB = [];
	l_p1sig_gB = [];
	l_p2sig_gB = [];
	for i in range(len(lA_exp)):
		l_obs_gB.append( math.sqrt(lA_obs[i]) );
		l_exp_gB.append( math.sqrt(lA_exp[i]) );
		l_m2sig_gB.append( math.sqrt(lA_m2sig[i]) );
		l_m1sig_gB.append( math.sqrt(lA_m1sig[i]) );
		l_p1sig_gB.append( math.sqrt(lA_p1sig[i]) );
		l_p2sig_gB.append( math.sqrt(lA_p2sig[i]) );

	gr_gB_exp    = makeAGraph( masses, l_exp_gB, 1, 2 );
	gr_gB_obs    = makeAGraph( masses, l_obs_gB, 1, 1 );
	gr_gB_1sigma = makeAFillGraph( masses, l_m1sig_gB, l_p1sig_gB, 0, 3, 1001 );
	gr_gB_2sigma = makeAFillGraph( masses, l_m2sig_gB, l_p2sig_gB, 0, 5, 1001 );

	#--------------------------------
	addFactor=True;
	gr_UA2 = csvToGraph( "externDat/UA2.csv",4,addFactor);
	gr_CDFRun1 = csvToGraph( "externDat/CDF_Run1.csv",2,addFactor );
	gr_CDFRun2 = csvToGraph( "externDat/CDF_Run2.csv",6,addFactor );
	gr_ATLAS = csvToGraph( "externDat/gBMZB_ATLAS_all_fbinv.csv",7,False );
	gr_CMS = csvToGraph( "externDat/CMS_Scouting.csv",8,False );
	#--------------------------------

	#--------------------------------
	#--------------------------------
	#--------------------------------
	#--------------------------------
	#--------------------------------
	#--------------------------------
	#--------------------------------
	# PLOTTING

	txta = ROOT.TLatex(0.11,0.92,"CMS");
	txta.SetNDC();
	txtb = ROOT.TLatex(0.18,0.92,"Preliminary");
	txtb.SetNDC(); txtb.SetTextFont(52);
	txtc = ROOT.TLatex(0.70,0.92,"0.46 fb^{-1} (13 TeV)");
	txtc.SetNDC(); txtc.SetTextFont(42); txtc.SetTextSize(0.04);
	txtd = ROOT.TLatex(0.60,0.80,"g_{B} = 1 or g_{q} = 1/6");
	txtd.SetNDC(); txtd.SetTextFont(42); txtd.SetTextSize(0.04);

	leg = ROOT.TLegend(0.20,0.65,0.4,0.85);
	leg.SetFillStyle(1001);
	leg.SetFillColor(0);    
	leg.SetBorderSize(1);  
	leg.AddEntry(grA_mu_exp,"expected","l")
	leg.AddEntry(grA_mu_obs,"observed","l")
	leg.AddEntry(grA_mu_2sigma,"expected 2#sigma","f")
	leg.AddEntry(grA_mu_1sigma,"expected 1#sigma","f")	
	leg.AddEntry(gr0_mu_obs,"1st jet cat. (obs)","l")	
	leg.AddEntry(gr1_mu_obs,"2nd jet cat. (obs)","l")		

	can_mu = ROOT.TCanvas("can_mu","can_mu",1200,800);
	hrl = can_mu.DrawFrame(40,0,320,15);
	hrl.GetYaxis().SetTitle("#mu = #sigma_{95%CL}/#sigma_{th}");
	hrl.GetXaxis().SetTitle("Z\' mass (GeV)");
	grA_mu_2sigma.Draw('f');
	grA_mu_1sigma.Draw('fsames');
	grA_mu_obs.Draw('lsames');
	grA_mu_exp.Draw('lsames');
	
	gr0_mu_obs.Draw('lsames');
	gr1_mu_obs.Draw('lsames');

	txta.Draw();
	txtb.Draw();
	txtc.Draw();
	txtd.Draw();
	leg.Draw();

	can_mu.SaveAs('limplots/mu.pdf');

	#--------------------------------

	leg2 = ROOT.TLegend(0.20,0.55,0.4,0.85);
	leg2.SetFillStyle(0);
	leg2.SetFillColor(0);    
	leg2.SetBorderSize(0);  
	leg2.AddEntry(gr_gB_exp,"expected","l")
	leg2.AddEntry(gr_gB_obs,"observed","l")
	leg2.AddEntry(gr_gB_2sigma,"expected 2#sigma","f")
	leg2.AddEntry(gr_gB_1sigma,"expected 1#sigma","f")	
	leg2.AddEntry(gr_UA2,"UA2","l")	
	leg2.AddEntry(gr_CDFRun1,"CDF Run 1","l")	
	leg2.AddEntry(gr_CDFRun2,"CDF Run 2","l")	
	leg2.AddEntry(gr_ATLAS,"ATLAS 13 #oplus 20.3  fb^{-1}","l")
	leg2.AddEntry(gr_CMS,"CMS 18.8 fb^{-1} (Data Scouting)","l")

	can_gB = ROOT.TCanvas("can_gB","can_gB",1200,800);
	hrl = can_gB.DrawFrame(40,0,650,4);
	hrl.GetYaxis().SetTitle("g_{B}");
	hrl.GetXaxis().SetTitle("Z\' mass (GeV)");
	
	gr_gB_2sigma.Draw('f');
	gr_gB_1sigma.Draw('fsames');
	gr_gB_obs.Draw('lsames');
	gr_gB_exp.Draw('lsames');

	gr_UA2.Draw("lsames")
	gr_CDFRun1.Draw("lsames")
	gr_CDFRun2.Draw("lsames")
	gr_ATLAS.Draw("lsames")
	gr_CMS.Draw("lsames")

	txta.Draw();
	txtb.Draw();
	txtc.Draw();
	# txtd.Draw();
	leg2.Draw();

	gr_gB_2sigma.GetXaxis().SetMoreLogLabels(True);	
	gr_gB_2sigma.GetXaxis().SetNdivisions(10);	

	# ROOT.gPad.SetLogx();
	can_gB.SaveAs('limplots/gB.pdf');

	#####

	can_gB2 = ROOT.TCanvas("can_gB2","can_gB2",1200,800);
	hrl2 = can_gB.DrawFrame(40,0,1000,4);
	hrl2.GetYaxis().SetTitle("g_{B}");
	hrl2.GetXaxis().SetTitle("Z\' mass (GeV)");
	
	gr_gB_2sigma.Draw('f');
	gr_gB_1sigma.Draw('fsames');
	gr_gB_obs.Draw('lsames');
	gr_gB_exp.Draw('lsames');

	gr_UA2.Draw("lsames")
	gr_CDFRun1.Draw("lsames")
	gr_CDFRun2.Draw("lsames")
	gr_ATLAS.Draw("lsames")
	gr_CMS.Draw("lsames")

	txta.Draw();
	txtb.Draw();
	txtc.Draw();
	# txtd.Draw();
	leg2.Draw();

	gr_gB_2sigma.GetXaxis().SetMoreLogLabels(True);	
	gr_gB_2sigma.GetXaxis().SetNdivisions(10);	

	ROOT.gPad.SetLogx();
	can_gB2.SaveAs('limplots/gB_logx.pdf');


########################################################################################
########################################################################################

def makeAGraph(listx,listy,linecolor = 1, linestyle = 1):

	a_m = array('d', []);
	a_g = array('d', []);

	for i in range(len(listx)):
		a_m.append(listx[i]);
		a_g.append(listy[i]);

	gr = ROOT.TGraph(len(listx),a_m,a_g);

	gr.SetLineColor(linecolor)
	gr.SetLineStyle(linestyle)
	gr.SetLineWidth(2)

	return gr

def makeAFillGraph(listx,listy1,listy2,linecolor = 1, fillcolor = 0, fillstyle = 0):

	a_m = array('d', []);
	a_g = array('d', []);

	for i in range(len(listx)):
		a_m.append(listx[i]);
		a_g.append(listy1[i]);
	
	for i in range(len(listx)-1,-1,-1):
		a_m.append(listx[i]);
		a_g.append(listy2[i]);

	gr = ROOT.TGraph(2*len(listx),a_m,a_g);

	gr.SetLineColor(linecolor)
	gr.SetFillColor(fillcolor)
	gr.SetFillStyle(fillstyle)

	return gr

def csvToGraph(fn, linecolor=1, addFactor=False):

	factor = 1.;
	if addFactor: factor = 6.;

	a_m = array('d', []);
	a_g = array('d', []);

	ifile = open(fn,'r');
	npoints = 0;
	for line in ifile: 
		lline = line.strip().split(',');
		a_m.append(float(lline[0]))
		a_g.append(float(lline[1])*factor)
		npoints += 1;

	gr = ROOT.TGraph(npoints,a_m,a_g);
	gr.SetLineColor(linecolor);
	gr.SetLineWidth(2);

	return gr


########################################################################################
########################################################################################
def getAsymLimits(fname,tag):

    lims = [0]*6;
    
    f = ROOT.TFile(fname);
    t = f.Get("limit");

    #if not t.GetListOfKeys().Contains("limit"): 
    if not t: 
        print "file is corrupted";
        return lims

    entries = t.GetEntries();
    for i in range(entries):

        t.GetEntry(i);
        t_quantileExpected = t.quantileExpected;
        t_limit = t.limit;

        #print "limit: ", t_limit, ", quantileExpected: ",t_quantileExpected;
        
        if t_quantileExpected == -1.: lims[0] = t_limit;
        elif t_quantileExpected >= 0.024 and t_quantileExpected <= 0.026: lims[1] = t_limit;
        elif t_quantileExpected >= 0.15 and t_quantileExpected <= 0.17: lims[2] = t_limit;            
        elif t_quantileExpected == 0.5: lims[3] = t_limit;            
        elif t_quantileExpected >= 0.83 and t_quantileExpected <= 0.85: lims[4] = t_limit;
        elif t_quantileExpected >= 0.974 and t_quantileExpected <= 0.976: lims[5] = t_limit;
        else: print "Unknown quantile!"

    print "lims = ", lims
    return lims;

if __name__ == '__main__':

	main();
