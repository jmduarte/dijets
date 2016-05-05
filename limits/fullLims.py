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
	results = [];
	for i in range(len(masses)): results.append( getAsymLimits(idir+'/higgsCombineZprime'+str(masses[i])+'.Asymptotic.mH120.root','Zp'+str(masses[i])) );

	names   = [];
	l_obs   = [];
	l_m2sig = [];
	l_m1sig = [];
	l_exp   = [];
	l_p1sig = [];
	l_p2sig = [];
	ctr = 0
	for r in results:
		names.append( "Zp"+str(masses[ctr]) );
		l_m2sig.append(r[1]);
		l_m1sig.append(r[2]);
		l_exp.append(r[3]);
		l_p1sig.append(r[4]);
		l_p2sig.append(r[5]);
		l_obs.append(r[0]);
		ctr+=1;

	gr_mu_exp    = makeAGraph( masses, l_exp, 1, 2 );
	gr_mu_obs    = makeAGraph( masses, l_obs, 1, 1 );
	gr_mu_1sigma = makeAFillGraph( masses, l_m1sig, l_p1sig, 0, 3, 1001 );
	gr_mu_2sigma = makeAFillGraph( masses, l_m2sig, l_p2sig, 0, 5, 1001 );

	# convert to g_B limits! 
	l_obs_gB = [];
	l_exp_gB = [];
	l_m2sig_gB = [];
	l_m1sig_gB = [];
	l_p1sig_gB = [];
	l_p2sig_gB = [];
	for i in range(len(l_exp)):
		l_obs_gB.append( math.sqrt(l_obs[i]) );
		l_exp_gB.append( math.sqrt(l_exp[i]) );
		l_m2sig_gB.append( math.sqrt(l_m2sig[i]) );
		l_m1sig_gB.append( math.sqrt(l_m1sig[i]) );
		l_p1sig_gB.append( math.sqrt(l_p1sig[i]) );
		l_p2sig_gB.append( math.sqrt(l_p2sig[i]) );

	gr_gB_exp    = makeAGraph( masses, l_exp_gB, 1, 2 );
	gr_gB_obs    = makeAGraph( masses, l_obs_gB, 1, 1 );
	gr_gB_1sigma = makeAFillGraph( masses, l_m1sig_gB, l_p1sig_gB, 0, 3, 1001 );
	gr_gB_2sigma = makeAFillGraph( masses, l_m2sig_gB, l_p2sig_gB, 0, 5, 1001 );

	#--------------------------------
	gr_UA2 = csvToGraph( "externDat/UA2.csv",4 );
	gr_CDFRun1 = csvToGraph( "externDat/CDF_Run1.csv",2 );
	gr_CDFRun2 = csvToGraph( "externDat/CDF_Run2.csv",6 );
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
	leg.AddEntry(gr_mu_exp,"expected","l")
	leg.AddEntry(gr_mu_obs,"observed","l")
	leg.AddEntry(gr_mu_2sigma,"expected 2#sigma","f")
	leg.AddEntry(gr_mu_1sigma,"expected 1#sigma","f")	

	can_mu = ROOT.TCanvas("can_mu","can_mu",1200,800);
	hrl = can_mu.DrawFrame(40,0,320,15);
	hrl.GetYaxis().SetTitle("#mu = #sigma_{95%CL}/#sigma_{th}");
	hrl.GetXaxis().SetTitle("Z\' mass (GeV)");
	gr_mu_2sigma.Draw('f');
	gr_mu_1sigma.Draw('fsames');
	gr_mu_obs.Draw('lsames');
	gr_mu_exp.Draw('lsames');

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
	leg2.AddEntry(gr_mu_exp,"expected","l")
	leg2.AddEntry(gr_mu_obs,"observed","l")
	leg2.AddEntry(gr_mu_2sigma,"expected 2#sigma","f")
	leg2.AddEntry(gr_mu_1sigma,"expected 1#sigma","f")	
	leg2.AddEntry(gr_UA2,"UA2","l")	
	leg2.AddEntry(gr_CDFRun1,"CDF Run 1","l")	
	leg2.AddEntry(gr_CDFRun2,"CDF Run 2","l")	

	can_gB = ROOT.TCanvas("can_gB","can_gB",1200,800);
	hrl = can_gB.DrawFrame(40,0,800,4);
	hrl.GetYaxis().SetTitle("g_{B}");
	hrl.GetXaxis().SetTitle("Z\' mass (GeV)");
	
	gr_gB_2sigma.Draw('f');
	gr_gB_1sigma.Draw('fsames');
	gr_gB_obs.Draw('lsames');
	gr_gB_exp.Draw('lsames');

	gr_UA2.Draw("lsames")
	gr_CDFRun1.Draw("lsames")
	gr_CDFRun2.Draw("lsames")

	txta.Draw();
	txtb.Draw();
	txtc.Draw();
	# txtd.Draw();
	leg2.Draw();

	gr_gB_2sigma.GetXaxis().SetMoreLogLabels(True);	
	gr_gB_2sigma.GetXaxis().SetNdivisions(10);	

	ROOT.gPad.SetLogx();
	can_gB.SaveAs('limplots/gB.pdf');

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

def csvToGraph(fn, linecolor=1):

	a_m = array('d', []);
	a_g = array('d', []);

	ifile = open(fn,'r');
	npoints = 0;
	for line in ifile: 
		lline = line.strip().split(',');
		a_m.append(float(lline[0]))
		a_g.append(float(lline[1])*6.)
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