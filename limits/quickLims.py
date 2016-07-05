import ROOT
from optparse import OptionParser
from operator import add
import math
import sys
import time
from array import array
import os
import glob

def main():
	
	sigXS = [11.,10.,10.,10.,10.,10.]; # in pb
	masses = [50,100,150,200,250,300];

	idir = "datacards";

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

	print "l_exp = ", l_exp
	print "l_obs = ", l_obs

	# gr_mu_exp    = MakeAGraph( masses, l_exp );
	# gr_mu_obs    = MakeAGraph( masses, l_obs );
	# gr_mu_1sigma = MakeAFillGraph( masses, l_m1sig, l_p1sig );
	# gr_mu_2sigma = MakeAFillGraph( masses, l_m2sig, l_p2sig );


	a_xax = array('d', []);
	a2_xax = array('d', []);
	a_exp = array('d', []);
	a_obs = array('d', []);
	a_1sig = array('d', []);
	a_2sig = array('d', []);

	for i in range(len(names)): a_xax.append( float(masses[i]) );
	for i in range(len(names)): a2_xax.append( float(masses[i]) );
	for i in range(len(names)-1,-1,-1): a2_xax.append( float(masses[i]) );
	for i in range(len(l_obs)): a_obs.append( float(l_obs[i]) );
	for i in range(len(l_exp)): a_exp.append( float(l_exp[i]) );
	
	for i in range(len(l_m2sig)): a_2sig.append( float(l_m2sig[i]) );
	for i in range(len(l_p2sig)-1,-1,-1): a_2sig.append( float(l_p2sig[i]) );
	
	for i in range(len(l_m1sig)): a_1sig.append( float(l_m1sig[i]) );
	for i in range(len(l_p1sig)-1,-1,-1): a_1sig.append( float(l_p1sig[i]) )	

	a_2sig.append(results[0][0])
	a2_xax.append(50)

	g_exp = ROOT.TGraph(len(a_xax), a_xax, a_exp)
	g_obs = ROOT.TGraph(len(a_xax), a_xax, a_obs)
	g_1sig = ROOT.TGraph(len(2*a_xax), a2_xax, a_1sig)
	g_2sig = ROOT.TGraph(len(2*a_xax)+1, a2_xax, a_2sig)

	print g_2sig;

	can = ROOT.TCanvas("can","can",1200,800);
	# hrl = ROOT.TH1F("hrl","hrl",6,0,6);
	hrl = can.DrawFrame(40,0,320,15);
	hrl.GetYaxis().SetTitle("#mu = #sigma_{95%CL}/#sigma_{SMS}");
	hrl.GetXaxis().SetTitle("Z\' mass (GeV)");
	# hrl.GetXaxis().SetBinLabel(1,names[0])
	# hrl.GetXaxis().SetBinLabel(2,names[1])
	# hrl.GetXaxis().SetBinLabel(3,names[2])
	# hrl.GetXaxis().SetBinLabel(4,names[3])
	# hrl.GetXaxis().SetBinLabel(5,names[4])
	# hrl.GetXaxis().SetBinLabel(6,names[5])
	hrl.SetMaximum(10.);
	hrl.Draw();

	can.SetGrid(); 

	txta = ROOT.TLatex(0.15,0.92,"CMS");
	txta.SetNDC();
	txtb = ROOT.TLatex(0.23,0.92,"Preliminary");
	txtb.SetNDC(); txtb.SetTextFont(52);
	txtc = ROOT.TLatex(0.70,0.92,"0.46 fb^{-1} (13 TeV)");
	txtc.SetNDC(); txtc.SetTextFont(42); txtc.SetTextSize(0.04);

	txtd = ROOT.TLatex(0.60,0.80,"g_{B} = 1 or g_{q} = 1/6");
	txtd.SetNDC(); txtd.SetTextFont(42); txtd.SetTextSize(0.06);

	leg = ROOT.TLegend(0.20,0.65,0.4,0.85);
	leg.SetFillStyle(1001);
	leg.SetFillColor(0);    
	leg.SetBorderSize(1);  
	# leg.SetNColumns(2);
	leg.AddEntry(g_exp,"expected","l")
	leg.AddEntry(g_obs,"observed","l")
	leg.AddEntry(g_2sig,"expected 2#sigma","f")
	leg.AddEntry(g_1sig,"expected 1#sigma","f")
   
	oneLine = ROOT.TF1("oneLine","1",0,6);
	oneLine.SetLineColor(ROOT.kRed+2);
	oneLine.SetLineWidth(2);
	oneLine.SetLineStyle(1);
	
	g_1sig.SetFillColor(ROOT.kGreen);
	g_1sig.SetFillStyle(3244);
	g_2sig.SetFillColor(ROOT.kYellow);
	g_2sig.SetFillStyle(3244);
	g_exp.SetLineStyle(2);
	g_exp.SetLineWidth(2);
	g_2sig.Draw('f');
	g_1sig.Draw('fsames');
	g_obs.Draw('lsames');
	g_exp.Draw('lsames');
	oneLine.Draw("LSAMES");
	txta.Draw();
	txtb.Draw();
	txtc.Draw();
	txtd.Draw();
	
	leg.Draw();

	can.SaveAs('brazilFullSim.pdf');

def getAsymLimits(fname,tag):

    lims = [0]*6;
    
    # if not os.path.isfile(file): 
    #     print "Warning (GetAsymLimits): "+file+" does not exist"
    #     return lims;

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