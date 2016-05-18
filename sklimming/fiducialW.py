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
	
	f = ROOT.TFile("zprimebits_17May/W.root")
	t = f.Get("Events");

	h_all = ROOT.TH1F("h_all",";pt;N",100,0,2000)
	h_all_nok = ROOT.TH1F("h_all_nok",";pt;N",100,0,2000)	
	h_cut_nok = ROOT.TH1F("h_cut_nok",";pt;N",100,0,2000)	
	h_cut = ROOT.TH1F("h_cut",";pt;N",100,0,2000)

	for i in range(t.GetEntries()):
		t.GetEntry(i);
		h_all.Fill( t.genVPt, t.kfactor )
		h_all_nok.Fill( t.genVPt );
		if t.genVPt > 500.: 
			h_cut.Fill( t.genVPt, t.kfactor )
			h_cut_nok.Fill( t.genVPt )

	
	f2 = ROOT.TFile("zprimebits_17May/ZPrimeToQQ_200GeV_v4_mc.root")
	t2 = f2.Get("Events");

	hsig_all = ROOT.TH1F("h_all",";pt;N",100,0,2000)
	hsig_all_nok = ROOT.TH1F("h_all_nok",";pt;N",100,0,2000)	
	hsig_cut_nok = ROOT.TH1F("h_cut_nok",";pt;N",100,0,2000)	
	hsig_cut = ROOT.TH1F("h_cut",";pt;N",100,0,2000)

	for i in range(t2.GetEntries()):
		t2.GetEntry(i);
		hsig_all.Fill( t2.genVPt, t2.kfactor )
		hsig_all_nok.Fill( t2.genVPt );
		if t2.genVPt > 500.: 
			hsig_cut.Fill( t2.genVPt, t2.kfactor )
			hsig_cut_nok.Fill( t2.genVPt )


	print h_all_nok.Integral();
	print h_all.Integral();
	print h_cut_nok.Integral();	
	print h_cut.Integral();

	print hsig_all_nok.Integral();
	print hsig_all.Integral();
	print hsig_cut_nok.Integral();	
	print hsig_cut.Integral();

	print "fraction of LO XS with pT > 500: ", h_cut_nok.Integral()/h_all_nok.Integral()
	print "NLO XS/LO XS with pT > 500: ", h_cut.Integral()/h_cut_nok.Integral()

	c = ROOT.TCanvas("c","c",1000,800);
	h_all_nok.SetLineColor(1);
	h_all.Draw("hist");
	h_cut.SetLineColor(2);
	h_cut.Draw("histsames");
	h_all_nok.SetLineColor(4);
	h_all_nok.Draw("histsames");
	h_cut_nok.SetLineColor(6);
	h_cut_nok.Draw("histsames");
	ROOT.gPad.SetLogy();
	c.SaveAs("WpT.pdf");

if __name__ == '__main__':

	main();