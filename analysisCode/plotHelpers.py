import ROOT
from ROOT import TFile, TTree, TChain, gPad, gDirectory, TVirtualFitter
import math
import sys
import time
import array

def makeCanvas(hists,normalize=False):

	color = [1,2,4,6,7,8,3,5]
	options = ["pe",
			   "pesames",
			   "pesames",
			   "pesames",
			   "pesames",
			   "pesames",
			   "pesames",
			   "pesames",
			   "pesames"]

	c = ROOT.TCanvas("c"+hists[0].GetName(),"c"+hists[0].GetName(),1000,800);

	max = -999;

	for i in range(len(hists)):
		hists[i].SetLineColor(color[i])
		hists[i].SetMarkerColor(color[i])
		hists[i].SetMarkerStyle(20)

		if hists[i].GetMaximum() > max: 
			max = hists[i].GetMaximum();
			hists[0].SetMaximum(max*1.25);
		if normalize and hists[i].Integral() > 0: hists[i].Scale(1./hists[i].Integral())

		hists[i].Draw(options[i]);

	c.SaveAs("plots_rhalph/"+hists[0].GetName()+".pdf")
	ROOT.gPad.SetLogy();
	c.SaveAs("plots_rhalph/"+hists[0].GetName()+"_log.pdf")

def makeCanvasDataMC(hd,hmcs,legname,name,pdir="plots",nodata=False):
	
	color = [2,4,6,7,8,3,5]
	for h in range(len(hmcs)): 
		hmcs[h].SetFillStyle(1001);
		hmcs[h].SetLineColor(1);
		hmcs[h].SetFillColor(color[h])

	hstack = ROOT.THStack("hstack","hstack");
	for h in hmcs: hstack.Add(h);
	fullmc = hstack.GetStack().Last();

	# normalize MC to data
	scalefactor = hd.Integral()/fullmc.Integral();
	for i in range(len(hmcs)): hmcs[i].Scale( scalefactor );

	xtitle = hmcs[0].GetXaxis().GetTitle();
	ytitle = hmcs[0].GetYaxis().GetTitle();
	hstack2 = ROOT.THStack("hstack2",";"+xtitle+";"+ytitle+";");
	for h in hmcs: hstack2.Add(h);

	maxval = 1.5*max(hstack2.GetStack().Last().GetMaximum(),hd.GetMaximum());
	# print maxval;
	leg = ROOT.TLegend(0.6,0.7,0.9,0.9);
	leg.SetFillStyle(0);
	leg.SetBorderSize(0);
	leg.SetTextSize(0.035);
	leg.AddEntry(hd,"data","pe");
	for i in range(len(hmcs)):
		leg.AddEntry(hmcs[i],legname[i],"f")
	# print hstack2.GetStack().Last().Integral(), hstack.GetStack().Last().Integral(),hd.Integral()
	# print hstack2.GetStack().Last().GetMaximum(),hd.GetMaximum())

	tag1 = ROOT.TLatex(0.7,0.95,"0.44 fb^{-1} (13 TeV)")
	tag1.SetNDC();
	tag1.SetTextSize(0.035);
	tag2 = ROOT.TLatex(0.17,0.95,"CMS preliminary")
	tag2.SetNDC();
	tag2.SetTextSize(0.035);

	c = ROOT.TCanvas("c"+name,"c"+name,1000,800);
	if nodata:
		hstack.SetMaximum(maxval);
		hstack.Draw("hist");
	else:
		hstack2.SetMaximum(maxval);
		hstack2.Draw("hist");
		hd.Draw("pesames");
	# ROOT.gPad.Update();
	# hstack2.GetXaxis.SetTitle( hmcs[0].GetXaxis().GetTitle() );
	# hstack2.GetYaxis.SetTitle( hmcs[0].GetYaxis().GetTitle() );	
	leg.Draw();
	tag1.Draw();
	tag2.Draw();
	c.SaveAs(pdir+"/"+name+".pdf");

	ROOT.gPad.SetLogy();
	hstack.SetMinimum(0.1);
	c.SaveAs(pdir+"/"+name+"_log.pdf")


def makeCanvasDataMC_wpred(hd,gpred,hmcs,legname,name,pdir="plots",blind=True):
	
	gpred.SetLineColor(2);
	gpred.SetFillColor(2);
	gpred.SetFillStyle(3001);

	color = [2,4,6,7,8,3,5]
	for h in range(len(hmcs)): 
		hmcs[h].SetFillStyle(0);
		hmcs[h].SetLineColor(4);
		hmcs[h].SetFillColor(0)

	hstack = ROOT.THStack("hstack","hstack");
	for h in hmcs: hstack.Add(h);
	fullmc = hstack.GetStack().Last();

	# normalize MC to data
	scalefactor = hd.Integral()/fullmc.Integral();
	for i in range(len(hmcs)): hmcs[i].Scale( scalefactor );

	xtitle = hmcs[0].GetXaxis().GetTitle();
	ytitle = hmcs[0].GetYaxis().GetTitle();
	hstack2 = ROOT.THStack("hstack2",";"+xtitle+";"+ytitle+";");
	for h in hmcs: hstack2.Add(h);

	maxval = 1.5*max(hstack2.GetStack().Last().GetMaximum(),hd.GetMaximum());
	# print maxval;
	leg = ROOT.TLegend(0.6,0.7,0.9,0.9);
	leg.SetFillStyle(0);
	leg.SetBorderSize(0);
	leg.SetTextSize(0.035);
	leg.AddEntry(hd,"data","pe");
	leg.AddEntry(gpred,"bkg pred.","f");
	for i in range(len(hmcs)):
		leg.AddEntry(hmcs[i],legname[i],"f")
	# print hstack2.GetStack().Last().Integral(), hstack.GetStack().Last().Integral(),hd.Integral()
	# print hstack2.GetStack().Last().GetMaximum(),hd.GetMaximum())

	tag1 = ROOT.TLatex(0.7,0.95,"0.44 fb^{-1} (13 TeV)")
	tag1.SetNDC();
	tag1.SetTextSize(0.035);
	tag2 = ROOT.TLatex(0.17,0.95,"CMS preliminary")
	tag2.SetNDC();
	tag2.SetTextSize(0.035);

	c = ROOT.TCanvas("c"+name,"c"+name,1000,800);
	mcall = hstack2.GetStack().Last()
	if not blind: 
		hd.SetMaximum(maxval);
		hd.Draw("pe");
		gpred.Draw("2");
		mcall.Draw("histsames");
		hd.Draw("pesames");
	if blind: 
		mcall.SetMaximum(maxval);
		mcall.Draw("hist");
		gpred.Draw("2");
		mcall.Draw("histsames");
	# ROOT.gPad.Update();
	# hstack2.GetXaxis.SetTitle( hmcs[0].GetXaxis().GetTitle() );
	# hstack2.GetYaxis.SetTitle( hmcs[0].GetYaxis().GetTitle() );	
	leg.Draw();
	tag1.Draw();
	tag2.Draw();
	c.SaveAs(pdir+"/"+name+".pdf");

	ROOT.gPad.SetLogy();
	hstack.SetMinimum(0.1);
	c.SaveAs(pdir+"/"+name+"_log.pdf")	

def makeCanvasShapeComparison(hs,legname,name,pdir="plots"):

	color = [2,4,6,7,8,3,5,2,4,6,7,8,3,5]
	style = [1,1,1,1,1,1,1,2,2,2,2,2,2,2]
	
	leg = ROOT.TLegend(0.6,0.7,0.9,0.9);
	leg.SetFillStyle(0);
	leg.SetBorderSize(0);
	leg.SetTextSize(0.035);

	maxval = -99;
	for h in range(len(hs)): 
		hs[h].SetLineColor(color[h]);
		hs[h].SetLineWidth(2);
		hs[h].SetFillStyle(0);
		hs[h].Scale(1./hs[h].Integral());
		if hs[h].GetMaximum() > maxval: maxval = hs[h].GetMaximum();
		leg.AddEntry(hs[h],legname[h],"l");

	c = ROOT.TCanvas("c"+name,"c"+name,1000,800);
	hs[0].SetMaximum(1.5*maxval);
	hs[0].Draw("hist");
	for h in range(1,len(hs)): hs[h].Draw("histsames"); 
	leg.Draw();
	c.SaveAs(pdir+"/"+name+".pdf");	
	ROOT.gPad.SetLogy();
	c.SaveAs(pdir+"/"+name+"_log.pdf")	

def	makeCanvas2D( TFMap, name, pdir='plots' ):

	c1 = ROOT.TCanvas("c1","c1",1000,800);
	TFMap.Draw("colz");
	c1.SetRightMargin(0.15);
	c1.SaveAs(pdir+"/"+name+".pdf");

	hxs = [];
	hys = [];
	
	for i in range(TFMap.GetNbinsX()):
		xnam = TFMap.GetYaxis().GetTitle();
		nbin = TFMap.GetNbinsY();
		ylo = TFMap.GetYaxis().GetBinLowEdge(1);
		yhi = TFMap.GetYaxis().GetBinUpEdge(nbin);
		hxs.append( ROOT.TH1F("hxs"+str(i),";"+xnam+";",nbin,ylo,yhi) );
	
	for i in range(TFMap.GetNbinsY()):
		xnam = TFMap.GetXaxis().GetTitle();
		nbin = TFMap.GetNbinsX();
		ylo = TFMap.GetXaxis().GetBinLowEdge(1);
		yhi = TFMap.GetXaxis().GetBinUpEdge(nbin);
		hys.append( ROOT.TH1F("hys"+str(i),";"+xnam+";",nbin,ylo,yhi) );		

	for i in range(TFMap.GetNbinsX()):
		for j in range(TFMap.GetNbinsY()):
			hxs[i].SetBinContent( j+1, TFMap.GetBinContent(i+1,j+1) );
			hxs[i].SetBinError( j+1, TFMap.GetBinError(i+1,j+1) );

	for i in range(TFMap.GetNbinsY()):
		for j in range(TFMap.GetNbinsX()):
			hys[i].SetBinContent( j+1, TFMap.GetBinContent(j+1,i+1) );			
			hys[i].SetBinError( j+1, TFMap.GetBinError(j+1,i+1) );			

	colors = [];
	for i in range(10):
		colors.append(1); colors.append(2); colors.append(4); colors.append(6);
	for i in range(len(hxs)): 
		hxs[i].SetLineColor(colors[i]);
		hxs[i].SetMarkerSize(0);
	for i in range(len(hys)): 
		hys[i].SetLineColor(colors[i]);
		hys[i].SetMarkerSize(0);

	cx = ROOT.TCanvas("cx","cx",1000,800);
	hxs[0].SetMaximum( 1.25*TFMap.GetMaximum() );
	hxs[0].SetMinimum( 0. );
	hxs[0].Draw("histe");
	for i in range(1,len(hxs)):
		hxs[i].Draw("histesames")
	cx.SaveAs(pdir+"/"+name+"_hxs.pdf");

	cy = ROOT.TCanvas("cy","cy",1000,800);
	hys[0].SetMaximum( 1.25*TFMap.GetMaximum() );
	hys[0].SetMinimum( 0. );
	hys[0].Draw("histe");
	for i in range(1,len(hys)):
		hys[i].Draw("histesames")
	cy.SaveAs(pdir+"/"+name+"_hys.pdf");

	return hys;

def dummy():
	print "hi";

