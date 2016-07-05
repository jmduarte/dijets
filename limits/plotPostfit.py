#! /usr/bin/env python
import ROOT as r,sys,math,array,os
from optparse import OptionParser
from ROOT import std,RooDataHist
import ROOT

import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.08);
ROOT.gStyle.SetPadLeftMargin(0.14);
ROOT.gStyle.SetPadTopMargin(0.10);

ROOT.gStyle.SetPalette(1);
parser = OptionParser()
parser.add_option('--input'  ,action='store',type='string',dest='input'    ,default='mlfit.root',help='input file')
parser.add_option('--data'   ,action='store',type='string',dest='data'     ,default='combine_Zprime150_0.root',help='data file')
parser.add_option('--name'   ,action='store',type='string',dest='name'     ,default='outname',help='data file')
parser.add_option('--cats'    ,action='store',type='string',dest='cats'    ,default='1',help='categories 0,1,2,...')
parser.add_option('--passfail',action='store_true',         dest='passfail',default=False,help='fit pass and failng') 
parser.add_option('--divide'  ,action='store_true',         dest='divide'  ,default=False,help='ratio') 
parser.add_option('--mass'    ,action='store',              dest='mass'    ,default=100  ,help='mass') 

(options,args) = parser.parse_args()


def main():

	lHFile = r.TFile(options.input)
	lDFile = r.TFile(options.data)
	lDSum=[]
	lSum=[]
	for cat in options.cats.split(','):
		lData  = loadData(lDFile,"pass_cat"+cat)
		lHists = loadHist(lHFile,"ch0",options.mass)
		if options.passfail:
			lData .extend(loadData(lDFile,"fail_cat"+cat))
			lHists.extend(loadHist(lHFile,"ch"+cat+"_fail_cat"+cat,options.mass,True,True))
		if len(lSum) == 0:
			lDSum = lData
			lSum  = lHists
		else:
			for i0 in range(0,len(lDSum)):
				lDSum[i0].Add(lData[i0])
			for i0 in range(0,len(lSum)):
				lSum[i0].Add(lHists[i0])
		draw(lData,lHists,options.name,options.passfail,options.divide)
	# draw(lDSum,lSum,"sum",options.passfail,options.divide)

# def end():
# 	#return
# 	if __name__ == '__main__':
# 		rep = ''
# 		while not rep in [ 'q', 'Q','a',' ' ]:
# 			rep = raw_input( 'enter "q" to quit: ' )
# 			if 1 < len(rep):
# 				rep = rep[0]

def divide(iData,iHists):
	lPass = -1
	lFail = -1
	for i0 in range(0,len(iHists)):
		iHists[i0].Sumw2()
		if iHists[i0].GetName().find("") > -1 and lPass == -1:
			lPass = i0
		if iHists[i0].GetName().find("fail") > -1 and lFail == -1:
			lFail = i0
	if lPass != -1:
		iData[0].Divide(iHists[lPass])
		for i0 in range(0,iData[0].GetNbinsX()+1):
			if iHists[lPass].GetBinContent(i0) > 0:
				iData[0].SetBinError(i0,iData[0].GetBinError(i0)/iHists[lPass].GetBinContent(i0))
		for i0 in range(0,len(iHists)):
			if i0 != lPass and iHists[i0].GetName().find("") > -1:
				iHists[i0].Divide(iHists[lPass])
		iHists[lPass].Divide(iHists[lPass])

	if lFail != -1:
		iData[1].Divide(iHists[lFail])
		for i0 in range(0,len(iHists)):
			if i0 != lFail and iHists[i0].GetName().find("fail") > -1:
				iHists[i0].Divide(iHists[lFail])
		iHists[lFail].Divide(iHists[lFail])

def draw(iData0,iHists0,iName="A",iPF=True,iDivide=False):

	iData = [];
	for i,h in enumerate(iData0): iData.append( reRangeHistogram(h,30,330) );
	iHists = [];
	for i,h in enumerate(iHists0): 
		iHists.append( reRangeHistogram(h,30,330) );
		iHists[i].SetLineWidth(2);

	txta = ROOT.TLatex(0.15,0.92,"CMS");
	txta.SetNDC();
	txtb = ROOT.TLatex(0.22,0.92,"Preliminary");
	txtb.SetNDC(); txtb.SetTextFont(52);
	txtc = ROOT.TLatex(0.73,0.92,"2.7 fb^{-1} (13 TeV)");
	txtc.SetNDC(); txtc.SetTextFont(42); txtc.SetTextSize(0.045);
	txtd = ROOT.TLatex(0.60,0.80,"g_{B} = 1 or g_{q} = 1/6");
	txtd.SetNDC(); txtd.SetTextFont(42); txtd.SetTextSize(0.04);

	if iDivide:
		divide(iData,iHists);

	lC0 = r.TCanvas("Can"+iName,"Can"+iName,900,800);
	# if iPF:
	# 	lC0.Divide(2)
	# lC0.cd(1)

	p12 = ROOT.TPad("p12","p12",0.0,0.3,1.0,1.0);
	p22 = ROOT.TPad("p22","p22",0.0,0.0,1.0,0.3);
	p12.SetBottomMargin(0.05);
	p22.SetTopMargin(0.05);
	p22.SetBottomMargin(0.3);

	lC0.cd();
	p12.Draw(); p12.cd();

	lHist=fix(iHists[0],iData[0])
	
	iHists[0].SetMarkerSize(0);
	iHists[0].SetLineWidth(0);
	iHists[0].SetFillColor(4);
	lHist.GetYaxis().SetTitleOffset(0.9);
	lHist.GetYaxis().SetTitle("N events");
	lHist.GetYaxis().SetTitleSize(0.06);
	iRatio = iData[0].Clone();
	iRatio.Divide(iHists[2]);
	iRatio.SetTitle("; soft drop mass (GeV); Data/Prediction");
	iRatio.GetYaxis().SetTitleSize(0.12);
	iRatio.GetYaxis().SetLabelSize(0.1);
	iRatio.GetYaxis().SetTitleOffset(0.4);
	iRatio.GetXaxis().SetTitleSize(0.12);
	iRatio.GetXaxis().SetLabelSize(0.1);
	iRatio.GetXaxis().SetTitleOffset(0.9);
	iRatio.GetYaxis().SetRangeUser(0.5,1.5);
	iOneWithErrors = iHists[2].Clone();
	iOneWithErrors.Divide(iHists[2].Clone());
	for i in range(iOneWithErrors.GetNbinsX()): iOneWithErrors.SetBinError( i+1, iHists[2].GetBinError(i+1)/iHists[2].GetBinContent(i+1) );
	iOneWithErrors.SetFillStyle(3001);
	iOneWithErrors.SetFillColor(4);
	iOneWithErrors.SetMarkerSize(0);
	iOneWithErrors.SetLineWidth(0);

	iW = iHists[1].Clone();
	iZ = iHists[2].Clone();
	iW.Add(iHists[0],-1);
	iZ.Add(iHists[1],-1);
	iZ.SetLineColor(ROOT.kGreen+2)

	iHists[2].Chi2Test(lHist,"P")

	lHist.Draw("ep")
	# pDraw=False
	# for pHist in iHists:
	# 	pHist.Draw("hist sames") if pDraw else pHist.Draw("e2 sames")
	# 	pDraw=True
	iHists[0].Draw("e2 sames");
	# iHists[1].Draw("hist sames");
	iHists[2].Draw("hist sames");
	iW.Draw("hist sames");
	iZ.Draw("hist sames");
	iHists[4].SetLineColor(6);
	iHists[4].SetLineStyle(2);
	iHists[4].Draw("hist sames");

	lHist.Draw("ep sames")
	lLegend = r.TLegend(0.63,0.53,0.88,0.84)
	lLegend.SetFillColor(0)
	lLegend.SetBorderSize(0)
	lLegend.AddEntry(iData [0],"data","lp")
	lLegend.AddEntry(iHists[0],"QCD Pred.","f")
	# lLegend.AddEntry(iHists[1],"W#rightarrow qq","l")
	lLegend.AddEntry(iHists[2],"Total SM Pred.","l")
	lLegend.AddEntry(iW,"W(qq)","l")
	lLegend.AddEntry(iZ,"Z(qq)","l")
	lLegend.AddEntry(iHists[4],"Z'(qq), g_{B} = 1","l")
	#lLegend.AddEntry(iHists[3],"Signal","lf")
	lLegend.Draw()
	txta.Draw();
	txtb.Draw();
	txtc.Draw();

	lC0.cd();
	p22.Draw(); p22.cd();	
	p22.SetGrid();
	iRatio.Draw();
	iOneWithErrors.Draw("e2 sames");
	iRatio.Draw("sames");

	if not iPF:
		lC0.SaveAs(iName+".png")
		lC0.SaveAs(iName+".pdf")
		# end()
		return

	# lC0.cd(2)
	# iData[1].Draw("ep")
	# pDraw=False
	# for pHist in iHists:
	# 	if pHist.GetName().find("fail") > -1:
	# 		pHist.Draw("hist sames") if pDraw else pHist.Draw("e2 sames")
	# 		pDraw=True
	# lC0.Modified()
	# lC0.Update()
	# iData[1].Draw("ep sames")
	# lC0.SaveAs(iName+".png")
	# lC0.SaveAs(iName+".pdf")
	# end()

def norm(iFile,iH,iName) :
	lNorms = iName.split("/")[0].replace("shapes","norm")
	print "---",lNorms
	lArg = iFile.Get(lNorms)
	lName = iName.replace(iName.split("/")[0]+"/","")
	lVal = r.RooRealVar(lArg.find(lName)).getVal();
	print "!!!",lName,lNorms,iH.Integral(),lVal
	iH.Scale(lVal/iH.Integral())

def load(iFile,iName):
	lHist = iFile.Get(iName)
	norm(iFile,lHist,iName)
	lHist.SetName(iName.replace("/","_"))
	lHist.SetTitle(iName.replace("/","_"))
	return lHist

def fix(iH0,iH1):
	lTHist = r.TH1F(iH1.GetName(),iH1.GetName(),iH0.GetNbinsX(),iH0.GetXaxis().GetXmin(),iH0.GetXaxis().GetXmax())
	for i0 in range(0,iH0.GetNbinsX()+1):
		lTHist.SetBinContent(i0,iH1.GetBinContent(i0))
		lTHist.SetBinError(i0,iH1.GetBinError(i0))
	lTHist.SetFillStyle(iH1.GetFillStyle())
	lTHist.SetFillColor(iH1.GetFillColor())
	lTHist.SetLineColor(iH1.GetLineColor())
	lTHist.SetLineStyle(iH1.GetLineStyle())
	return lTHist
	
def loadHist(iFile,iCat,iMass,iEWK=True,iS=True):
	lHists = []
	lFit = "shapes_fit_s/"+iCat+"/" if iS else "shapes_fit_b/"+iCat+"/"
	#lFit = "shapes_prefit/"+iCat+"/" if iS else "shapes_fit_b/"+iCat+"/"
	lHists.append(load(iFile,lFit+"qcd"))
	if iEWK:
		lId = len(lHists)
		lHists.append(load(iFile,lFit+"Winc"))
		lHists.append(load(iFile,lFit+"Zinc"))
		if iS:
			lHists.append(load(iFile,lFit+"sig"))
		lHists[lId].SetLineColor(46)
		lHists[lId+1].SetLineColor(9)
		if iS:
			lHists[lId+2].SetLineColor(8)
		lHists[lId].SetLineWidth(2)
		lHists[lId+1].SetLineWidth(2)
		if iS :
			lHists[lId+2].SetLineWidth(2)
		lHists[lId].Add(lHists[0])
		lHists[lId+1].Add(lHists[1])
		if iS :
			lHists[lId+2].Add(lHists[2])
	lHists[0].SetFillColor(16)
	lHists[0].SetFillStyle(3001)
	lHists[0].SetLineStyle(2)

	lHists.append(load(iFile,"shapes_prefit/ch0/sig"))	
	return lHists

def loadData(iDataFile,iCat):
	lData = iDataFile.Get("data_obs")
	lData.GetXaxis().SetTitle("m_{J} (GeV)")
	lData.SetMarkerStyle(20)
	lData.Sumw2()
	return [lData]

def reRangeHistogram(ih,xlo,xhi):

	nbins = ih.GetNbinsX();
	hnew = r.TH1F("hnew_"+ih.GetName(),";soft drop mass (GeV); N",nbins,xlo,xhi);
	for i in range(nbins):
		hnew.SetBinContent( i+1, ih.GetBinContent(i+1) );
		hnew.SetBinError( i+1, ih.GetBinError(i+1) );
		hnew.SetLineColor( ih.GetLineColor() );
		hnew.SetLineStyle( ih.GetLineStyle() );
		hnew.SetFillColor( ih.GetFillColor() );
		hnew.SetFillStyle( ih.GetFillStyle() );

	return hnew

if __name__ == "__main__":

	main();

