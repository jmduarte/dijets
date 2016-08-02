import ROOT
from optparse import OptionParser
from operator import add
from fullLims_1cat import getAsymLimits,makeAFillGraph,makeAGraph
import math
import sys
import time
from array import array
import os
import glob
import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.18);
ROOT.gStyle.SetPadLeftMargin(0.18);
ROOT.gStyle.SetPadTopMargin(0.10);
ROOT.gStyle.SetPalette(1);


def end():
    #if __name__ == '__main__':
        rep = ''
        while not rep in [ 'q', 'Q','a',' ' ]:
            rep = raw_input( 'enter "q" to quit: ' )
            if 1 < len(rep):
                rep = rep[0]

def load():
    idir = "datacards"
    #masses = [50,60,75,90,100,110, 125, 135, 150, 165, 180, 200, 250, 300]
    masses = [100,110, 125, 135, 150, 165, 180, 200, 250, 300]
    KFACTOR = 1.218
    idir = "datacards"
    results = [];
    for i in range(len(masses)):
        results.append( getAsymLimits(idir+'/higgsCombineZprime'+str(masses[i])+'.Asymptotic.mH120.root','Zp'+str(masses[i])) );

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
		l_m2sig.append(math.sqrt(r[1]));
		l_m1sig.append(math.sqrt(r[2]));
		l_exp.append  (math.sqrt(r[3]));
		l_p1sig.append(math.sqrt(r[4]));
		l_p2sig.append(math.sqrt(r[5]));
		l_obs.append  (math.sqrt(r[0]));
		ctr+=1;

    gr_mu_exp    = makeAGraph( masses, l_exp, 1, 2 );
    gr_mu_obs    = makeAGraph( masses, l_obs, 1, 1 );
    gr_mu_1sigma = makeAFillGraph( masses, l_m1sig, l_p1sig, 0, 3, 1001 );
    gr_mu_2sigma = makeAFillGraph( masses, l_m2sig, l_p2sig, 0, 5, 1001 );    
    return [gr_mu_exp,gr_mu_obs]

def avwidth(iType,g,med,mfm):
    front=g*g*(med*med+2*mfm*mfm)/(12*med*3.141159265)
    if abs(iType) == 2:
        front=g*g*(med*med-4*mfm*mfm)/(12*med*3.141159265)
    if 2.*mfm > med:
        return 0.001
    sqrtV=math.sqrt(1-4*(mfm/med)*(mfm/med))
    return front*sqrtV

def avtotwidth(iType,gdm,gsm,med,mdm):
    u=avwidth(iType,gsm,med,0.001)
    d=u
    s=avwidth(iType,gsm,med,0.135)
    c=avwidth(iType,gsm,med,1.5)
    b=avwidth(iType,gsm,med,5.1)
    t=0
    if med > 2.*172.5:
        t=avwidth(iType,gsm,med,172.5)
    quarks=3*(u+d+s+c+b+t)
    dm=avwidth(iType,gdm,med,mdm)
    #print u,d,s,c,b,t,dm,quarks
    return dm+quarks

def BRCorrGQ(iGQ,iGDM,iMed,iMDM):
    lNewWidth = avtotwidth(1,iGDM,iGQ,iMed,iMDM)
    lOldWidth = avtotwidth(1,0.  ,iGQ,iMed,1.)
    lDelta=lNewWidth-lOldWidth
    lCorr = 0.5+0.5*math.sqrt(1+lDelta/lOldWidth)
    lGQNew = math.sqrt(lCorr)*iGQ
    #if iMDM < 10 and iMed == 200:
    #    print "Scale:",iGQ,lGQNew,lDelta/lOldWidth,lCorr,iMed,iMDM
    return lGQNew

def convertgq(iGraph,iGDM,iMDM,iConvertGB=True):
    xnew = array('d', []);
    ynew = array('d', []);
    for i0 in range(0,iGraph.GetN()):
        xnew.append(iGraph.GetX()[i0])
        ynew.append(BRCorrGQ(iGraph.GetY()[i0]/6.,iGDM,iGraph.GetX()[i0],iMDM))
        if iConvertGB:
            ynew[i0]=ynew[i0]
    lGraph    = makeAGraph( xnew, ynew, 1, 3);
    return lGraph

def getCont(iGraph,iVal,color):
  #iGraph.SetContour(1,iVal)
  #iGraph.Draw("contz list")
  lContours = iGraph.GetContourList(iVal);
  print lContours,lContours.GetSize()
  xnew = array('d', []);
  ynew = array('d', []);
  for i0 in range(0,lContours.GetSize()):
    pCurv = lContours.At(i0)
    #pCurv = pCont.First()
    #for i1 in range(0,pCont.GetSize()):
    for i2 in range(0,pCurv.GetN()):
        xnew.append(pCurv.GetX()[i2])
        ynew.append(pCurv.GetY()[i2])
      #pCurv = pContLevel.fter(pCurv)
  lGraph    = makeAGraph( xnew, ynew, color, 3 )
  lGraph.GetXaxis().SetTitle("m_{med} (GeV)")
  lGraph.GetYaxis().SetTitle("m_{dm}  (GeV)")
  lGraph.SetLineWidth(4)
  return lGraph

def make2DGraph(gr,gdm,canv,leg,label,color):
    mdmarr = array('d', []);
    medarr = array('d', []);
    gqarr  = array('d', []);
    for i0 in range(0,100):
        gqgraph=convertgq(gr,gdm,i0*20,True)
        for i1 in range(0,gqgraph.GetN()):
            mdmarr .append(i0*20)
            medarr.append(gqgraph.GetX()[i1])
            gqarr .append(gqgraph.GetY()[i1])
    lOGraph = ROOT.TGraph2D(len(mdmarr),medarr,mdmarr,gqarr)
    lOGraph.GetXaxis().SetTitle("m_{med} (GeV)")
    lOGraph.GetYaxis().SetTitle("m_{dm}  (GeV)")
    lOGraph.GetZaxis().SetTitle("g_{q}")
    lOGraph.Draw("colz")
    canv.Update()
    lCont = getCont(lOGraph,0.25,color)
    leg.AddEntry(lCont,label,"l")
    return lOGraph,lCont

def main():
    grexp,grobs=load()
    gdm=1
    leg   = ROOT.TLegend(0.20,0.55,0.4,0.85)
    canv0 = ROOT.TCanvas("can0","can0",1200,800)
    leg.SetHeader("g_{DM}="+str(gdm))
    leg.SetFillColor(0)    
    leg.SetBorderSize(0)  
    lExp,lXExp=make2DGraph(grexp,gdm,canv0,leg,"expected 95% CL.",1)
    lObs,lXObs=make2DGraph(grobs,gdm,canv0,leg,"observed 95% CL.",ROOT.kOrange)
    lObs.Draw("colz")
    lXExp.Draw("l sames")
    lXObs.Draw("l sames")
    leg.Draw()
    lObs.GetXaxis().SetTitle("m_{med} (GeV)")
    lObs.GetYaxis().SetTitle("m_{dm}  (GeV)")
    lObs.GetZaxis().SetTitle("g_{q}")
    ROOT.gPad.Modified()
    ROOT.gPad.RedrawAxis()
    canv0.SaveAs("gq_mdm_mmed.png")
    canv0.SaveAs("gq_mdm_mmed.pdf")
    lFile = ROOT.TFile("MMedMDM.root","RECREATE")
    lObs.SetName("obs")
    lObs.SetTitle("obs")
    lExp.SetName("exp")
    lExp.SetTitle("exp")

    lXObs.SetName("obs_025")
    lXObs.SetTitle("obs_025")
    lXExp.SetName("exp_025")
    lXExp.SetTitle("exp_025")
    lObs.Write()
    lExp.Write()
    lXObs.Write()
    lXExp.Write()
    end()
        
if __name__ == '__main__':

	main();
