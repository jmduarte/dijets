#! /usr/bin/env python
import ROOT as r,sys,math,array,os
from optparse import OptionParser
from ROOT import std,RooDataHist
r.gROOT.Macro(os.path.expanduser('~/rootlogon.C'))
r.gInterpreter.GenerateDictionary("std::pair<std::string, RooDataHist*>", "map;string;RooDataHist.h")
r.gInterpreter.GenerateDictionary("std::map<std::string, RooDataHist*>", "map;string;RooDataHist.h")
#import rootpy.stl as stl
#from rootpy.plotting import Hist
#MapStrRootPtr = stl.map (stl.string, "TH1*")
#StrHist       = stl.pair(stl.string, "TH1*" )

fVars=[]
fDatas=[]
fFuncs=[]
fPars=[]

def end():
    if __name__ == '__main__':
        rep = ''
        while not rep in [ 'q', 'Q','a',' ' ]:
            rep = raw_input( 'enter "q" to quit: ' )
            if 1 < len(rep):
                rep = rep[0]

def parser():
    parser = OptionParser()
    parser.add_option('--card'    ,action='store_true',         dest='card'    ,default=False,       help='categorized bin by bin card fit')
    parser.add_option('--cat'     ,action='store_true',         dest='cat'     ,default=False,       help='categorized fit')
    parser.add_option('--1D'      ,action='store_true',         dest='fit1D'   ,default=False,       help='1D fit')
    parser.add_option('--2D'      ,action='store_true',         dest='fit2D'   ,default=False,       help='2D fit')
    parser.add_option('--fail'    ,action='store_true',         dest='fail'    ,default=False,       help='fit failing') 
    parser.add_option('--passfail',action='store_true',         dest='passfail',default=False,       help='fit pass and failng') 
    parser.add_option('--input',action='store',type='string',dest='input',default='hists.root',help='input file')
    (options,args) = parser.parse_args()
    return options

def drawFrame(iFrame,iData,iFuncs):
    iData.plotOn(iFrame)
    iColor=50
    for pFunc in iFuncs:
        #pFunc.plotOn(iFrame,r.RooFit.LineColor(iColor),r.RooFit.LineStyle(iColor != 50+1),r.RooFit.Components(pFunc.GetName()),r.RooFit.ProjWData(iData))
        pFunc.plotOn(iFrame,r.RooFit.LineColor(iColor),r.RooFit.LineStyle(iColor != 50+1),r.RooFit.ProjWData(iData))
        iColor+=10

def draw(iVar,iData,iFuncs,iLabel="A"):
    lCan   = r.TCanvas(str(iLabel),str(iLabel),800,600)
    lFrame = iVar.frame()
    drawFrame(lFrame,iData,iFuncs)
    lFrame.Draw()
    lCan.Modified()
    lCan.Update()
    end()

def drawPF(iVar,iData,iFuncs,iLabel="A"):
    lCan   = r.TCanvas(str(iLabel),str(iLabel),800,600)
    lCan.Divide(2)
    lPFrame = iVar.frame(r.RooFit.Bins(30),r.RooFit.Title("pass"))
    lFFrame = iVar.frame(r.RooFit.Bins(30),r.RooFit.Title("fail"))
    drawFrame(lPFrame,iData[0],iFuncs[0])
    drawFrame(lFFrame,iData[1],iFuncs[1])
    lCan.cd(1)
    lPFrame.Draw()
    lCan.cd(2)
    lFFrame.Draw()
    lCan.Modified()
    lCan.Update()
    end()

def shift(iVar,iDataHist,iShift=5.):
    lDM     = r.RooRealVar   ("Xdm","Xdm", 0.,-10,10)
    lShift  = r.RooFormulaVar("Xshift",iVar.GetName()+"-Xdm",r.RooArgList(iVar,lDM))  
    lSPdf   = r.RooHistPdf(iDataHist.GetName()+"P",iDataHist.GetName()+"P", r.RooArgList(lShift),r.RooArgList(iVar),iDataHist,0)
    lDM.setVal(iShift)
    lHUp   = lSPdf.createHistogram("x")
    lUp    = r.RooDataHist(iDataHist.GetName()+"_scaleUp",iDataHist.GetName()+"_scaleUp", r.RooArgList(iVar),lHUp)
    lDM.setVal(-iShift)
    lHDown = lSPdf.createHistogram("x")
    lDown  = r.RooDataHist(iDataHist.GetName()+"_scaleDown",iDataHist.GetName()+"_scaleDown", r.RooArgList(iVar),lHDown)
    return (lUp,lDown)
    
def workspace(iDatas,iFuncs,iVars,iCat="_cat0",iShift=False):
    lW = r.RooWorkspace("w_"+str(iCat))
    for var in iVars:
        try:
            var.setConstant(True)
        except:
            pass
    for pData in iDatas:
        getattr(lW,'import')(pData,r.RooFit.RecycleConflictNodes())
    
    for pFunc in iFuncs:
        getattr(lW,'import')(pFunc,r.RooFit.RecycleConflictNodes())
        if iShift:
            (pFUp,pFDown) = shift(iVars[0],pFunc,5.)
            getattr(lW,'import')(pFUp,  r.RooFit.RecycleConflictNodes())
            getattr(lW,'import')(pFDown,r.RooFit.RecycleConflictNodes())
       
    if iCat.find("pass_cat1") == -1:
        lW.writeToFile("simple-shapes-RooDataHist.root",False)
    else:
        lW.writeToFile("simple-shapes-RooDataHist.root")

def writeHist(iH,iShift,iTemps,iBlanks=[]):  
    lOut = r.TFile("simple-shapes-Hist.root","RECREATE")
    for pHist in iH:
        if pHist.GetName().find("data") > -1 and pHist.GetName().find("data_obs") < 0:
            pHist.SetName (pHist.GetName ().replace("data","data_obs"))
            pHist.SetTitle(pHist.GetTitle().replace("data","data_obs"))
        pHist.Write()
    iShift.setVal(5.)
    for pTemp in iTemps:
      lUp = pTemp.createHistogram("x")
      lUp.SetTitle(pTemp.GetName().replace("e_","_")+"_scaleUp")
      lUp.SetName (pTemp.GetName().replace("e_","_")+"_scaleUp")
      lUp.Write()
  
    iShift.setVal(-5.)
    for pTemp in iTemps:
      lDown = pTemp.createHistogram("x")
      lDown.SetTitle(pTemp.GetName().replace("e_","_")+"_scaleDown")
      lDown.SetName (pTemp.GetName().replace("e_","_")+"_scaleDown")
      lDown.Write()

    return
    for pBlank in iBlanks:
        for i0 in range(1,pBlank.GetNbinsX()+1):
            pLabel="pass" if pBlank.GetName().find("pass") > -1 else "fail"
            pOutBlankUp  =pBlank.Clone(pBlank.GetName()+"_"+pLabel+"Bin"+str(i0)+"Up")
            pOutBlankDown=pBlank.Clone(pBlank.GetName()+"_"+pLabel+"Bin"+str(i0)+"Down")
            for i1 in range(0,pOutBlankUp.GetNbinsX()+1):
                if i0 == i1:
                    pOutBlankUp.SetBinContent(i1,0)
                    pOutBlankDown.SetBinContent(i1,2)
            pOutBlankUp.Write()
            pOutBlankDown.Write()

def baseVars(i1D=False):
    lMSD    = r.RooRealVar("x","x",40,200)
    lPt  = r.RooRealVar   ("pt","pt",500,1500)
    lRho = r.RooFormulaVar("rho","log(x*x/pt)",r.RooArgList(lMSD,lPt))
    lEff    = r.RooRealVar("veff"      ,"veff"      ,0.5 ,0.,1.0)
    lEffQCD = r.RooRealVar("qcdeff"    ,"qcdeff"   ,0.03,0.,0.2)
    lDM     = r.RooRealVar("dm","dm", 0.,-10,10)
    lShift  = r.RooFormulaVar("shift",lMSD.GetName()+"-dm",r.RooArgList(lMSD,lDM))  
    lVars=[lMSD,lEff,lEffQCD,lDM,lShift]
    if not i1D:
        lVars.extend([lPt,lRho])
    fVars.extend([lMSD,lEff,lEffQCD,lDM])
    lPt .setBins(5)
    fPars.extend([lEffQCD,lDM,lEff])
    return lVars

def bestFit(iEffQCD,iDM,iVEff,iA0,iA1,iA2,iP1,iR1,iSigma):
    iA0.setVal(0.486)
    iA1.setVal(0.00326)
    iA2.setVal(-0.0000237)
    iP1.setVal(-0.00169)
    iR1.setVal(0.00939)
    iSigma.setVal(51.2)
    iEffQCD.setVal(0.0242)
    iDM.setVal(-2.3)
    iVEff.setVal(0.037)
    iA0.setConstant(True)
    iA1.setConstant(True)
    iA2.setConstant(True)
    iP1.setConstant(True)
    iSigma.setConstant(True)
    iEffQCD.setConstant(True)
    iDM.setConstant(True)
    iVEff.setConstant(True)
        
def command(iBase,iPt,iM,iRho,iPass=True,iNR=1,iNP=1):
    lFunc=iBase
    if iPt > 0 and iNP > 0:
        lFunc+=("*(1+p1*"+str(iPt-500))
        if iNP > 1:
            lFunc+=("+p2*"+str(iPt-500)+"*"+str(iPt-500))
        lFunc+=")"
    if iNR > 0:
        lRho = '(log(x*x/'+str(iPt)+')-2.5)' if iM < 0 else str(iRho-2.5)
        lFunc+=("*(1+r1*"+lRho)
        if iNR > 1:
            lFunc+=("+r2*"+lRho+"*"+lRho)            
        lFunc+=")"
    if not iPass:
        lFunc=lFunc.replace("+p","-p")
        lFunc=lFunc.replace("+r","-r")
    return lFunc

def qcdFunc(iH,iVars,iBin="_cat0",i1D=True,iPt=-1):
    lNTot   = r.RooRealVar   ("qcdnorm"+iBin,"qcdnorm"+iBin,(iH[0].Integral()+iH[1].Integral()),0.,3.*(iH[0].Integral()+iH[1].Integral()))
    lNPass  = r.RooFormulaVar("fqpass"+iBin ,"qcdnorm"+iBin+"*(qcdeff)"  ,r.RooArgList(lNTot,iVars[2]))
    lNFail  = r.RooFormulaVar("fqfail"+iBin ,"qcdnorm"+iBin+"*(1-qcdeff)",r.RooArgList(lNTot,iVars[2]))
    lA0     = r.RooRealVar   ("a0"          ,"a0"          ,0.486,0.,1)          
    lA1     = r.RooRealVar   ("a1"          ,"a1"          ,0.00325,-0.01,0.01)
    lA2     = r.RooRealVar   ("a2"          ,"a2"          ,-0.000025,-0.0001,0.0001)
    lSigma1 = r.RooRealVar   ("sigma1"      ,"sigma1"      ,50,0,100)
    lP1     = r.RooRealVar   ("p1" ,"p1", -0.001   ,-0.1  ,0.1)
    lP2     = r.RooRealVar   ("p2" ,"p2", 0.0   ,-0.001,0.001)
    lR1     = r.RooRealVar   ("r1" ,"r1", 0.0   ,-0.5 ,0.5)
    lR2     = r.RooRealVar   ("r2" ,"r2", 0.0   ,-0.1 ,0.1)
    if i1D:
        lFunc="(exp(-x*x/sigma1/sigma1)+a0+a1*x+a2*x*x)"
        if iPt > 0:
            iVars[5].setVal(iPt)
            lPFunc = command(lFunc,iPt,-1,-1,False,2,2)
            lFFunc = command(lFunc,iPt,-1,-1,True,2,2)
            lQFuncP = r.RooGenericPdf("qcdpass"+iBin,lPFunc,r.RooArgList(iVars[0],lSigma1,lA0,lA1,lA2,lP1,lR1,lP2,lR2))
            lQFuncF = r.RooGenericPdf("qcdfail"+iBin,lFFunc,r.RooArgList(iVars[0],lSigma1,lA0,lA1,lA2,lP1,lR1,lP2,lR2))
        else:
            lQFuncP = r.RooGenericPdf("qcdpass"+iBin,lPFunc,r.RooArgList(iVars[0],lSigma1,lA0,lA1,lA2))
            lQFuncF = r.RooGenericPdf("qcdfail"+iBin,lFFunc,r.RooArgList(iVars[0],lSigma1,lA0,lA1,lA2))
    else:
        lQFuncP = r.RooGenericPdf("qcdpass"+iBin,"(exp(-x*x/sigma1/sigma1)+a0+a1*x+a2*x*x)*(1+p1*(pt-500)+r1*(rho-2.5))",r.RooArgList(iVars[0],lSigma1,lA0,lA1,lA2,lP1,iVars[5],lR1,iVars[6]))
        lQFuncF = r.RooGenericPdf("qcdfail"+iBin,"(exp(-x*x/sigma1/sigma1)+a0+a1*x+a2*x*x)*(1-p1*(pt-500)-r1*(rho-2.5))",r.RooArgList(iVars[0],lSigma1,lA0,lA1,lA2,lP1,iVars[5],lR1,iVars[6]))
    lQCDP   = r.RooExtendPdf ("qcdpassE"+iBin,"qcdpass"+iBin,lQFuncP,lNPass)
    lQCDF   = r.RooExtendPdf ("qcdfailE"+iBin,"qcdfail"+iBin,lQFuncF,lNFail)
    lQCD    = [lQCDP,lQCDF,lQFuncP,lQFuncF]
    fVars.extend([lNTot,lA0,lA1,lA2,lR1,lSigma1,lNPass,lNFail,lP1,lR1,lR2,lP2])
    fFuncs.extend(lQCD)
    fPars.extend([lA0,lA1,lA2,lP1,lR1,lP2,lR2,lSigma1])
    return lQCD

def histFunc(iH,iVars,iLabel="w",iBin="_cat0",i1D=True):
    lNTot   = r.RooRealVar (iLabel+"norm"+iBin,iLabel+"norm"+iBin,(iH[0].Integral()+iH[1].Integral()),0.,5.*(iH[0].Integral()+iH[1].Integral()))
    lNPass  = r.RooFormulaVar(iLabel+"fpass"+iBin ,iLabel+"norm"+iBin+"*(veff)"  ,r.RooArgList(lNTot,iVars[1]))
    lNFail  = r.RooFormulaVar(iLabel+"fqail"+iBin ,iLabel+"norm"+iBin+"*(1-veff)",r.RooArgList(lNTot,iVars[1]))
    if i1D:
        lPData  = r.RooDataHist(iLabel+"_pass_"+iBin,iLabel+"_pass_"+iBin, r.RooArgList(iVars[0]),iH[0])
        lMData  = r.RooDataHist(iLabel+"_fail_"+iBin,iLabel+"_fail_"+iBin, r.RooArgList(iVars[0]),iH[1]) 
        lP      = r.RooHistPdf (iLabel+"passh"+iBin,iLabel+"passh"+iBin, r.RooArgList(iVars[4]),r.RooArgList(iVars[0]),lPData,1)
        lF      = r.RooHistPdf (iLabel+"failh"+iBin,iLabel+"failh"+iBin, r.RooArgList(iVars[4]),r.RooArgList(iVars[0]),lMData,1)
    else:
        lPData  = r.RooDataHist(iLabel+"_pass_" +iBin,iLabel+"_pass_"+iBin, r.RooArgList(iVars[0],iVars[5]),iH[0])
        lMData  = r.RooDataHist(iLabel+"_fail_" +iBin,iLabel+"_fail_"+iBin, r.RooArgList(iVars[0],iVars[5]),iH[1]) 
        lP      = r.RooHistPdf (iLabel+"passh" +iBin,iLabel+"passh"+iBin, r.RooArgList(iVars[4],iVars[5]),r.RooArgList(iVars[0],iVars[5]),lPData,0)
        lF      = r.RooHistPdf (iLabel+"failh" +iBin,iLabel+"failh"+iBin, r.RooArgList(iVars[4],iVars[5]),r.RooArgList(iVars[0],iVars[5]),lMData,0)
    lEP     = r.RooExtendPdf(iLabel+"_passe_" +iBin,iLabel+"pe" +iBin,lP,lNPass)
    lEF     = r.RooExtendPdf(iLabel+"_faile_" +iBin,iLabel+"fe" +iBin,lF,lNFail)
    lHist   = [lP,lF,lEP,lEF,lPData,lMData]
    fVars.extend([lNTot,lNPass,lNFail])
    fFuncs.extend(lHist)
    return lHist
         
def fit1D(iHP,iHF,iFail,iBin="_cat0"):
    lBase  = baseVars()
    if iFail:
        lData = r.RooDataHist("data_obs"+iBin,"data_obs"+iBin,r.RooArgList(lBase[0]),iHF[0])
    else:
        lData = r.RooDataHist("data_obs"+iBin,"data_obs"+iBin,r.RooArgList(lBase[0]),iHP[0])
    lW    = histFunc([iHP[1],iHF[1]],lBase,"wqq",iBin)
    lZ    = histFunc([iHP[2],iHF[2]],lBase,"zqq",iBin)
    lQCD  = qcdFunc ([iHP[3],iHF[3]],lBase,iBin)
    lTot  = r.RooAddPdf("tot"+iBin,"tot"+iBin,r.RooArgList(lQCD[0+iFail],lW[2+iFail],lZ[2+iFail]))
    lTot.fitTo(lData,r.RooFit.Extended(),r.RooFit.Minos())
    workspace([lData],fFuncs,fVars,iBin)
    writeHist(iHP,lBase[3],[lW[2+iFail],lZ[2+iFail]])
    lFuncs=[lTot,lQCD[0+iFail]]#,lW[2+iFail],lZ[2+iFail]]
    draw(lBase[0],lData,lFuncs)
    
def fit1DPF(iHP,iHF,iBin="_cat0"):
    lCats = r.RooCategory("sample","sample") 
    lCats.defineType("pass",1) 
    lCats.defineType("fail",0) 
    lBase  = baseVars()
    lPData = r.RooDataHist("data_obs_pass"+iBin,"data_obs_pass"+iBin,r.RooArgList(lBase[0]),iHP[0])
    lFData = r.RooDataHist("data_obs_fail"+iBin,"data_obs_fail"+iBin,r.RooArgList(lBase[0]),iHF[0])
    lData  = r.RooDataHist("comb_data_obs","comb_data_obs",r.RooArgList(lBase[0]),r.RooFit.Index(lCats),r.RooFit.Import("pass",lPData),r.RooFit.Import("fail",lFData))
    lW    = histFunc([iHP[1],iHF[1]],lBase,"wqq",iBin)
    lZ    = histFunc([iHP[2],iHF[2]],lBase,"zqq",iBin)
    lQCD  = qcdFunc ([iHP[3],iHF[3]],lBase,iBin)
    lTotP = r.RooAddPdf("tot_pass"+iBin,"tot_pass"+iBin,r.RooArgList(lQCD[0],lW[2],lZ[2]))
    lTotF = r.RooAddPdf("tot_fail"+iBin,"tot_fail"+iBin,r.RooArgList(lQCD[1],lW[3],lZ[3]))
    lTot  = r.RooSimultaneous("tot","tot",lCats) 
    lTot.addPdf(lTotP,"pass") 
    lTot.addPdf(lTotF,"fail")     
    lTot.fitTo(lData,r.RooFit.Extended(),r.RooFit.Minos())
    lDatas=[lPData,lFData]
    workspace(lDatas,fFuncs,fVars,iBin)
    lPFuncs=[lTotP,lQCD[0]]#,lW[2],lZ[2]]
    lFFuncs=[lTotF,lQCD[1]]#,lW[3],lZ[3]]
    drawPF(lBase[0],[lPData,lFData],[lPFuncs,lFFuncs])
        
def fit2D(iHP,iHF,iFail,iBin="_cat0"):
    lBase  = baseVars(False)
    if iFail:
        lData = r.RooDataHist("data_obs"+iBin,"data_obs"+iBin,r.RooArgList(lBase[0],lBase[5]),iHF[0])
    else:
        lData = r.RooDataHist("data_obs"+iBin,"data_obs"+iBin,r.RooArgList(lBase[0],lBase[5]),iHP[0])
    lW    = histFunc([iHP[1],iHF[1]],lBase,"wqq",iBin,False)
    lZ    = histFunc([iHP[2],iHF[2]],lBase,"zqq",iBin,False)
    lQCD  = qcdFunc ([iHP[3],iHF[3]],lBase,iBin,False)
    lTot  = r.RooAddPdf("tot"+iBin,"tot"+iBin,r.RooArgList(lQCD[0+iFail],lW[2+iFail],lZ[2+iFail]))
    lTot.fitTo(lData,r.RooFit.Extended(),r.RooFit.Minos(),r.RooFit.ConditionalObservables(r.RooArgSet(lBase[5])))
    workspace([lData],fFuncs,fVars,iBin)
    writeHist(iHP,lBase[3],[lW[2+iFail],lZ[2+iFail]])
    lFuncs=[lTot,lQCD[0+iFail]]#,lW[2+iFail],lZ[2+iFail]]
    draw(lBase[0],lData,lFuncs)

def fit2DPF(iHP,iHF,iBin="_cat0"):
    lBase  = baseVars(False)
    lCats = r.RooCategory("sample","sample") 
    lCats.defineType("pass",1) 
    lCats.defineType("fail",0) 
    lPData = r.RooDataHist("data_obs_pass"+iBin,"data_obs_pass"+iBin,r.RooArgList(lBase[0],lBase[5]),iHP[0])
    lFData = r.RooDataHist("data_obs_fail"+iBin,"data_obs_fail"+iBin,r.RooArgList(lBase[0],lBase[5]),iHF[0])
    lData  = r.RooDataHist("comb_data_obs","comb_data_obs",r.RooArgList(lBase[0],lBase[5]),r.RooFit.Index(lCats),r.RooFit.Import("pass",lPData),r.RooFit.Import("fail",lFData))
    lW    = histFunc([iHP[1],iHF[1]],lBase,"wqq",iBin,False)
    lZ    = histFunc([iHP[2],iHF[2]],lBase,"zqq",iBin,False)
    lQCD  = qcdFunc ([iHP[3],iHF[3]],lBase,iBin,False)
    lTotP = r.RooAddPdf("tot_pass"+iBin,"tot_pass"+iBin,r.RooArgList(lQCD[0],lW[2],lZ[2]))
    lTotF = r.RooAddPdf("tot_fail"+iBin,"tot_fail"+iBin,r.RooArgList(lQCD[1]))#,lW[3],lZ[3]))
    lTot  = r.RooSimultaneous("tot","tot",lCats) 
    lTot.addPdf(lTotP,"pass") 
    lTot.addPdf(lTotF,"fail")     
    bestFit(fPars[0],fPars[1],fPars[2],fPars[3],fPars[4],fPars[5],fPars[6],fPars[7],fPars[8])
    lTot.fitTo(lData,r.RooFit.Extended(),r.RooFit.ConditionalObservables(r.RooArgSet(lBase[5])))
    workspace([lPData,lFData],fFuncs,fVars,iBin)
    lPFuncs=[lTotP,lQCD[0]]
    lFFuncs=[lTotF,lQCD[1]]
    drawPF(lBase[0],[lPData,lFData],[lPFuncs,lFFuncs])

def proj(iLabel,iBin,iH):
    lH = r.TH1F(iH.GetName()+"_"+iLabel+iBin,iH.GetName()+"_"+iLabel+iBin,iH.GetNbinsX(),iH.GetXaxis().GetXmin(),iH.GetXaxis().GetXmax())
    for iM in range(1,iH.GetNbinsX()+1):
        lH.SetBinContent(iM,iH.GetBinContent(iM,int(iBin)))
    lH.SetDirectory(0)
    return lH
    
def cat(iCat,iHP,iHF,iBin,iBase,iPt):
    iCat.defineType("pass"+iBin)
    iCat.defineType("fail"+iBin)
    lPData = r.RooDataHist("data_obs_pass_"+iBin,"data_obs_pass_"+iBin,r.RooArgList(iBase[0]),iHP[0])
    lFData = r.RooDataHist("data_obs_fail_"+iBin,"data_obs_fail_"+iBin,r.RooArgList(iBase[0]),iHF[0])
    lW    = histFunc([iHP[1],iHF[1]],iBase,"wqq",iBin)
    lZ    = histFunc([iHP[2],iHF[2]],iBase,"zqq",iBin)
    lQCD  = qcdFunc ([iHP[3],iHF[3]],iBase,iBin,True,iPt)
    lTotP = r.RooAddPdf("tot_pass"+iBin,"tot_pass"+iBin,r.RooArgList(lQCD[0]))#,lW[2],lZ[2]))
    lTotF = r.RooAddPdf("tot_fail"+iBin,"tot_fail"+iBin,r.RooArgList(lQCD[1]))#,lW[3],lZ[3]))
    fDatas.extend([lPData,lFData])
    fFuncs.extend([lTotP,lTotF])
    return ([lPData,lFData],[lTotP,lTotF],[lW[4],lZ[4],lW[5],lZ[5]])
    
def fitCat(iHP,iHF,iBin="_cat0"):
    lBase = baseVars(False)
    lCats   = r.RooCategory("sample","sample") 
    lBase   = baseVars()
    lPtBins = iHP[0].GetNbinsY()
    lPdfs   = []
    lRDatas = []
    lDatas  = r.std.map ('string, RooDataHist*')()
    for pt in range(1,lPtBins+1):
        lPCat = []
        lFCat = []
        for pH in iHP:
            lHP = proj(iBin,str(pt),pH)
            lPCat.append(lHP)
        for pH in iHF:
            lHF = proj(iBin,str(pt),pH)
            lFCat.append(lHF)
        (pDatas,pPdfs,pHists) = cat(lCats,lPCat,lFCat,"_"+str(pt),lBase,iHP[0].GetYaxis().GetBinCenter(pt))
        #workspace(pDatas,fFuncs,fVars,pt)
        lPair0 = r.std.pair('string, RooDataHist*')("pass_"+str(pt),pDatas[0])
        lPair1 = r.std.pair('string, RooDataHist*')("fail_"+str(pt),pDatas[1])
        lDatas.insert(lPair0)
        lDatas.insert(lPair1)
        for i0 in range(0,len(pPdfs)):
            lPdfs .append(pPdfs[i0])
            lRDatas.append(pDatas[i0])
    print lCats.numTypes()
    lData  = r.RooDataHist("comb_data_obs","comb_data_obs",r.RooArgList(lBase[0]),lCats,lDatas)
    lPdf   = r.RooSimultaneous("tot","tot",lCats)
    for i0 in range(0,len(lPdfs)):
        lPdf .addPdf(lPdfs [i0],lPdfs[i0].GetName().replace("tot_",""))
    lPdf.fitTo(lData,r.RooFit.Extended(),r.RooFit.Minos())
    for i0 in range(0,lPtBins):
        drawPF(lBase[0],[lRDatas[2*i0],lRDatas[2*i0+1]],[[lPdfs[2*i0]],[lPdfs[2*i0+1]]],i0)

def dumpRalph(iH,iBase,iPt,iCat):
    lName=(((iH.GetName().replace("_fail_","")).replace("_pass_","")).replace(iCat,"")).replace("2D_","")
    iBase[5].setVal(iPt)
    lP1     = r.RooRealVar   ("p1" ,"p1", -0.001   ,-0.1  ,0.1)
    lP2     = r.RooRealVar   ("p2" ,"p2", 0.0   ,-0.001,0.001)
    lR1     = r.RooRealVar   ("r1" ,"r1", 0.0   ,-0.5 ,0.5)
    lR2     = r.RooRealVar   ("r2" ,"r2", 0.0   ,-0.1 ,0.1)
    lPassBins = r.RooArgList()
    lFailBins = r.RooArgList()
    for i0 in range(1,iH.GetNbinsX()+1):
        iBase[0].setVal(iH.GetXaxis().GetBinCenter(i0)) 
        lPass = command("qcdeff"    ,iPt,iBase[0].getVal(),iBase[6].getVal(),True)
        lFail = command("(1-qcdeff)",iPt,iBase[0].getVal(),iBase[6].getVal(),False)
        pFail = r.RooRealVar   (lName+"_fail_"+iCat+"_Bin"+str(i0),lName+"_fail_"+iCat+"_Bin"+str(i0),iH.GetBinContent(i0),0,50000000.);
        pPass = r.RooFormulaVar(lName+"_pass_"+iCat+"_Bin"+str(i0),lName+"_pass_"+iCat+"_Bin"+str(i0),lPass+"/("+lFail+")*@0",r.RooArgList(pFail,lP1,lR1,iBase[2]))
        lPassBins.add(pPass)
        lFailBins.add(pFail)
        fVars.extend([pPass,pFail])
        fPars.extend([pPass,pFail])
        print  pFail.GetName(),"flatParam"
        #lPassBin=lName+"_Bin"+str(i0)).replace("fail","pass")
        #lFailBin=lName+"_Bin"+str(i0)
        #print lPassBin,"rateParam",iCat,lPass+"/("+lFail+")*@0",lFailBin
        #print lFailBin,"rateParam",iCat,iH.GetBinContent(i0)
    lPass  = r.RooParametricHist(lName+"_pass_"+iCat,lName+"_pass_"+iCat,iBase[0],lPassBins,iH)
    lFail  = r.RooParametricHist(lName+"_fail_"+iCat,lName+"_fail_"+iCat,iBase[0],lFailBins,iH)
    lNPass = r.RooAddition(lName+"_pass_"+iCat+"_norm",lName+"_pass_"+iCat+"_norm",lPassBins)
    lNFail = r.RooAddition(lName+"_fail_"+iCat+"_norm",lName+"_fail_"+iCat+"_norm",lFailBins)
    fFuncs.extend([lPass,lFail,lNPass,lNFail])
    lWPass = r.RooWorkspace("w_pass_"+str(iCat))
    lWFail = r.RooWorkspace("w_fail_"+str(iCat))
    getattr(lWPass,'import')(lPass,r.RooFit.RecycleConflictNodes())
    getattr(lWPass,'import')(lNPass,r.RooFit.RecycleConflictNodes())
    getattr(lWFail,'import')(lFail,r.RooFit.RecycleConflictNodes())
    getattr(lWFail,'import')(lNFail,r.RooFit.RecycleConflictNodes())
    if iCat.find("1") > -1:
        lWPass.writeToFile("simple-shapes-paramHist.root")
    else:
        lWPass.writeToFile("simple-shapes-paramHist.root",False)
    lWFail.writeToFile("simple-shapes-paramHist.root",False)

def blank(iFront,iEnd,iH):
    lH = iH.Clone(iFront+iEnd)
    for i0 in range(0,iH.GetNbinsX()+1):
        lH.SetBinContent(i0,1.)
    return lH

#pt categories with rhalpha Fit constraints does on as RateParams
def rhoCard(iHP,iHF,iBin="cat"):
    lBase = baseVars(False)
    lCats   = r.RooCategory("sample","sample") 
    lPtBins = iHP[0].GetNbinsY()
    lHists  = []
    lRDatas = []
    lBlanks = []
    for pt in range(1,lPtBins+1):
        lPCat = []
        lFCat = []
        for pH in iHP:
            lHP = proj(iBin,str(pt),pH)
            lPCat.append(lHP)
        for pH in iHF:
            lHF = proj(iBin,str(pt),pH)
            lFCat.append(lHF)
        (pDatas,pPdfs,pHists) = cat(lCats,lPCat,lFCat,iBin+str(pt),lBase,iHP[0].GetYaxis().GetBinCenter(pt))
        for i0 in range(0,len(pPdfs)):
            lHists .extend(pHists)
            lRDatas.append(pDatas[i0])
        pPHists = [pHists[0],pHists[1]]
        pFHists = [pHists[2],pHists[3]]
        lPCat.append(blank("bkg_pass_",iBin+str(pt),lPCat[0]))
        lFCat.append(blank("bkg_fail_",iBin+str(pt),lFCat[0]))
        dumpRalph(lFCat[len(lFCat)-2],lBase,iHP[0].GetYaxis().GetBinCenter(pt),iBin+str(pt))
        #lSigs.extend([lW[2],lZ[2]])
        #lBlanks.extend([lPCat[len(lPCat)-1],lFCat[len(lFCat)-1]])
        workspace([pDatas[0]],pPHists,fVars,"pass_"+iBin+str(pt),True)
        workspace([pDatas[1]],pFHists,fVars,"fail_"+iBin+str(pt),True)
    #print lSigs[0].GetName()
    #writeHist(lRDatas,lBase[3],lSigs,lBlanks)  

def load(iFileName,iHP,iHF,i2D):
    end=""
    if i2D:
        end="_2D"
    lFile = r.TFile(iFileName)
    lHP0 = lFile.Get("data_pass"+end)
    lHF0 = lFile.Get("data_fail"+end)
    lHP1 = lFile.Get("wqq_pass" +end)
    lHF1 = lFile.Get("wqq_fail" +end)
    lHP2 = lFile.Get("zqq_pass" +end)
    lHF2 = lFile.Get("zqq_fail" +end)
    lHP3 = lFile.Get("qcd_pass" +end)
    lHF3 = lFile.Get("qcd_fail" +end)
    iHP.extend([lHP0,lHP1,lHP2,lHP3])
    iHF.extend([lHF0,lHF1,lHF2,lHF3])
    for lH in (iHP+iHF):
        lH.SetDirectory(0)
    return 
    
if __name__ == "__main__":
    options = parser()
    print options
    lHP =[]
    lHF =[]
    load(options.input,lHP,lHF,options.fit2D or (options.cat and not options.fit1D) or (options.card and not options.fit1D))
    if options.fit1D and not options.cat and not options.card:
        if options.passfail:
            fit1DPF(lHP,lHF)
        else:
            fit1D(lHP,lHF,options.fail)

    if options.fit2D:
        if options.passfail:
            fit2DPF(lHP,lHF)
        else:
            fit2D(lHP,lHF,options.fail)

    if options.cat:
        fitCat(lHP,lHF)
    
    if options.card:
        rhoCard(lHP,lHF)
