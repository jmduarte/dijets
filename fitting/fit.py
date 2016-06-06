#!/usr/bin/env python
import ROOT as r,sys,math,array,os
from optparse import OptionParser
from ROOT import std,RooDataHist
r.gROOT.Macro(os.path.expanduser('~/rootlogon.C'))
r.gInterpreter.GenerateDictionary("std::pair<std::string, RooDataHist*>", "map;string;RooDataHist.h")
r.gInterpreter.GenerateDictionary("std::map<std::string, RooDataHist*>", "map;string;RooDataHist.h")

fOutput="bern.root"
fNBins=68
#fNBins=30
fXMin=28
fXMax=240
f2D=False
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
    parser.add_option('--input'   ,action='store',type='string',dest='input'   ,default='hists.root',help='input file')
    parser.add_option('--func'    ,action='store',type='int'   ,dest='func'    ,default=0           ,help='n=ith order of Bernstein')
    parser.add_option('--output'  ,action='store',type='string',dest='output'  ,default='bern.root' ,help='workspace output')
    parser.add_option('--xMin'    ,action='store',type='float' ,dest='xmin'    ,default=28          ,help='x-min')
    parser.add_option('--xMax'    ,action='store',type='float' ,dest='xmax'    ,default=300         ,help='x-max')
    parser.add_option('--nBins'   ,action='store',type='int'   ,dest='nbins'    ,default=68          ,help='n-bins')
    (options,args) = parser.parse_args()
    return options

def drawFrame(iFrame,iData,iFuncs):
    iData.plotOn(iFrame)
    iColor=50
    for pFunc in iFuncs:
        #pFunc.plotOn(iFrame,r.RooFit.LineColor(iColor),r.RooFit.LineStyle(iColor != 50+1),r.RooFit.Components(pFunc.GetName()),r.RooFit.ProjWData(iData))
        #pFunc.plotOn(iFrame,r.RooFit.LineColor(iColor),r.RooFit.LineStyle(iColor != 50+1),r.RooFit.ProjWData(iData))
        pFunc.plotOn(iFrame,r.RooFit.LineColor(iColor),r.RooFit.LineStyle(iColor != 50+1))#,r.RooFit.ProjWData(iData))
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
    lInt    = iDataHist.createHistogram("x").Integral()
    lDM     = r.RooRealVar   ("Xdm","Xdm", 0.,-10,10)
    lShift  = r.RooFormulaVar("Xshift",iVar.GetName()+"-Xdm",r.RooArgList(iVar,lDM))  
    if f2D:
        lSPdf   = r.RooHistPdf(iDataHist.GetName()+"P",iDataHist.GetName()+"P", r.RooArgList(lShift,fVars[1]),r.RooArgList(iVar,fVars[1]),iDataHist,0)
    else:
        lSPdf   = r.RooHistPdf(iDataHist.GetName()+"P",iDataHist.GetName()+"P", r.RooArgList(lShift),r.RooArgList(iVar),iDataHist,0)
    lDM.setVal(iShift)
    lHUp   = lSPdf.createHistogram("x")
    lHUp.Scale(lInt)
    lUp    = r.RooDataHist(iDataHist.GetName()+"_scaleUp",iDataHist.GetName()+"_scaleUp", r.RooArgList(iVar),lHUp)
    lDM.setVal(-iShift)
    lHDown = lSPdf.createHistogram("x")
    lHDown.Scale(lInt)
    lDown  = r.RooDataHist(iDataHist.GetName()+"_scaleDown",iDataHist.GetName()+"_scaleDown", r.RooArgList(iVar),lHDown)
    return (lUp,lDown)

def smear(iVar,iDataHist,iScale=0.1):
    lDM     = r.RooRealVar("Xshift","Xshift", 1.,0.,2.)
    lVar    = iDataHist.createHistogram("x").GetMean()
    lInt    = iDataHist.createHistogram("x").Integral()
    lShift  = r.RooFormulaVar("Xsmear","("+iVar.GetName()+"-"+str(lVar)+")/Xshift+"+str(lVar),r.RooArgList(iVar,lDM))  
    if f2D:
        lHPdf   = r.RooHistPdf(iDataHist.GetName()+"S",iDataHist.GetName()+"S", r.RooArgList(lShift,fVars[1]),r.RooArgList(iVar,fVars[1]),iDataHist,0)
    else:
        lHPdf   = r.RooHistPdf(iDataHist.GetName()+"S",iDataHist.GetName()+"S", r.RooArgList(lShift),r.RooArgList(iVar),iDataHist,0)
    lDM.setVal(1.+iScale)
    lHUp = lHPdf.createHistogram("x")
    lHUp.Scale(lInt)
    lUp = r.RooDataHist(iDataHist.GetName()+"_smearUp",iDataHist.GetName()+"_smearUp", r.RooArgList(iVar),lHUp)    
    lDM.setVal(1.-iScale)
    lHDown = lHPdf.createHistogram("x")
    lHDown.Scale(lInt)
    lDown  = r.RooDataHist(iDataHist.GetName()+"_smearDown",iDataHist.GetName()+"_smearDown", r.RooArgList(iVar),lHDown)
    return [lUp,lDown]    

def workspace(iDatas,iFuncs,iVars,iCat="cat0",iShift=True):
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
        if iShift and pFunc.GetName().find("qq") > -1:
            (pFUp, pFDown)  = shift(iVars[0],pFunc,5.)
            (pSFUp,pSFDown) = smear(iVars[0],pFunc,0.05)
            getattr(lW,'import')(pFUp,  r.RooFit.RecycleConflictNodes())
            getattr(lW,'import')(pFDown,r.RooFit.RecycleConflictNodes())
            getattr(lW,'import')(pSFUp,  r.RooFit.RecycleConflictNodes())
            getattr(lW,'import')(pSFDown,r.RooFit.RecycleConflictNodes())
    
    if iCat.find("pass_cat1") == -1:
        lW.writeToFile(fOutput,False)
    else:
        lW.writeToFile(fOutput)

def writeHist(iH,iShift,iTemps,iBlanks=[]):  
    lOut = r.TFile(fOutput,"RECREATE")
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

def baseVars(i1D=False):
    global fNBins,fXMin,fXMax
    lMSD    = r.RooRealVar("x","x",fXMin,fXMax)
    lMSD.setBins(fNBins)
    lPt  = r.RooRealVar   ("pt","pt",500,1500)
    lRho = r.RooFormulaVar("rho","log(x*x/pt)",r.RooArgList(lMSD,lPt))
    lEff    = r.RooRealVar("veff"      ,"veff"      ,0.5 ,0.,1.0)
    lEffQCD = r.RooRealVar("qcdeff"    ,"qcdeff"   ,0.01,0.,10.)
    lDM     = r.RooRealVar("dm","dm", 0.,-10,10)
    lShift  = r.RooFormulaVar("shift",lMSD.GetName()+"-dm",r.RooArgList(lMSD,lDM))  
    lVars=[lMSD,lEff,lEffQCD,lDM,lShift]
    if not i1D:
        lVars.extend([lPt,lRho])
    fVars.extend([lMSD,lPt,lEff,lEffQCD,lDM])
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
        
def command(iBase,iPt,iM,iRho,iPass=True,iNR=1,iNP=1,iNRNP=1):
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
    
    if iPt > 0 and iNRNP > 0:
        lRho  = '(log(x*x/'+str(iPt)+')-2.5)'           if iM < 0 else str(iRho-2.5)
        lPRho = '(log(x*x/'+str(iPt)+')-2.5)*'+str(iPt) if iM < 0 else str((iRho-2.5)*(iPt-500))
        lFunc+=("*(1+pr1*"+str(lPRho))
        if iNRNP > 1:
            lFunc+=("+pr2*"+str(iPt-500)+"*"+lPRho)
        if iNRNP > 2:
            lFunc+=("+pr3*"+lPRho+"*"+lRho)
        if iNRNP > 3:
            lFunc+=("+pr4*"+lPRho+"*"+lPRho)
        lFunc+=")"

    if not iPass:
        lFunc=lFunc.replace("+p","-p")
        lFunc=lFunc.replace("+r","-r")
    return lFunc

def qcdFunc(iH,iVars,iBin="_cat0",iFunc=0,i1D=True,iPt=-1):
    lNTot   = r.RooRealVar   ("qcdnorm"+iBin,"qcdnorm"+iBin,(iH[0].Integral()+iH[1].Integral()),0.,3000.*(iH[0].Integral()+iH[1].Integral()))
    lNPass  = r.RooFormulaVar("fqpass"+iBin ,"qcdnorm"+iBin+"*(qcdeff)"  ,r.RooArgList(lNTot,iVars[2]))
    lNFail  = r.RooFormulaVar("fqfail"+iBin ,"qcdnorm"+iBin+"*(1-qcdeff)",r.RooArgList(lNTot,iVars[2]))
    lA0     = r.RooRealVar   ("a0"          ,"a0"          ,0.486,-1000.,1000.)          
    lA1     = r.RooRealVar   ("a1"          ,"a1"          ,0.0,-1  ,1.)
    lA2     = r.RooRealVar   ("a2"          ,"a2"          ,0.0,-0.1,0.1)
    lA3     = r.RooRealVar   ("a3"          ,"a3"          ,0.0,-0.1,0.1)
    lA4     = r.RooRealVar   ("a4"          ,"a4"          ,0.0,-0.1,0.1)
    lA5     = r.RooRealVar   ("a5"          ,"a5"          ,0.0,-0.1,0.1)
    lA6     = r.RooRealVar   ("a6"          ,"a6"          ,0.0,-0.1,0.1)
    lSigma1 = r.RooRealVar   ("sigma1"      ,"sigma1"      ,50,0,1000)
    lP1     = r.RooRealVar   ("p1" ,"p1", 0.0   ,-1.0  ,1.0)
    lP2     = r.RooRealVar   ("p2" ,"p2", 0.0   ,-0.001,0.001)
    lR1     = r.RooRealVar   ("r1" ,"r1", 0.0   ,-1.0 ,1.0)
    lR2     = r.RooRealVar   ("r2" ,"r2", 0.0   ,-0.1 ,0.1)
    if i1D:
        lFunc="1"
        lQFuncP1 = 0
        lQFuncF1 = 0
        if iFunc == 3:
            lQFuncP1 = r.RooBernstein("qcd_pass_"+iBin,"qcd_pass_"+iBin,iVars[0],r.RooArgList(lA0,lA1,lA2,lA3))
            lQFuncF1 = r.RooBernstein("qcd_fail_"+iBin,"qcd_fail_"+iBin,iVars[0],r.RooArgList(lA0,lA1,lA2,lA3))
        elif iFunc == 4:
            lQFuncP1 = r.RooBernstein("qcd_pass_"+iBin,"qcd_pass_"+iBin,iVars[0],r.RooArgList(lA0,lA1,lA2,lA3,lA4))
            lQFuncF1 = r.RooBernstein("qcd_fail_"+iBin,"qcd_fail_"+iBin,iVars[0],r.RooArgList(lA0,lA1,lA2,lA3,lA4))
        elif iFunc == 5:
            lQFuncP1 = r.RooBernstein("qcd_pass_"+iBin,"qcd_pass_"+iBin,iVars[0],r.RooArgList(lA0,lA1,lA2,lA3,lA4,lA5))
            lQFuncF1 = r.RooBernstein("qcd_fail_"+iBin,"qcd_fail_"+iBin,iVars[0],r.RooArgList(lA0,lA1,lA2,lA3,lA4,lA5))
        elif iFunc == 6:
            lQFuncP1 = r.RooBernstein("qcd_pass_"+iBin,"qcd_pass_"+iBin,iVars[0],r.RooArgList(lA0,lA1,lA2,lA3,lA4,lA5,lA6))
            lQFuncF1 = r.RooBernstein("qcd_fail_"+iBin,"qcd_fail_"+iBin,iVars[0],r.RooArgList(lA0,lA1,lA2,lA3,lA4,lA5,lA6))
        elif iFunc == 2:
            lQFuncP1 = r.RooBernstein("qcd_pass_"+iBin,"qcd_pass_"+iBin,iVars[0],r.RooArgList(lA0,lA1,lA2))
            lQFuncF1 = r.RooBernstein("qcd_fail_"+iBin,"qcd_fail_"+iBin,iVars[0],r.RooArgList(lA0,lA1,lA2))
        elif iFunc == 1:
            lQFuncP1 = r.RooBernstein("qcd_pass_"+iBin,"qcd_pass_"+iBin,iVars[0],r.RooArgList(lA0,lA1))
            lQFuncF1 = r.RooBernstein("qcd_fail_"+iBin,"qcd_fail_"+iBin,iVars[0],r.RooArgList(lA0,lA1))
        else:
            lQFuncP1 = r.RooBernstein("qcd_pass_"+iBin,"qcd_pass_"+iBin,iVars[0],r.RooArgList(lA0))
            lQFuncF1 = r.RooBernstein("qcd_fail_"+iBin,"qcd_fail_"+iBin,iVars[0],r.RooArgList(lA0))
        if iPt > 0:
            iVars[5].setVal(iPt)
            lPFunc = command(lFunc,iPt,-1,-1,False,2,2)
            lFFunc = command(lFunc,iPt,-1,-1,True ,2,2)
            lQFuncP1.SetName("qcd_pass1_"+iBin)
            lQFuncF1.SetName("qcd_fail1_"+iBin)
            lQFuncP1.SetTitle("qcd_pass1_"+iBin)
            lQFuncF1.SetTitle("qcd_fail1_"+iBin)
            lQFuncP2 = r.RooGenericPdf("qcd_pass2_"+iBin,lPFunc,r.RooArgList(iVars[0],lP1,lR1,lP2,lR2))
            lQFuncF2 = r.RooGenericPdf("qcd_fail2_"+iBin,lFFunc,r.RooArgList(iVars[0],lP1,lR1,lP2,lR2))
            lQFuncP  = r.RooProdPdf("qcd_pass_"+iBin,"qcd_pass_"+iBin,lQFuncP1,lQFuncP2)
            lQFuncF  = r.RooProdPdf("qcd_fail_"+iBin,"qcd_fail_"+iBin,lQFuncF1,lQFuncF2)
            fFuncs.extend([lQFuncP1,lQFuncF1,lQFuncP2,lQFuncF2])
        else:
            lQFuncP = lQFuncP1
            lQFuncF = lQFuncF1
    else:
        lQFuncP = r.RooGenericPdf("qcd_pass_"+iBin,"("+lFunc+"*(1+p1*(pt-500)+r1*(rho-2.5))",r.RooArgList(iVars[0],lSigma1,lA0,lA1,lA2,lP1,iVars[5],lR1,iVars[6]))
        lQFuncF = r.RooGenericPdf("qcd_fail_"+iBin,"("+lFunc+"*(1-p1*(pt-500)-r1*(rho-2.5))",r.RooArgList(iVars[0],lSigma1,lA0,lA1,lA2,lP1,iVars[5],lR1,iVars[6]))
    lQCDP   = r.RooExtendPdf ("qcdpassE"+iBin,"qcdpass"+iBin,lQFuncP,lNPass)
    lQCDF   = r.RooExtendPdf ("qcdfailE"+iBin,"qcdfail"+iBin,lQFuncF,lNFail)
    lQCD    = [lQCDP,lQCDF,lQFuncP,lQFuncF]
    fVars.extend([lNTot,lA0,lA1,lA2,lR1,lSigma1,lNPass,lNFail,lP1,lR1,lR2,lP2,lA3,lA4,lA5,lA6])
    fFuncs.extend(lQCD)
    fPars.extend([lA0,lA1,lA2,lP1,lR1,lP2,lR2,lSigma1])
    lP1.setConstant(True)
    lR1.setConstant(True)
    return lQCD

def histFunc(iH,iVars,iLabel="w",iBin="_cat0",i1D=True):
    lNTot   = r.RooRealVar (iLabel+"norm"+iBin,iLabel+"norm"+iBin,(iH[0].Integral()+iH[1].Integral()),0.,5.*(iH[0].Integral()+iH[1].Integral()))
    lNPass  = r.RooFormulaVar(iLabel+"fpass"+iBin ,iLabel+"norm"+iBin+"*(veff)"  ,r.RooArgList(lNTot,iVars[1]))
    lNFail  = r.RooFormulaVar(iLabel+"fqail"+iBin ,iLabel+"norm"+iBin+"*(1-veff)",r.RooArgList(lNTot,iVars[1]))
    if i1D:
        lPData  = r.RooDataHist(iLabel+"_pass_"+iBin,iLabel+"_pass_"+iBin, r.RooArgList(iVars[0]),iH[0])
        lMData  = r.RooDataHist(iLabel+"_fail_"+iBin,iLabel+"_fail_"+iBin, r.RooArgList(iVars[0]),iH[1]) 
        lP      = r.RooHistPdf (iLabel+"passh"+iBin,iLabel+"passh"+iBin, r.RooArgList(iVars[4]),r.RooArgList(iVars[0]),lPData,0)
        lF      = r.RooHistPdf (iLabel+"failh"+iBin,iLabel+"failh"+iBin, r.RooArgList(iVars[4]),r.RooArgList(iVars[0]),lMData,0)
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

def getSignals(iHP,iHF,iBin,iBase):
    lPSigs = []
    lFSigs = []
    masses=[50,75,100,125,150,200,250,300]
    for i0 in range(0,len(masses)):
        lSig = histFunc([iHP[i0+3],iHF[i0+3]],iBase,"zqq"+str(masses[i0]),iBin)
        lPSigs.append(lSig[4])
        lFSigs.append(lSig[5])
    return (lPSigs,lFSigs)

def fit1D(iHP,iHF,iFail,iFunc=0,iBin="cat0"):
    lBase  = baseVars()
    if iFail:
        lData = r.RooDataHist("data_obs_"+iBin,"data_obs_"+iBin,r.RooArgList(lBase[0]),iHF[0])
    else:
        lData = r.RooDataHist("data_obs_"+iBin,"data_obs_"+iBin,r.RooArgList(lBase[0]),iHP[0])
    lW    = histFunc([iHP[1],iHF[1]],lBase,"wqq",iBin)
    lZ    = histFunc([iHP[2],iHF[2]],lBase,"zqq",iBin)
    lQCD  = qcdFunc ([iHP[len(iHP)-1],iHF[len(iHF)-1]],lBase,iBin,iFunc)
    lTot  = r.RooAddPdf("tot"+iBin,"tot"+iBin,r.RooArgList(lQCD[0+iFail]))#,lW[2+iFail],lZ[2+iFail]))
    lTot.fitTo(lData,r.RooFit.Extended(),r.RooFit.Minos())
    lHists=[lW[4+iFail],lZ[4+iFail],lQCD[0+iFail]]
    lHists.extend(getSignals(iHP,iHF,iBin,lBase)[iFail])
    workspace([lData],lHists,fVars,iBin)
    lFuncs=[lTot,lQCD[0+iFail]]#,lW[2+iFail],lZ[2+iFail]]
    draw(lBase[0],lData,lFuncs)
    
def fit1DPF(iHP,iHF,iFunc=0,iBin="cat0"):
    lCats = r.RooCategory("sample","sample") 
    lCats.defineType("pass",1) 
    lCats.defineType("fail",0) 
    lBase  = baseVars()
    lPData = r.RooDataHist("data_obs_pass_"+iBin,"data_obs_pass_"+iBin,r.RooArgList(lBase[0]),iHP[0])
    lFData = r.RooDataHist("data_obs_fail_"+iBin,"data_obs_fail_"+iBin,r.RooArgList(lBase[0]),iHF[0])
    lData  = r.RooDataHist("comb_data_obs","comb_data_obs",r.RooArgList(lBase[0]),r.RooFit.Index(lCats),r.RooFit.Import("pass",lPData),r.RooFit.Import("fail",lFData))
    lW    = histFunc([iHP[1],iHF[1]],lBase,"wqq",iBin)
    lZ    = histFunc([iHP[2],iHF[2]],lBase,"zqq",iBin)
    lQCD  = qcdFunc ([iHP[len(iHP)-1],iHF[len(iHF)-1]],lBase,iBin,iFunc)
    lTotP = r.RooAddPdf("tot_pass"+iBin,"tot_pass"+iBin,r.RooArgList(lQCD[0],lW[2],lZ[2]))
    lTotF = r.RooAddPdf("tot_fail"+iBin,"tot_fail"+iBin,r.RooArgList(lQCD[1]))#,lW[3],lZ[3]))
    lTot  = r.RooSimultaneous("tot","tot",lCats) 
    lTot.addPdf(lTotP,"pass") 
    lTot.addPdf(lTotF,"fail")     
    lTot.fitTo(lData,r.RooFit.Extended(),r.RooFit.Minos())
    lHists=[lW[4],lZ[4],lQCD[0],lW[5],lZ[5],lQCD[1]]
    lHists.extend(getSignals(iHP,iHF,iBin,lBase)[0])
    lHists.extend(getSignals(iHP,iHF,iBin,lBase)[1])
    lDatas=[lPData,lFData]
    workspace(lDatas,lHists,fVars,iBin)
    lPFuncs=[lTotP,lQCD[0]]#,lW[2],lZ[2]]
    lFFuncs=[lTotF,lQCD[1]]#,lW[3],lZ[3]]
    drawPF(lBase[0],[lPData,lFData],[lPFuncs,lFFuncs])
        
def fit2D(iHP,iHF,iFail,iBin="cat0"):
    lBase  = baseVars(False)
    if iFail:
        lData = r.RooDataHist("data_obs_"+iBin,"data_obs_"+iBin,r.RooArgList(lBase[0],lBase[5]),iHF[0])
    else:
        lData = r.RooDataHist("data_obs_"+iBin,"data_obs_"+iBin,r.RooArgList(lBase[0],lBase[5]),iHP[0])
    lW    = histFunc([iHP[1],iHF[1]],lBase,"wqq",iBin,False)
    lZ    = histFunc([iHP[2],iHF[2]],lBase,"zqq",iBin,False)
    lQCD  = qcdFunc ([iHP[3],iHF[3]],lBase,iBin,False)
    lTot  = r.RooAddPdf("tot"+iBin,"tot"+iBin,r.RooArgList(lQCD[0+iFail],lW[2+iFail],lZ[2+iFail]))
    #lTot.fitTo(lData,r.RooFit.Extended(),r.RooFit.Minos(),r.RooFit.ConditionalObservables(r.RooArgSet(lBase[5])))
    workspace([lData],[lW[4+iFail],lZ[4+iFail],lQCD[2+iFail]],fVars,iBin)
    writeHist(iHP,lBase[3],[lW[4+iFail],lZ[4+iFail],lQCD[2+iFail]])
    lFuncs=[lTot,lQCD[0+iFail]]#,lW[2+iFail],lZ[2+iFail]]
    draw(lBase[0],lData,lFuncs)

def fit2DPF(iHP,iHF,iBin="_cat0"):
    lBase  = baseVars(False)
    lCats = r.RooCategory("sample","sample") 
    lCats.defineType("pass",1) 
    lCats.defineType("fail",0) 
    lPData = r.RooDataHist("data_obs_pass_"+iBin,"data_obs_pass"+iBin,r.RooArgList(lBase[0],lBase[5]),iHP[0])
    lFData = r.RooDataHist("data_obs_fail_"+iBin,"data_obs_fail"+iBin,r.RooArgList(lBase[0],lBase[5]),iHF[0])
    lData  = r.RooDataHist("comb_data_obs","comb_data_obs",r.RooArgList(lBase[0],lBase[5]),r.RooFit.Index(lCats),r.RooFit.Import("pass",lPData),r.RooFit.Import("fail",lFData))
    lW    = histFunc([iHP[1],iHF[1]],lBase,"wqq",iBin,False)
    lZ    = histFunc([iHP[2],iHF[2]],lBase,"zqq",iBin,False)
    lQCD  = qcdFunc ([iHP[3],iHF[3]],lBase,iBin,False)
    lTotP = r.RooAddPdf("tot_pass"+iBin,"tot_pass"+iBin,r.RooArgList(lQCD[0],lW[2],lZ[2]))
    lTotF = r.RooAddPdf("tot_fail"+iBin,"tot_fail"+iBin,r.RooArgList(lQCD[1]))#,lW[3],lZ[3]))
    lTot  = r.RooSimultaneous("tot","tot",lCats) 
    lTot.addPdf(lTotP,"pass") 
    lTot.addPdf(lTotF,"fail")     
    #bestFit(fPars[0],fPars[1],fPars[2],fPars[3],fPars[4],fPars[5],fPars[6],fPars[7],fPars[8])
    lTot.fitTo(lData,r.RooFit.Extended(),r.RooFit.ConditionalObservables(r.RooArgSet(lBase[5])))
    workspace([lPData,lFData],[lW[4],lZ[4],lQCD[2],lW[5],lZ[5],lQCD[3]],fVars,iBin)
    lPFuncs=[lTotP,lQCD[0]]
    lFFuncs=[lTotF,lQCD[1]]
    drawPF(lBase[0],[lPData,lFData],[lPFuncs,lFFuncs])

def proj(iLabel,iBin,iH):
    global fNBins,fXMin,fXMax
    lH = r.TH1F(iH.GetName()+"_"+iLabel+iBin,iH.GetName()+"_"+iLabel+iBin,fNBins,fXMin,fXMax)
    for iM in range(1,iH.GetNbinsX()+1):
        if iH.GetXaxis().GetBinCenter(iM) < lH.GetXaxis().GetXmin() or iH.GetXaxis().GetBinCenter(iM) > lH.GetXaxis().GetXmax():
            continue
        lH.SetBinContent(lH.FindBin(iH.GetXaxis().GetBinCenter(iM)),iH.GetBinContent(iM,int(iBin)))
    lH.SetDirectory(0)
    return lH
    
def cat(iCat,iHP,iHF,iBin,iBase,iPt,iFunc):
    iCat.defineType("pass"+iBin)
    iCat.defineType("fail"+iBin)
    lCats = r.RooCategory("Xsample","Xsample") 
    lCats.defineType("pass",1) 
    lCats.defineType("fail",0) 
    lPData = r.RooDataHist("data_obs_pass_"+iBin,"data_obs_pass_"+iBin,r.RooArgList(iBase[0]),iHP[0])
    lFData = r.RooDataHist("data_obs_fail_"+iBin,"data_obs_fail_"+iBin,r.RooArgList(iBase[0]),iHF[0])
    lData  = r.RooDataHist("comb_data_obs","comb_data_obs",r.RooArgList(iBase[0]),r.RooFit.Index(lCats),r.RooFit.Import("pass",lPData),r.RooFit.Import("fail",lFData))
    lW    = histFunc([iHP[1],iHF[1]],iBase,"wqq",iBin)
    lZ    = histFunc([iHP[2],iHF[2]],iBase,"zqq",iBin)
    lQCD  = qcdFunc ([iHP[len(iHP)-1],iHF[len(iHF)-1]],iBase,iBin,iFunc)#,True,iPt)
    lTotP = r.RooAddPdf("tot_pass"+iBin,"tot_pass"+iBin,r.RooArgList(lQCD[0]))#,lW[2],lZ[2]))
    lTotF = r.RooAddPdf("tot_fail"+iBin,"tot_fail"+iBin,r.RooArgList(lQCD[1]))#,lW[3],lZ[3]))
    lEWKP = r.RooAddPdf("ewk_pass"+iBin,"ewk_pass"+iBin,r.RooArgList(lW[2],lZ[2]))
    lEWKF = r.RooAddPdf("ewk_fail"+iBin,"ewk_fail"+iBin,r.RooArgList(lW[3],lZ[3]))
    lTot  = r.RooSimultaneous("tot","tot",lCats) 
    lTot.addPdf(lTotP,"pass") 
    lTot.addPdf(lTotF,"fail")     
    if iFunc > 0: lTot.fitTo(lData,r.RooFit.Extended(),r.RooFit.Minos()) 
    fDatas.extend([lPData,lFData])
    fFuncs.extend([lTotP,lTotF,lEWKP,lEWKF])
    return ([lPData,lFData],[lTotP,lTotF,lEWKP,lEWKF],[lW[4],lZ[4],lW[5],lZ[5]])
    
def fitCat(iHP,iHF,iFunc=0,iBin="_cat0"):
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
            print pH.GetName(),pH.Integral()
            lHP = proj(iBin,str(pt),pH)
            lPCat.append(lHP)
        for pH in iHF:
            print pH.GetName(),pH.Integral()
            lHF = proj(iBin,str(pt),pH)
            lFCat.append(lHF)
        (pDatas,pPdfs,pHists) = cat(lCats,lPCat,lFCat,"cat"+str(pt),lBase,iHP[0].GetYaxis().GetBinCenter(pt),iFunc)
        #for i0 in range(0,len(pDatas)):
        #    lPdfs .append(pPdfs[i0])
        #    lRDatas.append(pDatas[i0])
        #lPair0 = r.std.pair('string, RooDataHist*')("pass_"+str(pt),pDatas[0])
        #lPair1 = r.std.pair('string, RooDataHist*')("fail_"+str(pt),pDatas[1])
        #lDatas.insert(lPair0)
        #lDatas.insert(lPair1)
        lPHists=[pHists[0],pHists[1],pPdfs[0]]
        lFHists=[pHists[2],pHists[3],pPdfs[1]]
        lPHists.extend(getSignals(lPCat,lFCat,"cat"+str(pt),lBase)[0])
        lFHists.extend(getSignals(lPCat,lFCat,"cat"+str(pt),lBase)[1])
        workspace([pDatas[0]],lPHists,fVars,"pass_cat"+str(pt))
        workspace([pDatas[1]],lFHists,fVars,"fail_cat"+str(pt))
        print "Integral",lPCat[0].Integral(),lFCat[0].Integral(),"----->"
    return
    print len(lPdfs),lDatas["pass_1"]
    lData  = r.RooDataHist("comb_data_obs","comb_data_obs",r.RooArgList(lBase[0]),lCats,lDatas)
    lPdf   = r.RooSimultaneous("tot","tot",lCats)
    for i0 in range(0,len(lPdfs)):
        lPdf .addPdf(lPdfs [i0],lPdfs[i0].GetName().replace("tot_",""))
    lPdf.fitTo(lData,r.RooFit.Extended(),r.RooFit.Minos())
    #for i0 in range(0,lPtBins):
    #    drawPF(lBase[0],[lRDatas[2*i0],lRDatas[2*i0+1]],[[lPdfs[2*i0]],[lPdfs[2*i0+1]]],i0)

def dumpRalph(iH,iBase,iPt,iCat):
    #lName=(((iH.GetName().replace("_fail_","")).replace("_pass_","")).replace(iCat,"")).replace("2D_","")
    lName="qcd"
    iBase[5].setVal(iPt)
    lP1     = r.RooRealVar   ("p1" ,"p1", 0.0   ,-0.1  ,0.1)
    lP2     = r.RooRealVar   ("p2" ,"p2", 0.0   ,-0.001,0.001)
    lR1     = r.RooRealVar   ("r1" ,"r1", 0.0   ,-10.5 ,10.5)
    lR2     = r.RooRealVar   ("r2" ,"r2", 0.0   ,-10.5 ,10.5)
    lPR1    = r.RooRealVar   ("pr1","pr1", 0.0   ,-0.5 ,0.5)
    lPR2    = r.RooRealVar   ("pr2","pr2", 0.0   ,-0.1 ,0.1)
    lPR3    = r.RooRealVar   ("pr3","pr3", 0.0   ,-0.1 ,0.1)
    lPR4    = r.RooRealVar   ("pr4","pr4", 0.0   ,-0.1 ,0.1)
    lPassBins = r.RooArgList()
    lFailBins = r.RooArgList()
    iBase[2].setVal(0.02)
    iBase[2].setConstant(False)
    for i0 in range(1,iH.GetNbinsX()+1):
        iBase[0].setVal(iH.GetXaxis().GetBinCenter(i0)) 
        lPass = command("qcdeff"    ,iPt,iBase[0].getVal(),iBase[6].getVal(),True ,2,2,2)
        lFail = command("(1-qcdeff)",iPt,iBase[0].getVal(),iBase[6].getVal(),False,2,2,2)
        pFail = r.RooRealVar   (lName+"_fail_"+iCat+"_Bin"+str(i0),lName+"_fail_"+iCat+"_Bin"+str(i0),iH.GetBinContent(i0),-5,50000000.)
        #pPass = r.RooFormulaVar(lName+"_pass_"+iCat+"_Bin"+str(i0),lName+"_pass_"+iCat+"_Bin"+str(i0),lPass+"/("+lFail+")*"+pFail.GetName(),r.RooArgList(pFail,lP1,lR1,lP2,lR2,lPR1,lPR2,lPR3,iBase[2]))
        #pPass = r.RooFormulaVar(lName+"_pass_"+iCat+"_Bin"+str(i0),lName+"_pass_"+iCat+"_Bin"+str(i0),lPass+"*"+pFail.GetName(),r.RooArgList(pFail,lP1,lR1,iBase[2]))
        pPass = r.RooFormulaVar(lName+"_pass_"+iCat+"_Bin"+str(i0),lName+"_pass_"+iCat+"_Bin"+str(i0),lPass+"*"+pFail.GetName(),r.RooArgList(pFail,lP1,lR1,lP2,lR2,lPR1,lPR2,iBase[2]))
        #pFail.setConstant(True)
        lPassBins.add(pPass)
        lFailBins.add(pFail)
        fVars.extend([pPass,pFail])
        fPars.extend([pPass,pFail])
        print  pFail.GetName(),"flatParam"#,lPass+"/("+lFail+")*@0"
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
        lWPass.writeToFile("paramHist"+fOutput)
    else:
        lWPass.writeToFile("paramHist"+fOutput,False)
    lWFail.writeToFile("paramHist"+fOutput,False)
    return [lPass,lFail]

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
    lPdfs   = []
    lPPdfs  = []
    lHists  = []
    lRDatas = []
    lBlanks = []
    lDatas  = r.std.map ('string, RooDataHist*')()
    for pt in range(1,lPtBins+1):
        #fXmax = lXMax
        #if iHP[0].GetYaxis().GetBinCenter(pt) < 600: 
        #    fXMax = 270
        lPCat = []
        lFCat = []
        for pH in iHP:
            lHP = proj(iBin,str(pt),pH)
            lPCat.append(lHP)
        for pH in iHF:
            lHF = proj(iBin,str(pt),pH)
            lFCat.append(lHF)
        (pDatas,pPdfs,pHists) = cat(lCats,lPCat,lFCat,iBin+str(pt),lBase,iHP[0].GetYaxis().GetBinCenter(pt),0)
        for i0 in range(0,len(pPdfs)-2):
            lRDatas.append(pDatas[i0])
        lPPdfs.extend([pPdfs[0],pPdfs[1]])
        lParHists = dumpRalph(lFCat[0],lBase,iHP[0].GetYaxis().GetBinCenter(pt),iBin+str(pt))
        #pSPass,pSFail = sigcat(lCats,lPCat,lFCat,iBin+str(pt),lBase,iHP[0].GetYaxis().GetBinCenter(pt))
        lPHists=[pHists[0],pHists[1]]
        lFHists=[pHists[2],pHists[3]]
        lPHists.extend(getSignals(lPCat,lFCat,"cat"+str(pt),lBase)[0])
        lFHists.extend(getSignals(lPCat,lFCat,"cat"+str(pt),lBase)[1])
        workspace([pDatas[0]],lPHists,fVars,"pass_"+iBin+str(pt),True)
        workspace([pDatas[1]],lFHists,fVars,"fail_"+iBin+str(pt),True)
        lPair0 = r.std.pair('string, RooDataHist*')("pass_"+str(pt),pDatas[0])
        lPair1 = r.std.pair('string, RooDataHist*')("fail_"+str(pt),pDatas[1])
        lDatas.insert(lPair0)
        lDatas.insert(lPair1)
        pNCat   = r.RooRealVar ("bkg_"+iBin+str(pt)+"norm","bkg_"+iBin+str(pt)+"norm",iHF[3].Integral(),-1,500000.)
        pEPHist = r.RooExtendPdf(lParHists[0].GetName()+"E",lParHists[0].GetName()+"E",lParHists[0],pNCat)
        pEFHist = r.RooExtendPdf(lParHists[1].GetName()+"E",lParHists[1].GetName()+"E",lParHists[1],pNCat)
        pPass   = r.RooAddPdf("tot_pass"+iBin,"tot_pass"+iBin,r.RooArgList(pEPHist,pPdfs[2]))
        pFail   = r.RooAddPdf("tot_fail"+iBin,"tot_fail"+iBin,r.RooArgList(pEFHist))#,pPdfs[3]))
        lPdfs.extend([pPass,pFail])
    lData  = r.RooDataHist("comb_data_obs","comb_data_obs",r.RooArgList(lBase[0]),lCats,lDatas)
    lPdf   = r.RooSimultaneous("tot","tot",lCats)
    for i0 in range(0,len(lPdfs)):
        lPdf .addPdf(lPdfs [i0],lPdfs[i0].GetName().replace("tot_",""))
    lPdf.fitTo(lData,r.RooFit.Extended(),r.RooFit.Minos())
    #for i0 in range(0,lPtBins):
    #    drawPF(lBase[0],[lRDatas[2*i0],lRDatas[2*i0+1]],[[lPdfs[2*i0]],[lPdfs[2*i0+1]]],i0)

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
    iHP.extend([lHP0,lHP1,lHP2])
    iHF.extend([lHF0,lHF1,lHF2])
    masses=[50,75,100,125,150,200,250,300]
    for mass in masses:
        lHP.append(lFile.Get("zqq"+str(mass)+"_pass" +end))
        iHF.append(lFile.Get("zqq"+str(mass)+"_fail" +end))
    iHP.append(lHP3)
    iHF.append(lHF3)
    for lH in (iHP+iHF):
        lH.SetDirectory(0)
    return 
    
if __name__ == "__main__":
    options = parser()
    print options
    lHP =[]
    lHF =[]
    fOutput = options.output 
    fXMin   = options.xmin
    fXMax   = options.xmax
    fNBins  = options.nbins
    load(options.input,lHP,lHF,options.fit2D or (options.cat and not options.fit1D) or (options.card and not options.fit1D))
    if options.fit1D and not options.cat and not options.card:
        if options.passfail:
            fit1DPF(lHP,lHF,options.func)
        else:
            fit1D(lHP,lHF,options.fail,options.func)

    if options.fit2D:
        f2D=True
        if options.passfail:
            fit2DPF(lHP,lHF)
        else:
            fit2D(lHP,lHF,options.fail)

    if options.cat:
        fitCat(lHP,lHF,options.func)
    
    if options.card:
        rhoCard(lHP,lHF)
