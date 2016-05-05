#! /usr/bin/env python
import ROOT as r,sys,math,array
from optparse import OptionParser

fVars=[]
fFuncs=[]

def end():
    if __name__ == '__main__':
        rep = ''
        while not rep in [ 'q', 'Q','a',' ' ]:
            rep = raw_input( 'enter "q" to quit: ' )
            if 1 < len(rep):
                rep = rep[0]

def parser():
    parser = OptionParser()
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
        pFunc.plotOn(iFrame,r.RooFit.LineColor(iColor),r.RooFit.LineStyle(iColor != 50+1),r.RooFit.Components(pFunc.GetName()),r.RooFit.ProjWData(iData))
        iColor+=10

def draw(iVar,iData,iFuncs,iLabel="A"):
    lCan   = r.TCanvas(iLabel,iLabel,800,600)
    lFrame = iVar.frame()
    drawFrame(lFrame,iData,iFuncs)
    lFrame.Draw()
    lCan.Modified()
    lCan.Update()
    end()

def drawPF(iVar,iData,iFuncs,iLabel="A"):
    lCan   = r.TCanvas(iLabel,iLabel,800,600)
    lCan.Divide(2)
    iVar.setBins(20)
    lPFrame = iVar.frame(r.RooFit.Bins(20),r.RooFit.Title("pass"))
    lFFrame = iVar.frame(r.RooFit.Bins(20),r.RooFit.Title("fail"))
    drawFrame(lPFrame,iData[0],iFuncs[0])
    drawFrame(lFFrame,iData[1],iFuncs[1])
    lCan.cd(1)
    lPFrame.Draw()
    lCan.cd(2)
    lFFrame.Draw()
    lCan.Modified()
    lCan.Update()
    end()
    
def workspace(iDatas,iFuncs,iVars,iCat="_cat0"):
    lW = r.RooWorkspace("w"+iCat)
    for var in iVars:
        try:
            var.setConstant(True)
        except:
            pass
    for pData in iDatas:
        getattr(lW,'import')(pData,r.RooFit.RecycleConflictNodes())
    for pFunc in iFuncs:
        getattr(lW,'import')(pFunc,r.RooFit.RecycleConflictNodes())
    lW.writeToFile("simple-shapes-RooDataHist.root")

def writeHist(iH,iShift,iTemps):  
    lOut = r.TFile("simple-shapes-Hist.root","RECREATE")
    for pHist in iH:
      pHist.Write()
    iShift.setVal(5.)
    for pTemp in iTemps:
      lUp = pTemp.createHistogram("x")
      lUp.SetTitle(pTemp.GetName()+"_scaleUp")
      lUp.SetName (pTemp.GetName()+"_scaleUp")
      lUp.Write()
  
    iShift.setVal(-5.)
    for pTemp in iTemps:
      lDown = pTemp.createHistogram("x")
      lDown.SetTitle(pTemp.GetName()+"_scaleDown")
      lDown.SetName (pTemp.GetName()+"_scaleDown")
      lDown.Write()
    
def baseVars(i1D=False):
    lMSD    = r.RooRealVar("x","x",40,200)
    lPt  = r.RooRealVar   ("pt","pt",500,1500)
    lRho = r.RooFormulaVar("rho","log(x*x/pt)",r.RooArgList(lMSD,lPt))
    lEff    = r.RooRealVar("veff"      ,"veff"      ,0.5 ,0.,1.0)
    lEffQCD = r.RooRealVar("qcdeff"     ,"qcdeff"   ,0.03,0.,0.2)
    lDM     = r.RooRealVar("dm","dm",0.,-10,10)
    lShift  = r.RooFormulaVar("shift",lMSD.GetName()+"-dm",r.RooArgList(lMSD,lDM))  
    lVars=[lMSD,lEff,lEffQCD,lDM,lShift]
    if not i1D:
        lVars.extend([lPt,lRho])
    fVars.extend([lEff,lEffQCD,lDM])
    return lVars

def qcdFunc(iH,iVars,iBin="_cat0",i1D=True):
    lNTot   = r.RooRealVar   ("qcdnorm"+iBin,"qcdnorm"+iBin,(iH[0].Integral()+iH[1].Integral()),0.,3.*(iH[0].Integral()+iH[1].Integral()))
    lNPass  = r.RooFormulaVar("fqpass"+iBin ,"qcdnorm"+iBin+"*(qcdeff)"  ,r.RooArgList(lNTot,iVars[2]))
    lNFail  = r.RooFormulaVar("fqfail"+iBin ,"qcdnorm"+iBin+"*(1-qcdeff)",r.RooArgList(lNTot,iVars[2]))
    lA0     = r.RooRealVar   ("a0"          ,"a0"          ,0.521,0.,1)          
    lA1     = r.RooRealVar   ("a1"          ,"a1"          ,0.00357,-0.01,0.01)
    lA2     = r.RooRealVar   ("a2"          ,"a2"          ,-0.000025,-0.0001,0.0001)
    lSigma1 = r.RooRealVar   ("sigma1"      ,"sigma1"      ,50,0,100);
    lP1     = r.RooRealVar   ("p1" ,"p1", -0.001   ,-0.05  ,0.05)
    #lP2     = r.RooRealVar   ("p2" ,"p2", 0.0   ,-0.001,0.001)
    lR1     = r.RooRealVar   ("r1" ,"r1", 0.0   ,-0.2 ,0.2)
    #lR2     = r.RooRealVar("r2" ,"r2", 0.0   ,-0.001,0.001)
    if i1D:
        lQFuncP = r.RooGenericPdf("qcdpass"+iBin,"(exp(-x*x/sigma1/sigma1)+a0+a1*x+a2*x*x)",r.RooArgList(iVars[0],lSigma1,lA0,lA1,lA2))
        lQFuncF = r.RooGenericPdf("qcdfail"+iBin,"(exp(-x*x/sigma1/sigma1)+a0+a1*x+a2*x*x)",r.RooArgList(iVars[0],lSigma1,lA0,lA1,lA2))
    else:
        lQFuncP = r.RooGenericPdf("qcdpass"+iBin,"(exp(-x*x/sigma1/sigma1)+a0+a1*x+a2*x*x)*(1+p1*(pt-500)+r1*(rho-2.5))",r.RooArgList(iVars[0],lSigma1,lA0,lA1,lA2,lP1,iVars[5],lR1,iVars[6]))
        lQFuncF = r.RooGenericPdf("qcdfail"+iBin,"(exp(-x*x/sigma1/sigma1)+a0+a1*x+a2*x*x)*(1-p1*(pt-500)-r1*(rho-2.5))",r.RooArgList(iVars[0],lSigma1,lA0,lA1,lA2,lP1,iVars[5],lR1,iVars[6]))
    lQCDP   = r.RooExtendPdf ("qcdpassE"+iBin,"qcdpass"+iBin,lQFuncP,lNPass)
    lQCDF   = r.RooExtendPdf ("qcdfailE"+iBin,"qcdfail"+iBin,lQFuncF,lNFail)
    lQCD    = [lQCDP,lQCDF,lQFuncP,lQFuncF]
    fVars.extend([lNTot,lA0,lA1,lA2,lR1,lSigma1,lNPass,lNFail,lP1,lR1])
    fFuncs.extend(lQCD)
    #lA0.setConstant(True)
    #lA1.setConstant(True)
    #lP1.setConstant(True)
    #lR1.setConstant(True)
    #lA2.setConstant(True)
    return lQCD

def histFunc(iH,iVars,iLabel="w",iBin="_cat0",i1D=True):
    lNTot   = r.RooRealVar (iLabel+"norm"+iBin,iLabel+"norm"+iBin,(iH[0].Integral()+iH[1].Integral()),0.,5.*(iH[0].Integral()+iH[1].Integral()))
    lNPass  = r.RooFormulaVar(iLabel+"fpass"+iBin ,iLabel+"norm"+iBin+"*(veff)"  ,r.RooArgList(lNTot,iVars[1]))
    lNFail  = r.RooFormulaVar(iLabel+"fqail"+iBin ,iLabel+"norm"+iBin+"*(1-veff)",r.RooArgList(lNTot,iVars[1]))
    if i1D:
        lPData  = r.RooDataHist(iLabel+"dph" +iBin,iLabel+"dph"+iBin, r.RooArgList(iVars[0]),iH[0])
        lMData  = r.RooDataHist(iLabel+"dmh" +iBin,iLabel+"dmh"+iBin, r.RooArgList(iVars[0]),iH[1]) 
        lP      = r.RooHistPdf (iLabel+"pass"+iBin,iLabel+"pass"+iBin, r.RooArgList(iVars[4]),r.RooArgList(iVars[0]),lPData,1)
        lF      = r.RooHistPdf (iLabel+"fail"+iBin,iLabel+"fail"+iBin, r.RooArgList(iVars[4]),r.RooArgList(iVars[0]),lMData,1)
    else:
        lPData  = r.RooDataHist(iLabel+"dph" +iBin,iLabel+"dph"+iBin, r.RooArgList(iVars[0],iVars[5]),iH[0])
        lMData  = r.RooDataHist(iLabel+"dmh" +iBin,iLabel+"dmh"+iBin, r.RooArgList(iVars[0],iVars[5]),iH[1]) 
        lP      = r.RooHistPdf (iLabel+"pass"+iBin,iLabel+"pass"+iBin, r.RooArgList(iVars[4],iVars[5]),r.RooArgList(iVars[0],iVars[5]),lPData,1)
        lF      = r.RooHistPdf (iLabel+"fail"+iBin,iLabel+"fail"+iBin, r.RooArgList(iVars[4],iVars[5]),r.RooArgList(iVars[0],iVars[5]),lMData,1)
    lEP     = r.RooExtendPdf(iLabel+"pe" +iBin,iLabel+"pe" +iBin,lP,lNPass);
    lEF     = r.RooExtendPdf(iLabel+"fe" +iBin,iLabel+"fe" +iBin,lF,lNFail);
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
    lCats.defineType("pass",1) ;
    lCats.defineType("fail",0) ;
    lBase  = baseVars()
    lPData = r.RooDataHist("data_obs_pass"+iBin,"data_obs_pass"+iBin,r.RooArgList(lBase[0]),iHP[0])
    lFData = r.RooDataHist("data_obs_fail"+iBin,"data_obs_fail"+iBin,r.RooArgList(lBase[0]),iHF[0])
    lData  = r.RooDataHist("comb_data_obs","comb_data_obs",r.RooArgList(lBase[0]),r.RooFit.Index(lCats),r.RooFit.Import("pass",lPData),r.RooFit.Import("fail",lFData))
    lW    = histFunc([iHP[1],iHF[1]],lBase,"wqq",iBin)
    lZ    = histFunc([iHP[2],iHF[2]],lBase,"zqq",iBin)
    lQCD  = qcdFunc ([iHP[3],iHF[3]],lBase,iBin)
    lTotP = r.RooAddPdf("tot_pass"+iBin,"tot_pass"+iBin,r.RooArgList(lQCD[0],lW[2],lZ[2]))
    lTotF = r.RooAddPdf("tot_fail"+iBin,"tot_fail"+iBin,r.RooArgList(lQCD[1],lW[3],lZ[3]))
    lTot  = r.RooSimultaneous("tot","tot",lCats) ;
    lTot.addPdf(lTotP,"pass") ;
    lTot.addPdf(lTotF,"fail") ;    
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
    lCats.defineType("pass",1) ;
    lCats.defineType("fail",0) ;
    lPData = r.RooDataHist("data_obs_pass"+iBin,"data_obs_pass"+iBin,r.RooArgList(lBase[0],lBase[5]),iHP[0])
    lFData = r.RooDataHist("data_obs_fail"+iBin,"data_obs_fail"+iBin,r.RooArgList(lBase[0],lBase[5]),iHF[0])
    lData  = r.RooDataHist("comb_data_obs","comb_data_obs",r.RooArgList(lBase[0],lBase[5]),r.RooFit.Index(lCats),r.RooFit.Import("pass",lPData),r.RooFit.Import("fail",lFData))
    lW    = histFunc([iHP[1],iHF[1]],lBase,"wqq",iBin,False)
    lZ    = histFunc([iHP[2],iHF[2]],lBase,"zqq",iBin,False)
    lQCD  = qcdFunc ([iHP[3],iHF[3]],lBase,iBin,False)
    lTotP = r.RooAddPdf("tot_pass"+iBin,"tot_pass"+iBin,r.RooArgList(lQCD[0],lW[2],lZ[2]))
    lTotF = r.RooAddPdf("tot_fail"+iBin,"tot_fail"+iBin,r.RooArgList(lQCD[1],lW[3],lZ[3]))
    lTot  = r.RooSimultaneous("tot","tot",lCats) ;
    lTot.addPdf(lTotP,"pass") ;
    lTot.addPdf(lTotF,"fail") ;    
    lTot.fitTo(lData,r.RooFit.Extended(),r.RooFit.ConditionalObservables(r.RooArgSet(lBase[5])))
    workspace([lPData,lFData],fFuncs,fVars,iBin)
    lPFuncs=[lTotP,lQCD[0]]
    lFFuncs=[lTotF,lQCD[1]]
    drawPF(lBase[0],[lPData,lFData],[lPFuncs,lFFuncs])

#def rhoCard():
#pt categories with rhalpha Fit constraints does on as RateParams

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
    load(options.input,lHP,lHF,options.fit2D)
    if options.fit1D:
        if options.passfail:
            fit1DPF(lHP,lHF)
        else:
            fit1D(lHP,lHF,options.fail)

    if options.fit2D:
        if options.passfail:
            fit2DPF(lHP,lHF)
        else:
            fit2D(lHP,lHF,options.fail)
