import ROOT
import math
import ctypes
import numpy
import os
ROOT.gROOT.SetBatch(True)
ROOT.gDirectory.cd('PyROOT:/')
ROOT.gStyle.SetOptStat(0)

from class_ABCD import ABCD, FullABCD

def scanClosureOneDirection(scanDirection="x",
                             windowMinX=0.0, windowMaxX=1., windowMinY=0.0, windowMaxY=1.0, otherDirectionCut=0.0,
                             scanMin=0.5, scanMax=0.9,
                             histoFile="", histoString="",
                             dataMode="data", verbosity=1 ):
    theABCD=None
    if dataMode=="data":
        histo_data_obs=histoFile.Get(histoString.replace("$SAMPLE","data_obs"))
        histo_SM_Higgs = histoFile.Get(histoString.replace("$SAMPLE","SM_Higgs"))
        histo_TTbar = histoFile.Get(histoString.replace("$SAMPLE","TTbar"))
        histo_Others = histoFile.Get(histoString.replace("$SAMPLE","Others"))
        histo_nonQCD = histo_SM_Higgs.Clone()
        histo_nonQCD.Add(histo_TTbar)
        histo_nonQCD.Add(histo_Others)
        theABCD=ABCD()
        theABCD.setDataHisto(histo_data_obs)
        theABCD.setNonQCDHisto(histo_nonQCD)
        theABCD.subtractNonQCDFromData()
    elif dataMode=="simQCD":
        histo_QCD=histoFile.Get(histoString.replace("$SAMPLE","QCD"))
        theABCD=ABCD()
        theABCD.setQCDHisto(histo_QCD)
    else:
        print("?????")
        return None, None

    cutValue_list=[]
    nonClosure_list=[]
    non_Closure_error_list=[]

    for cutValue in numpy.arange(scanMin, scanMax, 0.01):
        if verbosity>0:
            print(cutValue)
        if scanDirection not in ["x","y","diag"]:
            print("DO NOT KNOW THIS DIRECTION")
            return None, None, None
        if scanDirection=="x":
            theABCD.setWindow("A", cutValue, windowMaxX, otherDirectionCut, windowMaxY)
            theABCD.setWindow("B", cutValue, windowMaxX, windowMinY, otherDirectionCut-0.0001)
            theABCD.setWindow("C", windowMinX, cutValue-0.0001, otherDirectionCut, windowMaxY)
            theABCD.setWindow("D", windowMinX, cutValue-0.0001, windowMinY, otherDirectionCut-0.0001)
        if scanDirection=="y":
            theABCD.setWindow("A", otherDirectionCut, windowMaxX, cutValue, windowMaxY)
            theABCD.setWindow("B", otherDirectionCut, windowMaxX, windowMinY, cutValue-0.0001)
            theABCD.setWindow("C", windowMinX, otherDirectionCut-0.0001, cutValue, windowMaxY)
            theABCD.setWindow("D", windowMinX, otherDirectionCut-0.0001, windowMinY, cutValue-0.0001)
        if scanDirection=="diag":
            theABCD.setWindow("A", cutValue, windowMaxX, cutValue, windowMaxY)
            theABCD.setWindow("B", cutValue, windowMaxX, windowMinY, cutValue-0.0001)
            theABCD.setWindow("C", windowMinX, cutValue-0.0001, cutValue, windowMaxY)
            theABCD.setWindow("D", windowMinX, cutValue-0.0001, windowMinY, cutValue-0.0001)    
        if verbosity>0:
            theABCD.printWindow("A")
            theABCD.printWindow("B")
            theABCD.printWindow("C")
            theABCD.printWindow("D")
            print("Integral of whole histo: "+str(theABCD.getIntegralAndErrorInWholeHisto()[0]))
            print("Sum of window Integrals: "+str(theABCD.getIntegralAndErrorInRegion("A")[0]+theABCD.getIntegralAndErrorInRegion("B")[0]+theABCD.getIntegralAndErrorInRegion("C")[0]+theABCD.getIntegralAndErrorInRegion("D")[0]))
            print("")
        if dataMode=="data":
            if verbosity>0:
                print("ABCD in these regions")
                obs_A, obs_error_A = theABCD.getIntegralAndErrorInRegion("A", "data")
                print("observed A:", obs_A, obs_error_A)
                pred_A, pred_error_A = theABCD.predictRegionA()
                print("predicted QCD A: ",pred_A, pred_error_A )
                nonQCD_A, nonQCD_error_A = theABCD.getIntegralAndErrorInRegion("A", "nonQCD")
                print("nonQCD A: ",nonQCD_A, nonQCD_error_A )
                print("-> prediction in A:", pred_A+nonQCD_A  )
            closure, closure_error=theABCD.getNonClosureWrtData()
            closure_alt, closure_alt_error=theABCD.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor=1., correctionFactor_unc=0.0, add_unc=0.0, )
            if verbosity>0:
                print("closure: "+str(round(closure,3))+" +-"+str(round(closure_error,3)))
                print("closure_ALT: "+str(round(closure_alt,3))+" +-"+str(round(closure_alt_error,3)))
            cutValue_list.append(cutValue)
            nonClosure_list.append(closure)
            non_Closure_error_list.append(closure_error)
        elif dataMode=="simQCD":
            if verbosity>0:
                print("Integral of whole histo: "+str(theABCD.getIntegralAndErrorInWholeHisto()[0]))
                print("Sum of window Integrals: "+str(theABCD.getIntegralAndErrorInRegion("A")[0]+theABCD.getIntegralAndErrorInRegion("B")[0]+theABCD.getIntegralAndErrorInRegion("C")[0]+theABCD.getIntegralAndErrorInRegion("D")[0]))
                print("")
                print("ABCD in these regions")
                obs_A, obs_error_A = theABCD.getIntegralAndErrorInRegion("A", "QCD")
                print("observed QCD A:", obs_A, obs_error_A)
                pred_A, pred_error_A = theABCD.predictRegionA()
                print("predicted QCD A: ",pred_A, pred_error_A )
            closure, closure_error=theABCD.getNonClosure()
            closure_alt, closure_alt_error=theABCD.getNonClosureWrtSimQCDWithCorrectionAndAddUncertainty(correctionFactor=1., correctionFactor_unc=0.0, add_unc=0.0, )
            if verbosity>0:
                print("closure: "+str(round(closure,3))+" +-"+str(round(closure_error,3)))
                print("closure: "+str(round(closure_alt,3))+" +-"+str(round(closure_alt_error,3)))
            cutValue_list.append(cutValue)
            nonClosure_list.append(closure)  
            non_Closure_error_list.append(closure_error)          

    return cutValue_list, nonClosure_list, non_Closure_error_list


def scan2D(windowMinX=0.0, windowMaxX=1., windowMinY=0.0, windowMaxY=1.0,
            scanMinX=0.6, scanMaxX=1.0, scanMinY=0.6, scanMaxY=1.0,
            histoFile="", histoString="",
            dataMode="data", verbosity=1 ):
    stepsize=0.01        

    theABCD=None
    if dataMode=="data":
        histo_data_obs=histoFile.Get(histoString.replace("$SAMPLE","data_obs"))
        histo_SM_Higgs = histoFile.Get(histoString.replace("$SAMPLE","SM_Higgs"))
        histo_TTbar = histoFile.Get(histoString.replace("$SAMPLE","TTbar"))
        histo_Others = histoFile.Get(histoString.replace("$SAMPLE","Others"))
        histo_nonQCD = histo_SM_Higgs.Clone()
        histo_nonQCD.Add(histo_TTbar)
        histo_nonQCD.Add(histo_Others)
        theABCD=ABCD()
        theABCD.setDataHisto(histo_data_obs)
        theABCD.setNonQCDHisto(histo_nonQCD)
        theABCD.subtractNonQCDFromData()
    elif dataMode=="simQCD":
        histo_QCD=histoFile.Get(histoString.replace("$SAMPLE","QCD"))
        theABCD=ABCD()
        theABCD.setQCDHisto(histo_QCD)
    else:
        print("?????")
        return None, None

    # set up histos
    nBinsX=int((1.0-scanMinX)/float(stepsize))
    nBinsY=int((1.0-scanMinY)/float(stepsize))
    histo2D_significance=ROOT.TH2D("histo2D_significance","histo2D_significance",nBinsX, scanMinX, 1.0, nBinsY, scanMinY, 1.0 )
    histo2D_significance.Sumw2()
    histo2D_closure=ROOT.TH2D("histo2D_closure","histo2D_closure",nBinsX, scanMinX, 1.0, nBinsY, scanMinY, 1.0 )
    histo2D_closure.Sumw2()
    histo2D_closure_error=ROOT.TH2D("histo2D_closure_error","histo2D_closure_error",nBinsX, scanMinX, 1.0, nBinsY, scanMinY, 1.0 )
    histo2D_closure_error.Sumw2()

    # do the scan
    for xCut in numpy.arange(scanMinX, scanMaxX, 0.01):
        for yCut in numpy.arange(scanMinY, scanMaxY, 0.01):
            if verbosity>0:
                print("Scan at x,y: "+str(xCut)+","+str(yCut))
            theABCD.setWindow("A", xCut, windowMaxX, yCut, windowMaxY)
            theABCD.setWindow("B", xCut, windowMaxX, windowMinY, yCut-0.0001)
            theABCD.setWindow("C", windowMinX, xCut-0.0001, yCut, windowMaxY)
            theABCD.setWindow("D", windowMinX, xCut-0.0001, windowMinY, yCut-0.0001)
            if verbosity>0:
                theABCD.printWindow("A")
                theABCD.printWindow("B")
                theABCD.printWindow("C")
                theABCD.printWindow("D")
                print("Integral of whole histo: "+str(theABCD.getIntegralAndErrorInWholeHisto()[0]))
                print("Sum of window Integrals: "+str(theABCD.getIntegralAndErrorInRegion("A")[0]+theABCD.getIntegralAndErrorInRegion("B")[0]+theABCD.getIntegralAndErrorInRegion("C")[0]+theABCD.getIntegralAndErrorInRegion("D")[0]))
                print("")
            # calculate stuff
            closure, closure_error, significance = 0,0,0
            if dataMode=="data":
                if verbosity>0:
                    print("ABCD in these regions")
                    obs_A, obs_error_A = theABCD.getIntegralAndErrorInRegion("A", "data")
                    print("observed A:", obs_A, obs_error_A)
                    pred_A, pred_error_A = theABCD.predictRegionA()
                    print("predicted QCD A: ",pred_A, pred_error_A )
                    nonQCD_A, nonQCD_error_A = theABCD.getIntegralAndErrorInRegion("A", "nonQCD")
                    print("nonQCD A: ",nonQCD_A, nonQCD_error_A )
                    print("-> prediction in A:", pred_A+nonQCD_A  )
                closure, closure_error=theABCD.getNonClosureWrtData()
                closure_alt, closure_alt_error=theABCD.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor=1., correctionFactor_unc=0.0, add_unc=0.0, )
                significance=theABCD.getSimpleClosureSignificance(dataMode=dataMode, verbosity=verbosity)
                if verbosity>0:
                    print("closure: "+str(round(closure,3))+" +-"+str(round(closure_error,3)))
                    print("closure_ALT: "+str(round(closure_alt,3))+" +-"+str(round(closure_alt_error,3)))
                    print("significance:"+str(round(significance,3)))
            elif dataMode=="simQCD":
                if verbosity>0:
                    print("Integral of whole histo: "+str(theABCD.getIntegralAndErrorInWholeHisto()[0]))
                    print("Sum of window Integrals: "+str(theABCD.getIntegralAndErrorInRegion("A")[0]+theABCD.getIntegralAndErrorInRegion("B")[0]+theABCD.getIntegralAndErrorInRegion("C")[0]+theABCD.getIntegralAndErrorInRegion("D")[0]))
                    print("")
                    print("ABCD in these regions")
                    obs_A, obs_error_A = theABCD.getIntegralAndErrorInRegion("A", "QCD")
                    print("observed QCD A:", obs_A, obs_error_A)
                    pred_A, pred_error_A = theABCD.predictRegionA()
                    print("predicted QCD A: ",pred_A, pred_error_A )
                closure, closure_error=theABCD.getNonClosure()
                closure_alt, closure_alt_error=theABCD.getNonClosureWrtSimQCDWithCorrectionAndAddUncertainty(correctionFactor=1., correctionFactor_unc=0.0, add_unc=0.0, )
                significance=theABCD.getSimpleClosureSignificance(dataMode=dataMode, verbosity=verbosity)
                if verbosity>0:
                    print("closure: "+str(round(closure,3))+" +-"+str(round(closure_error,3)))
                    print("closure: "+str(round(closure_alt,3))+" +-"+str(round(closure_alt_error,3)))
                    print("significance:"+str(round(significance,3)))
            # Fill histos
            binX=histo2D_significance.GetXaxis().FindBin(xCut)
            binY=histo2D_significance.GetYaxis().FindBin(yCut)
            histo2D_significance.SetBinContent(binX, binY, significance)
            # histo2D_significance.GetXaxis().SetTitle("x-Cut")
            # histo2D_significance.GetYaxis().SetTitle("y-Cut")
            # histo2D_significance.GetZaxis().SetTitle("significance = #frac{pred-obs}{total uncertainty}")            
            histo2D_closure.SetBinContent(binX, binY, closure)
            histo2D_closure.SetBinError(binX, binY, closure_error)
            # histo2D_significance.GetXaxis().SetTitle("x-Cut")
            # histo2D_significance.GetYaxis().SetTitle("y-Cut")
            # histo2D_significance.GetZaxis().SetTitle("pred/obs")
            histo2D_closure_error.SetBinContent(binX, binY, closure_error)
            # histo2D_significance.GetXaxis().SetTitle("x-Cut")
            # histo2D_significance.GetYaxis().SetTitle("y-Cut")
            # histo2D_significance.GetZaxis().SetTitle("uncertainty of pred/obs")
    return histo2D_significance, histo2D_closure, histo2D_closure_error

            

###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

# Scan closure = predicted/observed for different lower_bounds, upper bound and ABCD cuts


# create toplevel dir
toplevelDirName="closures_allRegions"
if not os.path.isdir(toplevelDirName):
    os.mkdir(toplevelDirName)

masses=["15","40","55"]
lifetimes=["0","0p05","0p1","1",]

histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
WP_tagName="M"+"15"+"_"+"1"
histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_VR"
histoString_VR=histoString_VR.replace("$WPTAG",WP_tagName)

histoString=histoString_VR
theABCD=None
dataMode="data"
if dataMode=="data":
    histo_data_obs=histoFile.Get(histoString.replace("$SAMPLE","data_obs"))
    histo_SM_Higgs = histoFile.Get(histoString.replace("$SAMPLE","SM_Higgs"))
    histo_TTbar = histoFile.Get(histoString.replace("$SAMPLE","TTbar"))
    histo_Others = histoFile.Get(histoString.replace("$SAMPLE","Others"))
    histo_nonQCD = histo_SM_Higgs.Clone()
    histo_nonQCD.Add(histo_TTbar)
    histo_nonQCD.Add(histo_Others)
    theABCD=ABCD()
    theABCD.setDataHisto(histo_data_obs)
    theABCD.setNonQCDHisto(histo_nonQCD)
    theABCD.subtractNonQCDFromData()
elif dataMode=="simQCD":
    histo_QCD=histoFile.Get(histoString.replace("$SAMPLE","QCD"))
    theABCD=ABCD()
    theABCD.setQCDHisto(histo_QCD)

xCut=0.91
yCut=0.91
windowMinX=0.3
windowMinY=0.3
windowMaxX=1.0
windowMaxY=1.0

theABCD.setWindow("A", xCut, windowMaxX, yCut, windowMaxY)
theABCD.setWindow("B", xCut, windowMaxX, windowMinY, yCut-0.0001)
theABCD.setWindow("C", windowMinX, xCut-0.0001, yCut, windowMaxY)
theABCD.setWindow("D", windowMinX, xCut-0.0001, windowMinY, yCut-0.0001)

print("ABCD in these regions")
obs_A, obs_error_A = theABCD.getIntegralAndErrorInRegion("A", "data")
print("observed A:", obs_A, obs_error_A)
pred_A, pred_error_A = theABCD.predictRegionA()
print("predicted QCD A: ",pred_A, pred_error_A )
nonQCD_A, nonQCD_error_A = theABCD.getIntegralAndErrorInRegion("A", "nonQCD")
print("nonQCD A: ",nonQCD_A, nonQCD_error_A )
print("-> prediction in A:", pred_A+nonQCD_A  )
closure, closure_error=theABCD.getNonClosureWrtData()
closure_alt, closure_alt_error=theABCD.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor=1., correctionFactor_unc=0.0, add_unc=0.0, )
significance=theABCD.getSimpleClosureSignificance(dataMode="data", verbosity=1)
print("closure: "+str(round(closure,3))+" +-"+str(round(closure_error,3)))
print("closure_ALT: "+str(round(closure_alt,3))+" +-"+str(round(closure_alt_error,3)))
print("significance:"+str(round(significance,3)))