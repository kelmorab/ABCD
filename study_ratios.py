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
            dataMode="data", verbosity=1):
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
    print(nBinsX, nBinsY)
    histo2D_significance=ROOT.TH2D("histo2D_significance","histo2D_significance",nBinsX, scanMinX, 1.0, nBinsY, scanMinY, 1.0 )
    histo2D_significance.Sumw2()
    histo2D_closure=ROOT.TH2D("histo2D_closure","histo2D_closure",nBinsX, scanMinX, 1.0, nBinsY, scanMinY, 1.0 )
    histo2D_closure.Sumw2()
    histo2D_closure_error=ROOT.TH2D("histo2D_closure_error","histo2D_closure_error",nBinsX, scanMinX, 1.0, nBinsY, scanMinY, 1.0 )
    histo2D_closure_error.Sumw2()
    histo2D_maxAllowedNonClosure=ROOT.TH2D("histo2D_closure_error","histo2D_closure_error",nBinsX, scanMinX, 1.0, nBinsY, scanMinY, 1.0 )
    histo2D_maxAllowedNonClosure.Sumw2()

    # do the scan
    for xCut in numpy.arange(scanMinX, scanMaxX, 0.01):
        for yCut in numpy.arange(scanMinY, scanMaxY, 0.01):
            if verbosity>0:
                print(" ")
                print("Scan at x,y: "+str(xCut)+","+str(yCut))
            xCut_upper=xCut
            yCut_upper=yCut
            # if exludeBlob:
            #     if xCut>0.7:
            #         xCut_upper=0.7
            #     if yCut>0.7:
            #         yCut_upper=0.7
            theABCD.setWindow("A", xCut, windowMaxX, yCut, windowMaxY)
            theABCD.setWindow("B", xCut, windowMaxX, windowMinY, yCut_upper-0.0001)
            theABCD.setWindow("C", windowMinX, xCut_upper-0.0001, yCut, windowMaxY)
            theABCD.setWindow("D", windowMinX, xCut_upper-0.0001, windowMinY, yCut_upper-0.0001)
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
                maxAllowedNonClosure=0
                nonclosure1, nonclosure2 = theABCD.getMaxAllowedNonClosures(sig=1, dataMode=dataMode, verbosity=verbosity)
                maxAllowedNonClosure=max(nonclosure1, nonclosure2)
                if verbosity>0:
                    print("closure: "+str(round(closure,3))+" +-"+str(round(closure_error,3)))
                    print("closure_ALT: "+str(round(closure_alt,3))+" +-"+str(round(closure_alt_error,3)))
                    print("significance: "+str(round(significance,3)))
                    print("maxAllowedNonClosure; "+str(round(maxAllowedNonClosure,3)))
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
                maxAllowedNonClosure=0
                nonclosure1, nonclosure2 = theABCD.getMaxAllowedNonClosures(sig=1, dataMode=dataMode, verbosity=verbosity)
                maxAllowedNonClosure=max(nonclosure1, nonclosure2)
                if verbosity>0:
                    print("closure: "+str(round(closure,3))+" +-"+str(round(closure_error,3)))
                    print("closure: "+str(round(closure_alt,3))+" +-"+str(round(closure_alt_error,3)))
                    print("significance:"+str(round(significance,3)))
                    print("maxAllowedNonClosure; "+str(round(maxAllowedNonClosure,3)))
            # Fill histos
            # First handle case that x,y lie outside of general window range
            if xCut<windowMinX or yCut<windowMinY:
                significance=0
                closure=0
                closure_error=0
                maxAllowedNonClosure=0
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
            histo2D_maxAllowedNonClosure.SetBinContent(binX, binY, maxAllowedNonClosure)
    return histo2D_significance, histo2D_closure, histo2D_closure_error, histo2D_maxAllowedNonClosure

            

###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

# Scan closure = predicted/observed for different lower_bounds, upper bound and ABCD cuts


# create toplevel dir
# toplevelDirName="closures_allRegions_jetComb12_V6"
toplevelDirName="closures_allRegions_jetComb02_V6"
if not os.path.isdir(toplevelDirName):
    os.mkdir(toplevelDirName)

# masses=["15","40","55"]
# lifetimes=["0","0p05","0p1","1",]
masses=["15", "40", "55"]
lifetimes=["0","0p05","0p1","1",]
# lifetimes=["0p1",]

histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")

do_VRs=False

rightmargin=0.2
markersize=0.5
ROOT.gStyle.SetPaintTextFormat("3.2f");
for mass in masses:
    for lifetime in lifetimes:
        WP_tagName="M"+mass+"_"+lifetime
        print("")
        print("At "+WP_tagName)
        outputDir_level1=toplevelDirName+"/"+"closures_allRegions_"+WP_tagName+"/"
        if not os.path.isdir(outputDir_level1):
            os.mkdir(outputDir_level1)

        verbosity=0
        # if mass=="15" and lifetime=="1":
        #     verbosity=1

        if mass=="15":
            histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG"
            histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_VR"
            histoString_SR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_invertedChi2"
            histoString_VR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_invertedChi2_VR"
        if mass=="40" or mass=="55":
            histoString_SR="resolved_Higgs_combi_chi2_A_three_phiSorted_02_$SAMPLE_sigprob_$WPTAG"
            histoString_VR="resolved_Higgs_combi_chi2_A_three_phiSorted_02_$SAMPLE_sigprob_$WPTAG_VR"
            histoString_SR_invChi2="resolved_Higgs_combi_chi2_A_three_phiSorted_02_$SAMPLE_sigprob_$WPTAG_invertedChi2"
            histoString_VR_invChi2="resolved_Higgs_combi_chi2_A_three_phiSorted_02_$SAMPLE_sigprob_$WPTAG_invertedChi2_VR"
        histoString_SR=histoString_SR.replace("$WPTAG",WP_tagName)
        histoString_VR=histoString_VR.replace("$WPTAG",WP_tagName)
        histoString_SR_invChi2=histoString_SR_invChi2.replace("$WPTAG",WP_tagName)
        histoString_VR_invChi2=histoString_VR_invChi2.replace("$WPTAG",WP_tagName)
        
        lowerbound_list=[0.0, 0.3, 0.4, 0.5, 0.6,0.7, 0.8]
        # lowerbound_list=[0.7, 0.8]
        for lowerbound in lowerbound_list:
            outputDir_level2=outputDir_level1+"/"+"closures_"+"_"+WP_tagName+"_lb_"+str(lowerbound).replace(".","p")+"/"
            if not os.path.isdir(outputDir_level2):
                os.mkdir(outputDir_level2)

            # Do Full range: SR_prime, VR and VR_prime and SR (discard dangerouse stuff)
            windowMinX=lowerbound
            windowMinY=lowerbound
            windowMaxY=1.0
            windowMaxX=1.0
            scanMinX=0.61
            scanMaxX=1.0
            scanMinY=0.61
            scanMaxY=1.0

            region_list=["SR_inv_chi2", "VR", "VR_inv_chi2","SR"]
            histoString_list=[histoString_SR_invChi2, histoString_VR, histoString_VR_invChi2, histoString_SR]
            dataMode_list=["data","data","data","data"]
            for theRegion, theHistoString, theDataMode in zip(region_list, histoString_list, dataMode_list):
                if theRegion=="SR":
                    verbosity=0
                if theRegion=="VR":
                    verbosity=0
                else:
                    verbosity=0        
                histo2D_significance, histo2D_closure, histo2D_closure_error, histo2D_maxAllowedNonClosure = scan2D(windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinY, windowMaxY=windowMaxY,
                                                                                    scanMinX=scanMinX, scanMaxX=scanMaxX, scanMinY=scanMinY, scanMaxY=scanMaxY,
                                                                                    histoFile=histoFile, histoString=theHistoString,
                                                                                    dataMode=theDataMode, verbosity=verbosity)
                histo_list=[histo2D_significance, histo2D_closure, histo2D_closure_error, histo2D_maxAllowedNonClosure]
                fom_list=["significance","closure","closure_error", "max_allowed_nonClosure"]
                # blind
                if theRegion=="SR":
                    histo_list=[histo2D_maxAllowedNonClosure]
                    fom_list=["max_allowed_nonClosure"]
                    
                for theHisto, theFom in zip(histo_list, fom_list):
                    canvas=ROOT.TCanvas("canvas_"+WP_tagName+"_"+theRegion+"_"+theDataMode,"canvas_"+WP_tagName+"_"+theRegion+"_"+theDataMode,1024,768)
                    canvas.SetRightMargin(rightmargin)
                    theHisto.SetTitle(theRegion+" "+theDataMode+" "+WP_tagName+" ["+str(windowMinX)+","+str(windowMaxX)+"]["+str(windowMinY)+","+str(windowMaxY)+"]")
                    theHisto.SetMarkerSize(markersize)
                    theHisto.GetXaxis().SetTitle("x-cut")
                    theHisto.GetYaxis().SetTitle("y-cut")
                    if theFom=="significance":
                        theHisto.GetZaxis().SetTitle("significance = #frac{pred-obs}{total uncertainty}")
                        theHisto.GetZaxis().SetRangeUser(0,5)
                    if theFom=="closure":
                        theHisto.GetZaxis().SetRangeUser(0,5)
                        theHisto.GetZaxis().SetTitle("pred/obs")
                    if theFom=="closure_error":
                        theHisto.GetZaxis().SetRangeUser(0,1)
                        theHisto.GetZaxis().SetTitle("uncertainty(pred/obs)")
                    if theFom=="max_allowed_nonClosure":
                        theHisto.GetZaxis().SetRangeUser(0,1)
                        theHisto.GetZaxis().SetTitle("max allowed nonClosure (1 sigma)")
                        
                    theHisto.Draw("colz text")
                    windowString="["+str(windowMinX)+","+str(windowMaxX)+"]["+str(windowMinY)+","+str(windowMaxY)+"]"
                    windowString=windowString.replace(".","p")
                    canvas.SaveAs(outputDir_level2+theRegion+"_"+theDataMode+"_"+windowString+"_"+"unrestricted_"+theFom+".png")
                    canvas.SaveAs(outputDir_level2+theRegion+"_"+theDataMode+"_"+windowString+"_"+"unrestricted_"+theFom+".pdf")


 
 
            # Do y-restricted range: SR, SR_prime, VR and VR_prime
            windowMinX=lowerbound
            windowMinY=lowerbound
            windowMaxX=0.7
            windowMaxY=1.0
            scanMinX=0.61
            scanMaxX=0.7
            scanMinY=0.61
            scanMaxY=1.0

            region_list=["SR","SR_inv_chi2", "VR", "VR_inv_chi2"]
            histoString_list=[histoString_SR, histoString_SR_invChi2, histoString_VR, histoString_VR_invChi2]
            dataMode_list=["data","data","data", "data"]
            for theRegion, theHistoString, theDataMode in zip(region_list, histoString_list, dataMode_list):
                histo2D_significance, histo2D_closure, histo2D_closure_error, histo2D_maxAllowedNonClosure = scan2D(windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinY, windowMaxY=windowMaxY,
                                                                                    scanMinX=scanMinX, scanMaxX=scanMaxX, scanMinY=scanMinY, scanMaxY=scanMaxY,
                                                                                    histoFile=histoFile, histoString=theHistoString,
                                                                                    dataMode=theDataMode, verbosity=verbosity )
                histo_list=[histo2D_significance, histo2D_closure, histo2D_closure_error, histo2D_maxAllowedNonClosure]
                fom_list=["significance","closure","closure_error", "max_allowed_nonClosure"]
                for theHisto, theFom in zip(histo_list, fom_list):
                    canvas=ROOT.TCanvas("canvas_"+WP_tagName+"_"+theRegion+"_"+theDataMode,"canvas_"+WP_tagName+"_"+theRegion+"_"+theDataMode,1024,768)
                    canvas.SetRightMargin(rightmargin)
                    theHisto.SetTitle(theRegion+" "+theDataMode+" "+WP_tagName+" ["+str(windowMinX)+","+str(windowMaxX)+"]["+str(windowMinY)+","+str(windowMaxY)+"]")
                    theHisto.SetMarkerSize(markersize)
                    theHisto.GetXaxis().SetTitle("x-cut")
                    theHisto.GetYaxis().SetTitle("y-cut")
                    if theFom=="significance":
                        theHisto.GetZaxis().SetTitle("significance = #frac{pred-obs}{total uncertainty}")
                        theHisto.GetZaxis().SetRangeUser(0,5)
                    if theFom=="closure":
                        theHisto.GetZaxis().SetTitle("pred/obs")
                        theHisto.GetZaxis().SetRangeUser(0,5)
                    if theFom=="closure_error":
                        theHisto.GetZaxis().SetTitle("uncertainty(pred/obs)")
                        theHisto.GetZaxis().SetRangeUser(0,1)
                    if theFom=="max_allowed_nonClosure":
                        theHisto.GetZaxis().SetRangeUser(0,1)
                        theHisto.GetZaxis().SetTitle("max allowed nonClosure (1 sigma)")    
                    theHisto.Draw("colz text")
                    windowString="["+str(windowMinX)+","+str(windowMaxX)+"]["+str(windowMinY)+","+str(windowMaxY)+"]"
                    windowString=windowString.replace(".","p")
                    canvas.SaveAs(outputDir_level2+theRegion+"_"+theDataMode+"_"+windowString+"_"+"Xrestricted_"+theFom+".png")
                    canvas.SaveAs(outputDir_level2+theRegion+"_"+theDataMode+"_"+windowString+"_"+"Xrestricted_"+theFom+".pdf")

            # Do y-restricted range: SR, SR_prime, VR and VR_prime
            windowMinX=lowerbound
            windowMinY=lowerbound
            windowMaxX=1.0
            windowMaxY=0.7
            scanMinX=0.61
            scanMaxX=1.0
            scanMinY=0.61
            scanMaxY=0.7

            region_list=["SR","SR_inv_chi2", "VR", "VR_inv_chi2"]
            histoString_list=[histoString_SR, histoString_SR_invChi2, histoString_VR, histoString_VR_invChi2]
            dataMode_list=["data","data","data", "data"]
            for theRegion, theHistoString, theDataMode in zip(region_list, histoString_list, dataMode_list):
                histo2D_significance, histo2D_closure, histo2D_closure_error, histo2D_maxAllowedNonClosure = scan2D(windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinY, windowMaxY=windowMaxY,
                                                                                    scanMinX=scanMinX, scanMaxX=scanMaxX, scanMinY=scanMinY, scanMaxY=scanMaxY,
                                                                                    histoFile=histoFile, histoString=theHistoString,
                                                                                    dataMode=theDataMode, verbosity=0 )
                histo_list=[histo2D_significance, histo2D_closure, histo2D_closure_error, histo2D_maxAllowedNonClosure]
                fom_list=["significance","closure","closure_error", "max_allowed_nonClosure"]
                for theHisto, theFom in zip(histo_list, fom_list):
                    canvas=ROOT.TCanvas("canvas_"+WP_tagName+"_"+theRegion+"_"+theDataMode,"canvas_"+WP_tagName+"_"+theRegion+"_"+theDataMode,1024,768)
                    canvas.SetRightMargin(rightmargin)
                    theHisto.SetTitle(theRegion+" "+theDataMode+" "+WP_tagName+" ["+str(windowMinX)+","+str(windowMaxX)+"]["+str(windowMinY)+","+str(windowMaxY)+"]")
                    theHisto.SetMarkerSize(markersize)
                    theHisto.GetXaxis().SetTitle("x-cut")
                    theHisto.GetYaxis().SetTitle("y-cut")
                    if theFom=="significance":
                        theHisto.GetZaxis().SetTitle("significance = #frac{pred-obs}{total uncertainty}")
                        theHisto.GetZaxis().SetRangeUser(0,5)
                    if theFom=="closure":
                        theHisto.GetZaxis().SetTitle("pred/obs")
                        theHisto.GetZaxis().SetRangeUser(0,5)
                    if theFom=="closure_error":
                        theHisto.GetZaxis().SetTitle("uncertainty(pred/obs)")
                        theHisto.GetZaxis().SetRangeUser(0,1)
                    if theFom=="max_allowed_nonClosure":
                        theHisto.GetZaxis().SetRangeUser(0,1)
                        theHisto.GetZaxis().SetTitle("max allowed nonClosure (1 sigma)")    
                    theHisto.Draw("colz text")
                    windowString="["+str(windowMinX)+","+str(windowMaxX)+"]["+str(windowMinY)+","+str(windowMaxY)+"]"
                    windowString=windowString.replace(".","p")
                    canvas.SaveAs(outputDir_level2+theRegion+"_"+theDataMode+"_"+windowString+"_"+"Yrestricted_"+theFom+".png")
                    canvas.SaveAs(outputDir_level2+theRegion+"_"+theDataMode+"_"+windowString+"_"+"Yrestricted_"+theFom+".pdf")

            # # Do  both restricted range: SR, SR_prime, VR and VR_prime
            # windowMinX=lowerbound
            # windowMinY=lowerbound
            # windowMaxX=0.7
            # windowMaxY=0.7
            # scanMinX=0.61
            # scanMaxX=0.7
            # scanMinY=0.61
            # scanMaxY=0.7

            # region_list=["SR","SR_inv_chi2", "VR", "VR_inv_chi2"]
            # histoString_list=[histoString_SR, histoString_SR_invChi2, histoString_VR, histoString_VR_invChi2]
            # dataMode_list=["data","data","data", "data"]
            # for theRegion, theHistoString, theDataMode in zip(region_list, histoString_list, dataMode_list):
            #     histo2D_significance, histo2D_closure, histo2D_closure_error, histo2D_maxAllowedNonClosure = scan2D(windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinY, windowMaxY=windowMaxY,
            #                                                                         scanMinX=scanMinX, scanMaxX=scanMaxX, scanMinY=scanMinY, scanMaxY=scanMaxY,
            #                                                                         histoFile=histoFile, histoString=theHistoString,
            #                                                                         dataMode=theDataMode, verbosity=0 )
            #     histo_list=[histo2D_significance, histo2D_closure, histo2D_closure_error, histo2D_maxAllowedNonClosure]
            #     fom_list=["significance","closure","closure_error", "max_allowed_nonClosure"]
            #     for theHisto, theFom in zip(histo_list, fom_list):
            #         canvas=ROOT.TCanvas("canvas_"+WP_tagName+"_"+theRegion+"_"+theDataMode,"canvas_"+WP_tagName+"_"+theRegion+"_"+theDataMode,1024,768)
            #         canvas.SetRightMargin(rightmargin)
            #         theHisto.SetTitle(theRegion+" "+theDataMode+" "+WP_tagName+" ["+str(windowMinX)+","+str(windowMaxX)+"]["+str(windowMinY)+","+str(windowMaxY)+"]")
            #         theHisto.SetMarkerSize(markersize)
            #         theHisto.GetXaxis().SetTitle("x-cut")
            #         theHisto.GetYaxis().SetTitle("y-cut")
            #         if theFom=="significance":
            #             theHisto.GetZaxis().SetTitle("significance = #frac{pred-obs}{total uncertainty}")
            #             theHisto.GetZaxis().SetRangeUser(0,5)
            #         if theFom=="closure":
            #             theHisto.GetZaxis().SetTitle("pred/obs")
            #             theHisto.GetZaxis().SetRangeUser(0,5)
            #         if theFom=="closure_error":
            #             theHisto.GetZaxis().SetTitle("uncertainty(pred/obs)")
            #             theHisto.GetZaxis().SetRangeUser(0,1)
            #         if theFom=="max_allowed_nonClosure":
            #             theHisto.GetZaxis().SetRangeUser(0,1)
            #             theHisto.GetZaxis().SetTitle("max allowed nonClosure (1 sigma)")    
            #         theHisto.Draw("colz text")
            #         windowString="["+str(windowMinX)+","+str(windowMaxX)+"]["+str(windowMinY)+","+str(windowMaxY)+"]"
            #         windowString=windowString.replace(".","p")
            #         canvas.SaveAs(outputDir_level2+theRegion+"_"+theDataMode+"_"+windowString+"_"+"BOTHrestricted_"+theFom+".png")
            #         canvas.SaveAs(outputDir_level2+theRegion+"_"+theDataMode+"_"+windowString+"_"+"BOTHrestricted_"+theFom+".pdf")



    #     break
    # break
              


# histoFile.Close()