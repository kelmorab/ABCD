import ROOT
import math
import ctypes
import numpy
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
from class_ABCD import ABCD

def scanClosure2D(windowMinX=0.0, windowMaxX=1., windowMinY=0.0, windowMaxY=1.0,
                    scanMin=0.7, scanMax=0.98,
                    histoFile="", histoString="",
                    dataMode="data", verbosity=0 ):

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

    stepX=0.01
    stepY=stepX
    nBins=int((scanMax+0.5*stepX-scanMin-0.5*stepX)/stepX)
    closureHisto=ROOT.TH2D("h_closure_"+dataMode, "h_closure_"+dataMode, nBins,scanMin-0.5*stepX,scanMax-0.5*stepX,nBins,scanMin-0.5*stepY, scanMax-0.5*stepY)
    closuresList=[]
    for cutX in numpy.arange(scanMin, scanMax,stepX):
        for cutY in numpy.arange(scanMin, scanMax, stepY):
            if verbosity>0:
                print("-----------------------------")
            theABCD.setWindow("A", cutX, windowMaxX, cutY, windowMaxY)
            theABCD.setWindow("B", cutX, windowMaxX, windowMinY, cutY-0.0001)
            theABCD.setWindow("C", windowMinX, cutX-0.0001, cutY, windowMaxY)
            theABCD.setWindow("D", windowMinX, cutX-0.0001, windowMinY, cutY-0.0001)
            if verbosity>0:
                theABCD.printWindow("A")
                theABCD.printWindow("B")
                theABCD.printWindow("C")
                theABCD.printWindow("D")
                print("Integral of whole histo: "+str(theABCD.getIntegralAndErrorInWholeHisto()[0]))
                print("Sum of window Integrals: "+str(theABCD.getIntegralAndErrorInRegion("A")[0]+theABCD.getIntegralAndErrorInRegion("B")[0]+theABCD.getIntegralAndErrorInRegion("C")[0]+theABCD.getIntegralAndErrorInRegion("D")[0]))
                print("")

            closure=0
            closure_error=0
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
                if verbosity>0: print("non-closure: "+str(round(closure,3))+" +- "+str(round(closure_error,3)))

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
                if verbosity>0: print("non-closure: "+str(round(closure,3))+" +- "+str(round(closure_error,3)))

            binX=closureHisto.GetXaxis().FindBin(cutX)
            binY=closureHisto.GetYaxis().FindBin(cutY)
            closureHisto.SetBinContent(binX, binY, round(closure,2))
            closureHisto.SetBinError(binX, binY, round(closure_error,2))
            
            closuresList.append([round(cutX,2), round(cutY,2),round(closure,2), round(closure_error,2)])
    return closureHisto, closuresList



##############################################################################################################################################################################
histoString_SR_nom="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_0p1"
histoString_VR_nom="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_0p1_VR"
histoString_SR_chi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_0p1_invertedChi2"
histoString_VR_chi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_0p1_invertedChi2_VR"

outputDir="closureScanPlots_M15_0p1/"

for bmPoint in [0.7, 0.8, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91]:

    histoString_SR=histoString_SR_nom
    histoString_VR=histoString_VR_nom
    histoFile=ROOT.TFile("all_histos_datacards_main_v6.root","r")

    verb=0
    allLL=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    global_max=0

    bmGraph_data_VR_0p89_0p89=ROOT.TGraphErrors()
    bmGraph_data_VR_0p89_0p89.SetLineColor(ROOT.kMagenta+3)
    bmGraph_data_VR_0p89_0p89.SetTitle("data VR ("+str(bmPoint)+","+str(bmPoint)+")")
    bmClosureList=[]
    bmClosureErrorList=[]
    bmLowerLimitList=[]
    bmGraphXOffset=0.0
    for lowerLimit in allLL:
        histoString=histoString_VR
        closureHisto, dum = scanClosure2D(
                                    windowMinX=lowerLimit, windowMaxX=1., windowMinY=lowerLimit, windowMaxY=1.0,
                                    scanMin=0.7, scanMax=0.96,
                                    histoFile=histoFile, histoString=histoString,
                                    dataMode="data", verbosity=verb )
        canvas=ROOT.TCanvas("c_closure_data_VR","c_closure_data_VR",1024,768)
        ROOT.gStyle.SetOptStat(0)
        closureHisto.SetTitle("closure data VR, ("+str(lowerLimit)+", 1.0) X ("+str(lowerLimit)+", 1.0)")
        closureHisto.SetMarkerSize(0.6)
        closureHisto.GetXaxis().SetTitle("AB vs CD cut")
        closureHisto.GetYaxis().SetTitle("AC vs BD cut")
        closureHisto.Draw("colz text")
        # closureHisto.Draw("text same")
        lowerLimitString=str(lowerLimit).replace(".","p")
        canvas.SaveAs(outputDir+"2D_closure_data_VR_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.png")
        canvas.SaveAs(outputDir+"2D_closure_data_VR_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.pdf")
        bmClosureBin=closureHisto.GetXaxis().FindBin(bmPoint)
        bmClosure=closureHisto.GetBinContent(bmClosureBin, bmClosureBin)
        bmClosure_error=closureHisto.GetBinError(bmClosureBin, bmClosureBin)
        bmClosureList.append(bmClosure)
        bmClosureErrorList.append(bmClosure_error)
        bmLowerLimitList.append(lowerLimit)
    for i in range(len(bmClosureList)):
        bmGraph_data_VR_0p89_0p89.SetPoint(i,bmLowerLimitList[i]+bmGraphXOffset, bmClosureList[i])
        bmGraph_data_VR_0p89_0p89.SetPointError(i,0, bmClosureErrorList[i] )
        global_max=max(global_max, max(bmClosureList))

    bmGraph_QCD_VR_0p89_0p89=ROOT.TGraphErrors()
    bmGraph_QCD_VR_0p89_0p89.SetLineColor(ROOT.kMagenta+3)
    bmGraph_QCD_VR_0p89_0p89.SetLineStyle(2)
    bmGraph_QCD_VR_0p89_0p89.SetTitle("QCD VR ("+str(bmPoint)+","+str(bmPoint)+")")
    bmClosureList=[]
    bmClosureErrorList=[]
    bmLowerLimitList=[]
    bmGraphXOffset=0.005
    for lowerLimit in allLL:
        histoString=histoString_VR
        closureHisto, dum = scanClosure2D(
                                    windowMinX=lowerLimit, windowMaxX=1., windowMinY=lowerLimit, windowMaxY=1.0,
                                    scanMin=0.7, scanMax=0.96,
                                    histoFile=histoFile, histoString=histoString,
                                    dataMode="simQCD", verbosity=verb )
        canvas=ROOT.TCanvas("c_closure_QCD_VR","c_closure_QCD_VR",1024,768)
        ROOT.gStyle.SetOptStat(0)
        closureHisto.SetTitle("closure QCD VR, ("+str(lowerLimit)+", 1.0) X ("+str(lowerLimit)+", 1.0)")
        closureHisto.SetMarkerSize(0.6)
        closureHisto.GetXaxis().SetTitle("AB vs CD cut")
        closureHisto.GetYaxis().SetTitle("AC vs BD cut")
        closureHisto.Draw("colz text")
        # closureHisto.Draw("text same")
        lowerLimitString=str(lowerLimit).replace(".","p")
        canvas.SaveAs(outputDir+"2D_closure_QCD_VR_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.png")
        canvas.SaveAs(outputDir+"2D_closure_QCD_VR_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.pdf")
        bmClosureBin=closureHisto.GetXaxis().FindBin(bmPoint)
        bmClosure=closureHisto.GetBinContent(bmClosureBin, bmClosureBin)
        bmClosure_error=closureHisto.GetBinError(bmClosureBin, bmClosureBin)
        bmClosureErrorList.append(bmClosure_error)
        bmClosureList.append(bmClosure)
        bmLowerLimitList.append(lowerLimit)
    for i in range(len(bmClosureList)):
        bmGraph_QCD_VR_0p89_0p89.SetPoint(i,bmLowerLimitList[i]+bmGraphXOffset, bmClosureList[i])
        bmGraph_QCD_VR_0p89_0p89.SetPointError(i,0, bmClosureErrorList[i] )
        global_max=max(global_max, max(bmClosureList))

    bmGraph_QCD_SR_0p89_0p89=ROOT.TGraphErrors()
    bmGraph_QCD_SR_0p89_0p89.SetLineColor(ROOT.kRed)
    bmGraph_QCD_SR_0p89_0p89.SetLineStyle(2)
    bmGraph_QCD_SR_0p89_0p89.SetTitle("QCD SR ("+str(bmPoint)+","+str(bmPoint)+")")
    bmClosureList=[]
    bmClosureErrorList=[]
    bmLowerLimitList=[]
    bmGraphXOffset=0.005*2
    for lowerLimit in allLL:
        histoString=histoString_SR
        closureHisto, dum = scanClosure2D(
                                    windowMinX=lowerLimit, windowMaxX=1., windowMinY=lowerLimit, windowMaxY=1.0,
                                    scanMin=0.7, scanMax=0.96,
                                    histoFile=histoFile, histoString=histoString,
                                    dataMode="simQCD", verbosity=verb )
        canvas=ROOT.TCanvas("c_closure_QCD_SR","c_closure_QCD_SR",1024,768)
        ROOT.gStyle.SetOptStat(0)
        closureHisto.SetTitle("closure QCD SR, ("+str(lowerLimit)+", 1.0) X ("+str(lowerLimit)+", 1.0)")
        closureHisto.SetMarkerSize(0.6)
        closureHisto.GetXaxis().SetTitle("AB vs CD cut")
        closureHisto.GetYaxis().SetTitle("AC vs BD cut")
        closureHisto.Draw("colz text")
        
        # closureHisto.Draw("text same")
        lowerLimitString=str(lowerLimit).replace(".","p")
        canvas.SaveAs(outputDir+"2D_closure_QCD_SR_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.png")
        canvas.SaveAs(outputDir+"2D_closure_QCD_SR_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.pdf")
        bmClosureBin=closureHisto.GetXaxis().FindBin(bmPoint)
        bmClosure=closureHisto.GetBinContent(bmClosureBin, bmClosureBin)
        bmClosureList.append(bmClosure)
        bmClosure_error=closureHisto.GetBinError(bmClosureBin, bmClosureBin)
        bmClosureErrorList.append(bmClosure_error)
        bmLowerLimitList.append(lowerLimit)
    for i in range(len(bmClosureList)):
        bmGraph_QCD_SR_0p89_0p89.SetPoint(i,bmLowerLimitList[i]+bmGraphXOffset, bmClosureList[i])
        bmGraph_QCD_SR_0p89_0p89.SetPointError(i,0, bmClosureErrorList[i] )
        global_max=max(global_max, max(bmClosureList))


    #################################
    probeHisto=None

    histoString_SR=histoString_SR_chi2
    histoString_VR=histoString_VR_chi2
    histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")


    bmGraph_data_VR_invertedChi2_0p89_0p89=ROOT.TGraphErrors()
    bmGraph_data_VR_invertedChi2_0p89_0p89.SetLineColor(ROOT.kCyan)
    bmGraph_data_VR_invertedChi2_0p89_0p89.SetTitle("data VR inv Chi2 ("+str(bmPoint)+","+str(bmPoint)+")")
    bmClosureList=[]
    bmClosureErrorList=[]
    bmLowerLimitList=[]
    bmGraphXOffset=0.005*3
    for lowerLimit in allLL:
        histoString=histoString_VR
        closureHisto, dum = scanClosure2D(
                                    windowMinX=lowerLimit, windowMaxX=1., windowMinY=lowerLimit, windowMaxY=1.0,
                                    scanMin=0.7, scanMax=0.96,
                                    histoFile=histoFile, histoString=histoString,
                                    dataMode="data", verbosity=verb )
        canvas=ROOT.TCanvas("c_closure_data_VR_invertedChi2","c_closure_data_VR_invertedChi2",1024,768)
        ROOT.gStyle.SetOptStat(0)
        closureHisto.SetTitle("closure data VR inverted Chi2, ("+str(lowerLimit)+", 1.0) X ("+str(lowerLimit)+", 1.0)")
        closureHisto.SetMarkerSize(0.6)
        closureHisto.GetXaxis().SetTitle("AB vs CD cut")
        closureHisto.GetYaxis().SetTitle("AC vs BD cut")
        closureHisto.Draw("colz text")
        # closureHisto.Draw("text same")
        lowerLimitString=str(lowerLimit).replace(".","p")
        canvas.SaveAs(outputDir+"2D_closure_data_VR_invertedChi2_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.png")
        canvas.SaveAs(outputDir+"2D_closure_data_VR_invertedChi2_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.pdf")
        bmClosureBin=closureHisto.GetXaxis().FindBin(bmPoint)
        bmClosure=closureHisto.GetBinContent(bmClosureBin, bmClosureBin)
        bmClosure_error=closureHisto.GetBinError(bmClosureBin, bmClosureBin)
        bmClosureErrorList.append(bmClosure_error)
        bmClosureList.append(bmClosure)
        bmLowerLimitList.append(lowerLimit)
    for i in range(len(bmClosureList)):
        bmGraph_data_VR_invertedChi2_0p89_0p89.SetPoint(i,bmLowerLimitList[i], bmClosureList[i])
        bmGraph_data_VR_invertedChi2_0p89_0p89.SetPointError(i,0, bmClosureErrorList[i])
        
        global_max=max(global_max, max(bmClosureList))

    bmGraph_QCD_VR_invertedChi2_0p89_0p89=ROOT.TGraphErrors()
    bmGraph_QCD_VR_invertedChi2_0p89_0p89.SetLineColor(ROOT.kCyan)
    bmGraph_QCD_VR_invertedChi2_0p89_0p89.SetLineStyle(2)
    bmGraph_QCD_VR_invertedChi2_0p89_0p89.SetTitle("QCD VR inv Chi2 ("+str(bmPoint)+","+str(bmPoint)+")")
    bmClosureList=[]
    bmClosureErrorList=[]
    bmLowerLimitList=[]
    for lowerLimit in allLL:
        histoString=histoString_VR
        closureHisto, dum = scanClosure2D(
                                    windowMinX=lowerLimit, windowMaxX=1., windowMinY=lowerLimit, windowMaxY=1.0,
                                    scanMin=0.7, scanMax=0.96,
                                    histoFile=histoFile, histoString=histoString,
                                    dataMode="simQCD", verbosity=verb )
        canvas=ROOT.TCanvas("c_closure_QCD_VR_invertedChi2","c_closure_QCD_VR_invertedChi2",1024,768)
        ROOT.gStyle.SetOptStat(0)
        closureHisto.SetTitle("closure QCD VR inverted Chi2, ("+str(lowerLimit)+", 1.0) X ("+str(lowerLimit)+", 1.0)")
        closureHisto.SetMarkerSize(0.6)
        closureHisto.GetXaxis().SetTitle("AB vs CD cut")
        closureHisto.GetYaxis().SetTitle("AC vs BD cut")
        closureHisto.Draw("colz text")
        probeHisto=closureHisto
        # closureHisto.Draw("text same")
        lowerLimitString=str(lowerLimit).replace(".","p")
        canvas.SaveAs(outputDir+"2D_closure_QCD_VR_invertedChi2_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.png")
        canvas.SaveAs(outputDir+"2D_closure_QCD_VR_invertedChi2_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.pdf")
        bmClosureBin=closureHisto.GetXaxis().FindBin(bmPoint)
        bmClosure=closureHisto.GetBinContent(bmClosureBin, bmClosureBin)
        bmClosureList.append(bmClosure)
        bmClosure_error=closureHisto.GetBinError(bmClosureBin, bmClosureBin)
        bmClosureErrorList.append(bmClosure_error)
        bmLowerLimitList.append(lowerLimit)
    for i in range(len(bmClosureList)):
        bmGraph_QCD_VR_invertedChi2_0p89_0p89.SetPoint(i,bmLowerLimitList[i]+bmGraphXOffset, bmClosureList[i])
        bmGraph_QCD_VR_invertedChi2_0p89_0p89.SetPointError(i,0, bmClosureErrorList[i])

    bmGraph_data_SR_invertedChi2_0p89_0p89=ROOT.TGraphErrors()
    bmGraph_data_SR_invertedChi2_0p89_0p89.SetLineColor(ROOT.kBlue)
    bmGraph_data_SR_invertedChi2_0p89_0p89.SetTitle("data SR inv. Chi2 ("+str(bmPoint)+","+str(bmPoint)+")")
    bmClosureList=[]
    bmLowerLimitList=[]
    bmClosureErrorList=[]
    bmGraphXOffset=0.005*4
    for lowerLimit in allLL:
        histoString=histoString_SR
        closureHisto, dum = scanClosure2D(
                                    windowMinX=lowerLimit, windowMaxX=1., windowMinY=lowerLimit, windowMaxY=1.0,
                                    scanMin=0.7, scanMax=0.96,
                                    histoFile=histoFile, histoString=histoString,
                                    dataMode="data", verbosity=verb )
        canvas=ROOT.TCanvas("c_closure_data_SR_invertedChi2","c_closure_data_SR_invertedChi2",1024,768)
        ROOT.gStyle.SetOptStat(0)
        closureHisto.SetTitle("closure data SR inverted Chi2, ("+str(lowerLimit)+", 1.0) X ("+str(lowerLimit)+", 1.0)")
        closureHisto.SetMarkerSize(0.6)
        closureHisto.GetXaxis().SetTitle("AB vs CD cut")
        closureHisto.GetYaxis().SetTitle("AC vs BD cut")
        closureHisto.Draw("colz text")
        # closureHisto.Draw("text same")
        lowerLimitString=str(lowerLimit).replace(".","p")
        canvas.SaveAs(outputDir+"2D_closure_data_SR_invertedChi2_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.png")
        canvas.SaveAs(outputDir+"2D_closure_data_SR_invertedChi2_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.pdf")
        bmClosureBin=closureHisto.GetXaxis().FindBin(bmPoint)
        bmClosure=closureHisto.GetBinContent(bmClosureBin, bmClosureBin)
        bmClosureList.append(bmClosure)
        bmClosure_error=closureHisto.GetBinError(bmClosureBin, bmClosureBin)
        bmClosureErrorList.append(bmClosure_error)
        bmLowerLimitList.append(lowerLimit)
    for i in range(len(bmClosureList)):
        bmGraph_data_SR_invertedChi2_0p89_0p89.SetPoint(i,bmLowerLimitList[i]+bmGraphXOffset, bmClosureList[i])
        bmGraph_data_SR_invertedChi2_0p89_0p89.SetPointError(i,0, bmClosureErrorList[i])
        global_max=max(global_max, max(bmClosureList))


    bmGraph_QCD_SR_invertedChi2_0p89_0p89=ROOT.TGraphErrors()
    bmGraph_QCD_SR_invertedChi2_0p89_0p89.SetLineColor(ROOT.kBlue)
    bmGraph_QCD_SR_invertedChi2_0p89_0p89.SetLineStyle(2)
    bmGraph_QCD_SR_invertedChi2_0p89_0p89.SetTitle("QCD SR inv. Chi2 ("+str(bmPoint)+","+str(bmPoint)+")")
    bmClosureList=[]
    bmClosureErrorList=[]
    bmLowerLimitList=[]
    bmGraphXOffset=0.005*5
    for lowerLimit in allLL:
        histoString=histoString_SR
        closureHisto, dum = scanClosure2D(
                                    windowMinX=lowerLimit, windowMaxX=1., windowMinY=lowerLimit, windowMaxY=1.0,
                                    scanMin=0.7, scanMax=0.96,
                                    histoFile=histoFile, histoString=histoString,
                                    dataMode="simQCD", verbosity=verb )
        canvas=ROOT.TCanvas("c_closure_QCD_SR_invertedChi2","c_closure_QCD_SR_invertedChi2",1024,768)
        ROOT.gStyle.SetOptStat(0)
        closureHisto.SetTitle("closure QCD SR inverted Chi2, ("+str(lowerLimit)+", 1.0) X ("+str(lowerLimit)+", 1.0)")
        closureHisto.SetMarkerSize(0.6)
        closureHisto.GetXaxis().SetTitle("AB vs CD cut")
        closureHisto.GetYaxis().SetTitle("AC vs BD cut")
        closureHisto.Draw("colz text")
        # closureHisto.Draw("text same")
        lowerLimitString=str(lowerLimit).replace(".","p")
        canvas.SaveAs(outputDir+"2D_closure_QCD_SR_invertedChi2_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.png")
        canvas.SaveAs(outputDir+"2D_closure_QCD_SR_invertedChi2_"+lowerLimitString+"_1p0_"+lowerLimitString+"_1p0.pdf")
        bmClosureBin=closureHisto.GetXaxis().FindBin(bmPoint)
        bmClosure=closureHisto.GetBinContent(bmClosureBin, bmClosureBin)
        bmClosureList.append(bmClosure)
        bmClosure_error=closureHisto.GetBinError(bmClosureBin, bmClosureBin)
        bmClosureErrorList.append(bmClosure_error)
        bmLowerLimitList.append(lowerLimit)
    for i in range(len(bmClosureList)):
        bmGraph_QCD_SR_invertedChi2_0p89_0p89.SetPoint(i,bmLowerLimitList[i]+bmGraphXOffset, bmClosureList[i])
        bmGraph_QCD_SR_invertedChi2_0p89_0p89.SetPointError(i,0, bmClosureErrorList[i])
        global_max=max(global_max, max(bmClosureList))


    canvas=ROOT.TCanvas("c_benchmark_0p89_0p89","c_benchmark_0p89_0p89", 1024, 768)
    bmGraph_data_VR_0p89_0p89.SetMaximum(global_max*1.1)
    bmGraph_data_VR_0p89_0p89.SetMinimum(0)
    bmGraph_data_VR_0p89_0p89.SetLineWidth(2)
    bmGraph_data_VR_0p89_0p89.GetXaxis().SetTitle("lower BCD cut")
    bmGraph_data_VR_0p89_0p89.GetYaxis().SetTitle("closure")
    bmGraph_data_VR_0p89_0p89.Draw("ALP 0")

    bmGraph_QCD_VR_0p89_0p89.SetLineWidth(2)
    bmGraph_QCD_VR_0p89_0p89.Draw("sameLP 0")

    bmGraph_QCD_SR_0p89_0p89.SetLineWidth(2)
    bmGraph_QCD_SR_0p89_0p89.Draw("sameLP 0")

    bmGraph_QCD_VR_invertedChi2_0p89_0p89.SetLineWidth(2)
    bmGraph_QCD_VR_invertedChi2_0p89_0p89.Draw("sameLP 0")

    bmGraph_data_VR_invertedChi2_0p89_0p89.SetLineWidth(2)
    bmGraph_data_VR_invertedChi2_0p89_0p89.Draw("sameLP 0")

    bmGraph_QCD_SR_invertedChi2_0p89_0p89.SetLineWidth(2)
    bmGraph_QCD_SR_invertedChi2_0p89_0p89.Draw("sameLP 0")

    bmGraph_data_SR_invertedChi2_0p89_0p89.SetLineWidth(2)
    bmGraph_data_SR_invertedChi2_0p89_0p89.Draw("sameLP 0")
    canvas.BuildLegend()
    canvas.SaveAs(outputDir+"leftRestricted_benchmark_"+str(bmPoint).replace(".","p")+"_"+str(bmPoint).replace(".","p")+".png")
    canvas.SaveAs(outputDir+"leftRestricted_benchmark_"+str(bmPoint).replace(".","p")+"_"+str(bmPoint).replace(".","p")+".pdf")