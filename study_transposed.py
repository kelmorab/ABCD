import ROOT
import math
import ctypes
import numpy
ROOT.gROOT.SetBatch(True)
from class_ABCD import ABCD



def scanClosureOneDirection(scanDirection="x",
                             windowMinX=0.0, windowMaxX=1., windowMinY=0.0, windowMaxY=1.0, otherDirectionCut=0.0,
                             scanMin=0.5, scanMax=0.9,
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

    cutValue_list=[]
    nonClosure_list=[]

    for cutValue in numpy.arange(scanMin, scanMax, 0.01):
        if verbosity>0:
            print(cutValue)
        if scanDirection=="x":
            theABCD.setWindow("A", cutValue, windowMaxX, otherDirectionCut, windowMaxY)
            theABCD.setWindow("B", cutValue, windowMaxX, windowMinY, otherDirectionCut-0.0001)
            theABCD.setWindow("C", windowMinX, cutValue-0.0001, otherDirectionCut, windowMaxY)
            theABCD.setWindow("D", windowMinX, cutValue-0.0001, windowMinY, otherDirectionCut-0.0001)
        if scanDirection=="y":
            print("???? not implemented yet")
            return None, None
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
            nonClosure=theABCD.getNonClosureWrtData()
            if verbosity>0: print("non-closure: "+str(round(nonClosure,3)))
            cutValue_list.append(cutValue)
            nonClosure_list.append(nonClosure)
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
            nonClosure=theABCD.getNonClosure()
            if verbosity>0: print("non-closure: "+str(round(nonClosure,3)))
            cutValue_list.append(cutValue)
            nonClosure_list.append(nonClosure)            

    return cutValue_list, nonClosure_list


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

            nonClosure=0
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
                nonClosure=theABCD.getNonClosureWrtData()
                if verbosity>0: print("non-closure: "+str(round(nonClosure,3)))

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
                nonClosure=theABCD.getNonClosure()
                if verbosity>0: print("non-closure: "+str(round(nonClosure,3)))

            binX=closureHisto.GetXaxis().FindBin(cutX)
            binY=closureHisto.GetYaxis().FindBin(cutY)
            closureHisto.SetBinContent(binX, binY, round(nonClosure,2))
            closuresList.append([round(cutX,2), round(cutY,2),round(nonClosure,2)])
    return closureHisto, closuresList

def scan2DTransposedABCDBinWise(histoFile=None, dataMode="data", 
                         blindSR=True):
    histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1"
    histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR"
    histoString_SR_invertedChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2"
    histoString_VR_invertedChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR"

    # histoFile=ROOT.TFile(histoFileName,"r")
    
    theABCD_SR=ABCD()
    theABCD_VR=ABCD()
    theABCD_SR_invertedChi2=ABCD()
    theABCD_VR_invertedChi2=ABCD()


    if dataMode=="data":
        histo_data_obs_Chi2_SR=histoFile.Get(histoString_SR.replace("$SAMPLE","data_obs"))
        histo_SM_Higgs_Chi2_SR = histoFile.Get(histoString_SR.replace("$SAMPLE","SM_Higgs"))
        histo_TTbar_Chi2_SR = histoFile.Get(histoString_SR.replace("$SAMPLE","TTbar"))
        histo_Others_Chi2_SR = histoFile.Get(histoString_SR.replace("$SAMPLE","Others"))
        histo_nonQCD_Chi2_SR = histo_SM_Higgs_Chi2_SR.Clone()
        histo_nonQCD_Chi2_SR.Add(histo_TTbar_Chi2_SR)
        histo_nonQCD_Chi2_SR.Add(histo_Others_Chi2_SR)
        theABCD_SR.setDataHisto(histo_data_obs_Chi2_SR)
        theABCD_SR.setNonQCDHisto(histo_nonQCD_Chi2_SR)
        theABCD_SR.subtractNonQCDFromData()

        histo_data_obs_Chi2_VR=histoFile.Get(histoString_VR.replace("$SAMPLE","data_obs"))
        histo_SM_Higgs_Chi2_VR = histoFile.Get(histoString_VR.replace("$SAMPLE","SM_Higgs"))
        histo_TTbar_Chi2_VR = histoFile.Get(histoString_VR.replace("$SAMPLE","TTbar"))
        histo_Others_Chi2_VR = histoFile.Get(histoString_VR.replace("$SAMPLE","Others"))
        histo_nonQCD_Chi2_VR = histo_SM_Higgs_Chi2_VR.Clone()
        histo_nonQCD_Chi2_VR.Add(histo_TTbar_Chi2_VR)
        histo_nonQCD_Chi2_VR.Add(histo_Others_Chi2_VR)
        theABCD_VR.setDataHisto(histo_data_obs_Chi2_VR)
        theABCD_VR.setNonQCDHisto(histo_nonQCD_Chi2_VR)
        theABCD_VR.subtractNonQCDFromData()

        histo_data_obs_Chi2_SR_invertedChi2=histoFile.Get(histoString_SR_invertedChi2.replace("$SAMPLE","data_obs"))
        histo_SM_Higgs_Chi2_SR_invertedChi2 = histoFile.Get(histoString_SR_invertedChi2.replace("$SAMPLE","SM_Higgs"))
        histo_TTbar_Chi2_SR_invertedChi2 = histoFile.Get(histoString_SR_invertedChi2.replace("$SAMPLE","TTbar"))
        histo_Others_Chi2_SR_invertedChi2 = histoFile.Get(histoString_SR_invertedChi2.replace("$SAMPLE","Others"))
        histo_nonQCD_Chi2_SR_invertedChi2 = histo_SM_Higgs_Chi2_SR_invertedChi2.Clone()
        histo_nonQCD_Chi2_SR_invertedChi2.Add(histo_TTbar_Chi2_SR_invertedChi2)
        histo_nonQCD_Chi2_SR_invertedChi2.Add(histo_Others_Chi2_SR_invertedChi2)
        theABCD_SR_invertedChi2.setDataHisto(histo_data_obs_Chi2_SR_invertedChi2)
        theABCD_SR_invertedChi2.setNonQCDHisto(histo_nonQCD_Chi2_SR_invertedChi2)
        theABCD_SR_invertedChi2.subtractNonQCDFromData()

        histo_data_obs_Chi2_VR_invertedChi2=histoFile.Get(histoString_VR_invertedChi2.replace("$SAMPLE","data_obs"))
        histo_SM_Higgs_Chi2_VR_invertedChi2 = histoFile.Get(histoString_VR_invertedChi2.replace("$SAMPLE","SM_Higgs"))
        histo_TTbar_Chi2_VR_invertedChi2 = histoFile.Get(histoString_VR_invertedChi2.replace("$SAMPLE","TTbar"))
        histo_Others_Chi2_VR_invertedChi2 = histoFile.Get(histoString_VR_invertedChi2.replace("$SAMPLE","Others"))
        histo_nonQCD_Chi2_VR_invertedChi2 = histo_SM_Higgs_Chi2_VR_invertedChi2.Clone()
        histo_nonQCD_Chi2_VR_invertedChi2.Add(histo_TTbar_Chi2_VR_invertedChi2)
        histo_nonQCD_Chi2_VR_invertedChi2.Add(histo_Others_Chi2_VR_invertedChi2)
        theABCD_VR_invertedChi2.setDataHisto(histo_data_obs_Chi2_VR_invertedChi2)
        theABCD_VR_invertedChi2.setNonQCDHisto(histo_nonQCD_Chi2_VR_invertedChi2)
        theABCD_VR_invertedChi2.subtractNonQCDFromData()

    elif dataMode=="simQCD":
        histo_QCD_Chi2_SR=histoFile.Get(histoString_SR.replace("$SAMPLE","QCD"))
        theABCD_SR.setQCDHisto(histo_QCD_Chi2_SR)

        histo_QCD_Chi2_VR=histoFile.Get(histoString_VR.replace("$SAMPLE","QCD"))
        theABCD_VR.setQCDHisto(histo_QCD_Chi2_VR)

        histo_QCD_Chi2_SR_invertedChi2=histoFile.Get(histoString_SR_invertedChi2.replace("$SAMPLE","QCD"))
        theABCD_SR_invertedChi2.setQCDHisto(histo_QCD_Chi2_SR_invertedChi2)

        histo_QCD_Chi2_VR_invertedChi2=histoFile.Get(histoString_VR_invertedChi2.replace("$SAMPLE","QCD"))
        theABCD_VR_invertedChi2.setQCDHisto(histo_QCD_Chi2_VR_invertedChi2)


    # for now binwise
    # closureHistoSR=ROOT.TH2D("clousure_SR","closure_SR",100,0,1,100,0,1)
    closureHistoSR=theABCD_SR.histo.Clone()
    closureHistoSR.Reset()
    # closureHistoSR.SetDirectory(0)
    print(closureHistoSR)
    for xbin in range(closureHistoSR.GetXaxis().GetNbins()):
        for ybin in range(closureHistoSR.GetYaxis().GetNbins()):
            closureHistoSR.SetBinContent(xbin,ybin,0.0)
            x = closureHistoSR.GetXaxis().GetBinCenter(xbin)
            y = closureHistoSR.GetYaxis().GetBinCenter(ybin)
            if dataMode=="data" and blindSR and x>0.7 and y>0.7:
                continue
            A=theABCD_SR.histo.GetBinContent(xbin,ybin)
            B=theABCD_SR_invertedChi2.histo.GetBinContent(xbin,ybin)
            C=theABCD_VR.histo.GetBinContent(xbin,ybin)
            D=theABCD_VR_invertedChi2.histo.GetBinContent(xbin,ybin)
            predA=0
            if D>0:
                predA=B*C/D
            closure=0
            if A>0:
                closure=predA/A
            closureHistoSR.SetBinContent(xbin,ybin,round(closure,2))    
    return closureHistoSR


def scan2DTransposedABCD(histoFile=None, dataMode="data", windowSize=0.1,
                         windowMinX=0.0, windowMaxX=1., windowMinY=0.0, windowMaxY=1.0,
                         blindSR=True):
    histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1"
    histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR"
    histoString_SR_invertedChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2"
    histoString_VR_invertedChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR"

    # histoFile=ROOT.TFile(histoFileName,"r")
    
    theABCD_SR=ABCD()
    theABCD_VR=ABCD()
    theABCD_SR_invertedChi2=ABCD()
    theABCD_VR_invertedChi2=ABCD()


    if dataMode=="data":
        histo_data_obs_Chi2_SR=histoFile.Get(histoString_SR.replace("$SAMPLE","data_obs"))
        histo_SM_Higgs_Chi2_SR = histoFile.Get(histoString_SR.replace("$SAMPLE","SM_Higgs"))
        histo_TTbar_Chi2_SR = histoFile.Get(histoString_SR.replace("$SAMPLE","TTbar"))
        histo_Others_Chi2_SR = histoFile.Get(histoString_SR.replace("$SAMPLE","Others"))
        histo_nonQCD_Chi2_SR = histo_SM_Higgs_Chi2_SR.Clone()
        histo_nonQCD_Chi2_SR.Add(histo_TTbar_Chi2_SR)
        histo_nonQCD_Chi2_SR.Add(histo_Others_Chi2_SR)
        theABCD_SR.setDataHisto(histo_data_obs_Chi2_SR)
        theABCD_SR.setNonQCDHisto(histo_nonQCD_Chi2_SR)
        theABCD_SR.subtractNonQCDFromData()

        histo_data_obs_Chi2_VR=histoFile.Get(histoString_VR.replace("$SAMPLE","data_obs"))
        histo_SM_Higgs_Chi2_VR = histoFile.Get(histoString_VR.replace("$SAMPLE","SM_Higgs"))
        histo_TTbar_Chi2_VR = histoFile.Get(histoString_VR.replace("$SAMPLE","TTbar"))
        histo_Others_Chi2_VR = histoFile.Get(histoString_VR.replace("$SAMPLE","Others"))
        histo_nonQCD_Chi2_VR = histo_SM_Higgs_Chi2_VR.Clone()
        histo_nonQCD_Chi2_VR.Add(histo_TTbar_Chi2_VR)
        histo_nonQCD_Chi2_VR.Add(histo_Others_Chi2_VR)
        theABCD_VR.setDataHisto(histo_data_obs_Chi2_VR)
        theABCD_VR.setNonQCDHisto(histo_nonQCD_Chi2_VR)
        theABCD_VR.subtractNonQCDFromData()

        histo_data_obs_Chi2_SR_invertedChi2=histoFile.Get(histoString_SR_invertedChi2.replace("$SAMPLE","data_obs"))
        histo_SM_Higgs_Chi2_SR_invertedChi2 = histoFile.Get(histoString_SR_invertedChi2.replace("$SAMPLE","SM_Higgs"))
        histo_TTbar_Chi2_SR_invertedChi2 = histoFile.Get(histoString_SR_invertedChi2.replace("$SAMPLE","TTbar"))
        histo_Others_Chi2_SR_invertedChi2 = histoFile.Get(histoString_SR_invertedChi2.replace("$SAMPLE","Others"))
        histo_nonQCD_Chi2_SR_invertedChi2 = histo_SM_Higgs_Chi2_SR_invertedChi2.Clone()
        histo_nonQCD_Chi2_SR_invertedChi2.Add(histo_TTbar_Chi2_SR_invertedChi2)
        histo_nonQCD_Chi2_SR_invertedChi2.Add(histo_Others_Chi2_SR_invertedChi2)
        theABCD_SR_invertedChi2.setDataHisto(histo_data_obs_Chi2_SR_invertedChi2)
        theABCD_SR_invertedChi2.setNonQCDHisto(histo_nonQCD_Chi2_SR_invertedChi2)
        theABCD_SR_invertedChi2.subtractNonQCDFromData()

        histo_data_obs_Chi2_VR_invertedChi2=histoFile.Get(histoString_VR_invertedChi2.replace("$SAMPLE","data_obs"))
        histo_SM_Higgs_Chi2_VR_invertedChi2 = histoFile.Get(histoString_VR_invertedChi2.replace("$SAMPLE","SM_Higgs"))
        histo_TTbar_Chi2_VR_invertedChi2 = histoFile.Get(histoString_VR_invertedChi2.replace("$SAMPLE","TTbar"))
        histo_Others_Chi2_VR_invertedChi2 = histoFile.Get(histoString_VR_invertedChi2.replace("$SAMPLE","Others"))
        histo_nonQCD_Chi2_VR_invertedChi2 = histo_SM_Higgs_Chi2_VR_invertedChi2.Clone()
        histo_nonQCD_Chi2_VR_invertedChi2.Add(histo_TTbar_Chi2_VR_invertedChi2)
        histo_nonQCD_Chi2_VR_invertedChi2.Add(histo_Others_Chi2_VR_invertedChi2)
        theABCD_VR_invertedChi2.setDataHisto(histo_data_obs_Chi2_VR_invertedChi2)
        theABCD_VR_invertedChi2.setNonQCDHisto(histo_nonQCD_Chi2_VR_invertedChi2)
        theABCD_VR_invertedChi2.subtractNonQCDFromData()

    elif dataMode=="simQCD":
        histo_QCD_Chi2_SR=histoFile.Get(histoString_SR.replace("$SAMPLE","QCD"))
        theABCD_SR.setQCDHisto(histo_QCD_Chi2_SR)

        histo_QCD_Chi2_VR=histoFile.Get(histoString_VR.replace("$SAMPLE","QCD"))
        theABCD_VR.setQCDHisto(histo_QCD_Chi2_VR)

        histo_QCD_Chi2_SR_invertedChi2=histoFile.Get(histoString_SR_invertedChi2.replace("$SAMPLE","QCD"))
        theABCD_SR_invertedChi2.setQCDHisto(histo_QCD_Chi2_SR_invertedChi2)

        histo_QCD_Chi2_VR_invertedChi2=histoFile.Get(histoString_VR_invertedChi2.replace("$SAMPLE","QCD"))
        theABCD_VR_invertedChi2.setQCDHisto(histo_QCD_Chi2_VR_invertedChi2)


    # for now binwise
    # closureHistoSR=ROOT.TH2D("clousure_SR","closure_SR",100,0,1,100,0,1)
    nBinsX=int((windowMaxX-windowMinX)/windowSize)
    nBinsY=int((windowMaxY-windowMinY)/windowSize)
    closureHistoSR=ROOT.TH2D("closure","closure",nBinsX,windowMinX,windowMaxX,nBinsY,windowMinY,windowMaxY)
    # closureHistoSR.SetDirectory(0)
    print(closureHistoSR)
    for xbin in range(closureHistoSR.GetXaxis().GetNbins()):
        for ybin in range(closureHistoSR.GetYaxis().GetNbins()):
            closureHistoSR.SetBinContent(xbin,ybin,0.0)
            # xlow = closureHistoSR.GetXaxis().GetBinLowEdge(xbin)
            # ylow = closureHistoSR.GetYaxis().GetBinLowEdge(ybin)
            if dataMode=="data" and blindSR and xbin>=closureHistoSR.GetXaxis().FindBin(0.701) and ybin>=closureHistoSR.GetYaxis().FindBin(0.701):
                continue
            theABCD_SR.setWindow("A",closureHistoSR.GetXaxis().GetBinLowEdge(xbin),closureHistoSR.GetXaxis().GetBinUpEdge(xbin), closureHistoSR.GetYaxis().GetBinLowEdge(ybin), closureHistoSR.GetYaxis().GetBinUpEdge(ybin))
            A, Aerr=theABCD_SR.getIntegralAndErrorInRegion("A","QCD")
            theABCD_SR.printWindow("A")
            theABCD_SR_invertedChi2.setWindow("A",closureHistoSR.GetXaxis().GetBinLowEdge(xbin),closureHistoSR.GetXaxis().GetBinUpEdge(xbin), closureHistoSR.GetYaxis().GetBinLowEdge(ybin), closureHistoSR.GetYaxis().GetBinUpEdge(ybin))
            B, Berr=theABCD_SR_invertedChi2.getIntegralAndErrorInRegion("A","QCD")
            theABCD_SR_invertedChi2.printWindow("A")
            theABCD_VR.setWindow("A",closureHistoSR.GetXaxis().GetBinLowEdge(xbin),closureHistoSR.GetXaxis().GetBinUpEdge(xbin), closureHistoSR.GetYaxis().GetBinLowEdge(ybin), closureHistoSR.GetYaxis().GetBinUpEdge(ybin))
            C,Cerr=theABCD_VR.getIntegralAndErrorInRegion("A","QCD")
            theABCD_VR.printWindow("A")
            theABCD_VR_invertedChi2.setWindow("A",closureHistoSR.GetXaxis().GetBinLowEdge(xbin),closureHistoSR.GetXaxis().GetBinUpEdge(xbin), closureHistoSR.GetYaxis().GetBinLowEdge(ybin), closureHistoSR.GetYaxis().GetBinUpEdge(ybin))
            D, Derr=theABCD_VR_invertedChi2.getIntegralAndErrorInRegion("A","QCD")
            theABCD_VR_invertedChi2.printWindow("A")
            
            predA=0
            if D>0:
                predA=B*C/D
            closure=0
            if A>0:
                closure=predA/A
            print(A,B,C,D,predA,closure)    
            closureHistoSR.SetBinContent(xbin,ybin,round(closure,2))    
    return closureHistoSR

def analyzeClosureFrequency(histo):
    nBinsX=histo.GetXaxis().GetNbins()
    nBinsY=histo.GetYaxis().GetNbins()
    freqMax=10
    h_freq=ROOT.TH1D("freq","freq",100,0,freqMax)
    for binX in range(nBinsX):
        for binY in range(nBinsY):
            closure=histo.GetBinContent(binX,binY)
            h_freq.Fill(closure)
    h_freq.SetBinContent(100,h_freq.GetBinContent(100)+h_freq.GetBinContent(101))        
    h_freq.SetBinContent(101,0)
    return h_freq

    




ROOT.gStyle.SetOptStat(0)
histoFileName="all_histos_main_v6_invertedChi2.root"
histoFile=ROOT.TFile(histoFileName,"r")
outputDir="transposed_ABCD_M15_1/"

closureHisto_data=scan2DTransposedABCDBinWise(histoFile=histoFile, dataMode="data",
                         blindSR=True)
canvas=ROOT.TCanvas("closure_data","closure_data",1024,786)
closureHisto_data.SetTitle("binwise closure, transposed ABCD, data")
closureHisto_data.SetMarkerSize(0.6)
closureHisto_data.Draw("colz")
canvas.SaveAs(outputDir+"2D_binwise_trans_closure_M15_1_data.pdf")
canvas.SaveAs(outputDir+"2D_binwise_trans_closure_M15_1_data.png")
h_freq_data=analyzeClosureFrequency(closureHisto_data)
canvas=ROOT.TCanvas("closure_binwise_frequency_data","closure_binwise_frequency_data",1024,786)
h_freq_data.Draw("histo")
canvas.SaveAs(outputDir+"freq_binwise_trans_closure_M15_1_data.pdf")
canvas.SaveAs(outputDir+"freq_binwise_trans_closure_M15_1_data.png")

closureHisto_QCD=scan2DTransposedABCDBinWise(histoFile=histoFile, dataMode="simQCD", 
                         blindSR=True)
canvas=ROOT.TCanvas("closure_QCD","closure_QCD",1024,786)
closureHisto_QCD.SetTitle("binwise closure, transposed ABCD, QCD")
closureHisto_QCD.SetMarkerSize(0.6)
closureHisto_QCD.Draw("colz")
canvas.SaveAs(outputDir+"2D_binwise_trans_closure_M15_1_QCD.pdf")
canvas.SaveAs(outputDir+"2D_binwise_trans_closure_M15_1_QCD.png")


print("-------------------")
print("data 0.05")
closureHisto_data=scan2DTransposedABCD(histoFile=histoFile, dataMode="data", windowSize=0.05,
                         windowMinX=0.0, windowMaxX=1., windowMinY=0.0, windowMaxY=1.0,
                         blindSR=True)
canvas=ROOT.TCanvas("closure_data","closure_data",1024,786)
closureHisto_data.SetTitle("0.05w closure, transposed ABCD, data")
closureHisto_data.SetMarkerSize(0.6)
closureHisto_data.Draw("colz text")
canvas.SaveAs(outputDir+"2D_0p05_trans_closure_M15_1_data.pdf")
canvas.SaveAs(outputDir+"2D_0p05_trans_closure_M15_1_data.png")

print("-------------------")
print("QCD 0.05")
closureHisto_QCD=scan2DTransposedABCD(histoFile=histoFile, dataMode="simQCD",  windowSize=0.05,
                         windowMinX=0.0, windowMaxX=1., windowMinY=0.0, windowMaxY=1.0,
                         blindSR=True)
canvas=ROOT.TCanvas("closure_QCD","closure_QCD",1024,786)
closureHisto_QCD.SetTitle("0.05w closure, transposed ABCD, QCD")
closureHisto_QCD.SetMarkerSize(0.6)
closureHisto_QCD.Draw("colz text")
canvas.SaveAs(outputDir+"2D_0p05_trans_closure_M15_1_QCD.pdf")
canvas.SaveAs(outputDir+"2D_0p05_trans_closure_M15_1_QCD.png")

print("-------------------")
print("data 0.1")
closureHisto_data=scan2DTransposedABCD(histoFile=histoFile, dataMode="data", windowSize=0.1,
                         windowMinX=0.0, windowMaxX=1., windowMinY=0.0, windowMaxY=1.0,
                         blindSR=True)
canvas=ROOT.TCanvas("closure_data","closure_data",1024,786)
closureHisto_data.SetTitle("0.1w closure, transposed ABCD, data")
closureHisto_data.SetMarkerSize(0.6)
closureHisto_data.Draw("colz text")
canvas.SaveAs(outputDir+"2D_0p1_trans_closure_M15_1_data.pdf")
canvas.SaveAs(outputDir+"2D_0p1_trans_closure_M15_1_data.png")

print("-------------------")
print("QCD 0.1")
closureHisto_QCD=scan2DTransposedABCD(histoFile=histoFile, dataMode="simQCD",  windowSize=0.1,
                         windowMinX=0.0, windowMaxX=1., windowMinY=0.0, windowMaxY=1.0,
                         blindSR=True)
canvas=ROOT.TCanvas("closure_QCD","closure_QCD",1024,786)
closureHisto_QCD.SetTitle("0.1w closure, transposed ABCD, QCD")
closureHisto_QCD.SetMarkerSize(0.6)
closureHisto_QCD.Draw("colz text")
canvas.SaveAs(outputDir+"2D_0p1_trans_closure_M15_1_QCD.pdf")
canvas.SaveAs(outputDir+"2D_0p1_trans_closure_M15_1_QCD.png")
