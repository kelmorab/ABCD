import ROOT
import math
import ctypes
import numpy
import os
ROOT.gROOT.SetBatch(True)
ROOT.gDirectory.cd('PyROOT:/')
from class_ABCD import ABCD, FullABCD


def scanClosureOneDirection_NewScheme(scanDirection="x",
                             windowMinX=0.0, windowMaxX=1., windowMinY=0.0, windowMaxY=1.0, otherDirectionCut=0.0,
                             scanMin=0.5, scanMax=0.9,
                             filename="",
                             histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1",
                             histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR",
                             histoString_SR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2",
                             histoString_VR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR",
                             dataMode="data", verbosity=0, applyAdditionalUncertainty=True ):
    
    histoString_SR=histoString_SR
    histoString_VR=histoString_VR
    histoString_SR_invChi2=histoString_SR_invChi2
    histoString_VR_invChi2=histoString_VR_invChi2

    theFullABCD = FullABCD()
    theFullABCD.setup_ABCD("SR",filename,histoString_SR, dataMode )
    theFullABCD.setup_ABCD("VR",filename,histoString_VR, dataMode )
    theFullABCD.setup_ABCD("SR_prime",filename,histoString_SR_invChi2, dataMode )
    theFullABCD.setup_ABCD("VR_prime",filename,histoString_VR_invChi2, dataMode )

    cutValue_list=[]
    nonClosure_list=[]
    non_Closure_error_list=[]

    for cutValue in numpy.arange(scanMin, scanMax, 0.01):
        if verbosity>0:
            print(cutValue)
        if scanDirection=="x":
            theFullABCD.setWindow("A", cutValue, windowMaxX, otherDirectionCut, windowMaxY)
            theFullABCD.setWindow("B", cutValue, windowMaxX, windowMinY, otherDirectionCut-0.0001)
            theFullABCD.setWindow("C", windowMinX, cutValue-0.0001, otherDirectionCut, windowMaxY)
            theFullABCD.setWindow("D", windowMinX, cutValue-0.0001, windowMinY, otherDirectionCut-0.0001)
        if scanDirection=="y":
            print("???? not implemented yet")
            return None, None, None
        if verbosity>0:
            theFullABCD.theABCD_SR.printWindow("A")
            theFullABCD.theABCD_SR.printWindow("B")
            theFullABCD.theABCD_SR.printWindow("C")
            theFullABCD.theABCD_SR.printWindow("D")
            print("Integral of whole histo: "+str(theFullABCD.theABCD_SR.getIntegralAndErrorInWholeHisto()[0]))
            print("Sum of window Integrals: "+str(theFullABCD.theABCD_SR.getIntegralAndErrorInRegion("A")[0]+theFullABCD.theABCD_SR.getIntegralAndErrorInRegion("B")[0]+theFullABCD.theABCD_SR.getIntegralAndErrorInRegion("C")[0]+theFullABCD.theABCD_SR.getIntegralAndErrorInRegion("D")[0]))
            print("")
        
        nonClosure, nonClosure_error = theFullABCD.getFullClosureAndUncertainty(dataMode, verbosity, applyAdditionalUncertainty=applyAdditionalUncertainty)
        cutValue_list.append(cutValue)
        nonClosure_list.append(nonClosure)
        non_Closure_error_list.append(nonClosure_error)
        if verbosity>0:
            print(nonClosure, nonClosure_error)         

    return cutValue_list, nonClosure_list, non_Closure_error_list


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



###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

# Scan closure = predicted/observed for different lower_bounds, upper bound and ABCD cuts


# create toplevel dir
toplevelDirName="outputs_comparison_sideband_closure_allRegions_New"
if not os.path.isdir(toplevelDirName):
    os.mkdir(toplevelDirName)

masses=["15","40","55"]
lifetimes=["0","0p05","0p1","1",]

histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
filename="all_histos_main_v6_invertedChi2.root"

do_VRs=False

for mass in masses:
    for lifetime in lifetimes:
        WP_tagName="M"+mass+"_"+lifetime
        print("")
        print("At "+WP_tagName)
        outputDir=toplevelDirName+"/"+"comparison_sideband_x-scan_"+WP_tagName+"/"
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)

        if mass=="15":
            histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG"
            histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_VR"
            histoString_SR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_invertedChi2"
            histoString_VR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_invertedChi2_VR"
        if mass=="40" or mass=="55":
            histoString_SR="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG"
            histoString_VR="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_VR"
            histoString_SR_invChi2="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_invertedChi2"
            histoString_VR_invChi2="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_invertedChi2_VR"
        histoString_SR=histoString_SR.replace("$WPTAG",WP_tagName)
        histoString_VR=histoString_VR.replace("$WPTAG",WP_tagName)
        histoString_SR_invChi2=histoString_SR_invChi2.replace("$WPTAG",WP_tagName)
        histoString_VR_invChi2=histoString_VR_invChi2.replace("$WPTAG",WP_tagName)
        
        lowerbound_list=[0.0, 0.3, 0.4, 0.5, 0.6,]
        for lowerbound in lowerbound_list:
            global_max=0.0
            global_min=1.0

            windowMinX=lowerbound
            windowMinY=lowerbound
            windowMaxY=1.0
            scanDirection="x"
            scanMin=lowerbound+0.01

            scanMaxList=[ 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69]
            otherDirectionCutList=[ 0.7,  0.83, 0.85, 0.87, 0.89, 0.91, 0.93, 0.95]
            windowMaxXList=[  0.7 for s in scanMaxList]

            for i in range(len(scanMaxList)):
                print(i)
                verb=0

                graphXOffset=0.0000*i
                graphsList=[]
                global_max=0.
                global_min=10.
                colorList=[   ROOT.kRed,  ROOT.kMagenta,  ]
                modes_list=["data", "data", ]
                region_list=["SR","SR inv Chi2", ]
                histoStringList=[histoString_SR, histoString_SR_invChi2,]
                if do_VRs:
                    colorList=[   ROOT.kRed,  ROOT.kMagenta, ROOT.kBlue, ROOT.kCyan,  ]
                    modes_list=["data", "data", "data", "data",]
                    region_list=["SR","SR inv Chi2", "VR", "VR inv Chi2" ]                
                    histoStringList=[histoString_SR, histoString_SR_invChi2, histoString_VR, histoString_VR_invChi2]

                for theRegion, theHistoString, theColor, theMode in zip(region_list, histoStringList, colorList, modes_list):
                    print("At", mass, lifetime, "loop counter", i, theRegion, theHistoString, theColor, theMode)
                    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                                            windowMinX=windowMinX, windowMaxX=windowMaxXList[i], windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                                            scanMin=scanMin, scanMax=scanMaxList[i],
                                            histoFile=histoFile,
                                            histoString=theHistoString, 
                                            dataMode=theMode, verbosity=verb )
                    graphsList.append(ROOT.TGraphErrors())
                    for k in range(len(cutValue_list)):
                        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
                        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
                    graphsList[-1].SetLineColor(theColor)
                    graphsList[-1].SetName(theRegion+" "+theMode+" window=["+str(windowMinX)+","+str(windowMaxXList[i])+"]x["+str(windowMinY)+","+str(windowMaxY)+"], yCut="+str(otherDirectionCutList[i]))
                    graphsList[-1].GetXaxis().SetTitle("ABCD x Cut")
                    graphsList[-1].GetYaxis().SetTitle("prediction/observed")
                    if theMode=="simQCD":
                        graphsList[-1].SetLineStyle(2)
                    else:
                        graphsList[-1].SetLineStyle(1)    
                    global_max=max(global_max,max(nonClosure_list))
                    global_min=min(global_min, min(nonClosure_list))

                canvas=ROOT.TCanvas("canvas_M15_1_xWindow","canvas_M15_1_xWindow",1024,768)
                global_max=min(2.5, global_max)
                global_min=min(0.9, global_min)
                # global_min=0.8
                # print(graphsList)
                # print(graphsList[0])
                # graphsList[0].Print()
                print("global max", global_max)
                graphsList[0].SetMaximum(global_max*1.1)
                graphsList[0].SetMinimum(global_min*0.9)
                graphsList[0].GetXaxis().SetLimits(scanMin-0.05, windowMaxXList[i])
                graphsList[0].Draw("ALP0")
                for gr in graphsList[1:]:
                    gr.Draw("same LP0")
                canvas.BuildLegend()
                canvas.SetGrid()
                windowString="["+str(windowMinX)+","+str(windowMaxXList[i])+"]x["+str(windowMinY)+","+str(windowMaxY)+"]"
                windowString=windowString.replace(".","p")
                canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_"+windowString+"_yCut_"+str(otherDirectionCutList[i])+"_"+WP_tagName+".png")
                canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_"+windowString+"_yCut_"+str(otherDirectionCutList[i])+"_"+WP_tagName+".pdf")


            # Now restrict y<0.7 and scan whole x-range

            global_max=0.0
            global_min=1.0

            windowMinX=lowerbound
            windowMaxX=1.
            windowMinY=lowerbound
            windowMaxY=0.7
            scanDirection="x"
            scanMin=lowerbound+0.01

            scanMaxList=[  0.95, 0.95, ]
            otherDirectionCutList=[  0.65, 0.67,  ]

            graphXOffset=0.0000*1

            for i in range(len(scanMaxList)):
                print(i)
                verb=0
                graphXOffset=0.0000*i
                graphsList=[]
                global_max=0.
                global_min=10.
                colorList=[   ROOT.kRed,  ROOT.kMagenta,  ]
                modes_list=["data", "data", ]
                region_list=["SR","SR inv Chi2", ]
                histoStringList=[histoString_SR, histoString_SR_invChi2,]
                if do_VRs:
                    colorList=[   ROOT.kRed,  ROOT.kMagenta, ROOT.kBlue, ROOT.kCyan,  ]
                    modes_list=["data", "data", "data", "data",]
                    region_list=["SR","SR inv Chi2", "VR", "VR inv Chi2" ]                
                    histoStringList=[histoString_SR, histoString_SR_invChi2, histoString_VR, histoString_VR_invChi2]
                
                for theRegion, theHistoString, theColor, theMode in zip(region_list, histoStringList, colorList, modes_list):
                    print("At", mass, lifetime, "loop counter", i, theRegion, theHistoString, theColor, theMode)
                    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                                            windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                                            scanMin=scanMin, scanMax=scanMaxList[i],
                                            histoFile=histoFile,
                                            histoString=theHistoString, 
                                            dataMode=theMode, verbosity=verb )
                    graphsList.append(ROOT.TGraphErrors())
                    for k in range(len(cutValue_list)):
                        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
                        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
                    graphsList[-1].SetLineColor(theColor)
                    graphsList[-1].SetName(theRegion+" "+theMode+" window=["+str(windowMinX)+","+str(windowMaxX)+"]x["+str(windowMinY)+","+str(windowMaxY)+"], yCut="+str(otherDirectionCutList[i]))
                    graphsList[-1].GetXaxis().SetTitle("ABCD x Cut")
                    graphsList[-1].GetYaxis().SetTitle("predicted/observed")
                    if theMode=="simQCD":
                        graphsList[-1].SetLineStyle(2)
                    else:
                        graphsList[-1].SetLineStyle(1)    
                    global_max=max(global_max,max(nonClosure_list))
                    global_min=min(global_min, min(nonClosure_list))

                canvas=ROOT.TCanvas("canvas_M15_1_yWindow","canvas_M15_1_yWindow",1024,768)
                global_max=min(2.5, global_max)
                print("global max", global_max)
                global_min=min(0.9, global_min)
                # global_min=0.8
                # print(graphsList)
                # print(graphsList[0])
                # graphsList[0].Print()
                graphsList[0].SetMaximum(global_max*1.1)
                graphsList[0].SetMinimum(global_min*0.9)
                graphsList[0].GetXaxis().SetLimits(scanMin-0.05, 1.0)
                graphsList[0].Draw("ALP0")
                for gr in graphsList[1:]:
                    gr.Draw("same LP0")
                canvas.BuildLegend()
                canvas.SetGrid()
                windowString="["+str(windowMinX)+","+str(windowMaxX)+"]x["+str(windowMinY)+","+str(windowMaxY)+"]"
                windowString=windowString.replace(".","p")
                canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_"+windowString+"_yCut_"+str(otherDirectionCutList[i])+"_"+WP_tagName+".png")
                canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_"+windowString+"_yCut_"+str(otherDirectionCutList[i])+"_"+WP_tagName+".pdf")

###########################################################################################################################################################
###########################################################################################################################################################

        # now scan in y direction
        # restrict in y-direction
        print("")
        print("At "+WP_tagName)
        outputDir=toplevelDirName+"/"+"comparison_sideband_y-scan_"+WP_tagName+"/"
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)

        if mass=="15":
            histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG"
            histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_VR"
            histoString_SR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_invertedChi2"
            histoString_VR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_invertedChi2_VR"
        if mass=="40" or mass=="55":
            histoString_SR="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG"
            histoString_VR="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_VR"
            histoString_SR_invChi2="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_invertedChi2"
            histoString_VR_invChi2="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_invertedChi2_VR"
        histoString_SR=histoString_SR.replace("$WPTAG",WP_tagName)
        histoString_VR=histoString_VR.replace("$WPTAG",WP_tagName)
        histoString_SR_invChi2=histoString_SR_invChi2.replace("$WPTAG",WP_tagName)
        histoString_VR_invChi2=histoString_VR_invChi2.replace("$WPTAG",WP_tagName)
        
        lowerbound_list=[0.0, 0.3, 0.4, 0.5, 0.6,]
        for lowerbound in lowerbound_list:
            global_max=0.0
            global_min=1.0

            windowMinX=lowerbound
            windowMinY=lowerbound
            windowMaxX=1.0
            scanDirection="y"
            scanMin=lowerbound+0.01

            scanMaxList=[ 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69]
            otherDirectionCutList=[ 0.7,  0.83, 0.85, 0.87, 0.89, 0.91, 0.93, 0.95]
            windowMaxYList=[  0.7 for s in scanMaxList]

            for i in range(len(scanMaxList)):
                print(i)
                verb=0
                graphXOffset=0.0000*i
                graphsList=[]
                global_max=0.
                global_min=10.
                colorList=[   ROOT.kRed,  ROOT.kMagenta,  ]
                modes_list=["data", "data", ]
                region_list=["SR","SR inv Chi2", ]
                histoStringList=[histoString_SR, histoString_SR_invChi2,]
                if do_VRs:
                    colorList=[   ROOT.kRed,  ROOT.kMagenta, ROOT.kBlue, ROOT.kCyan,  ]
                    modes_list=["data", "data", "data", "data",]
                    region_list=["SR","SR inv Chi2", "VR", "VR inv Chi2" ]                
                    histoStringList=[histoString_SR, histoString_SR_invChi2, histoString_VR, histoString_VR_invChi2]

                for theRegion, theHistoString, theColor, theMode in zip(region_list, histoStringList, colorList, modes_list):
                    print("At", mass, lifetime, "loop counter", i, theRegion, theHistoString, theColor, theMode)
                    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                                            windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinY, windowMaxY=windowMaxYList[i], otherDirectionCut=otherDirectionCutList[i],
                                            scanMin=scanMin, scanMax=scanMaxList[i],
                                            histoFile=histoFile,
                                            histoString=theHistoString, 
                                            dataMode=theMode, verbosity=verb )
                    graphsList.append(ROOT.TGraphErrors())
                    for k in range(len(cutValue_list)):
                        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
                        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
                    graphsList[-1].SetLineColor(theColor)
                    graphsList[-1].SetName(theRegion+" "+theMode+" window=["+str(windowMinX)+","+str(windowMaxX)+"]x["+str(windowMinY)+","+str(windowMaxYList[i])+"], xCut="+str(otherDirectionCutList[i]))
                    graphsList[-1].GetXaxis().SetTitle("ABCD y Cut")
                    graphsList[-1].GetYaxis().SetTitle("prediction/observed")
                    if theMode=="simQCD":
                        graphsList[-1].SetLineStyle(2)
                    else:
                        graphsList[-1].SetLineStyle(1)    
                    global_max=max(global_max,max(nonClosure_list))
                    global_min=min(global_min, min(nonClosure_list))

                canvas=ROOT.TCanvas("canvas_M15_1_xWindow","canvas_M15_1_xWindow",1024,768)
                global_max=min(2.5, global_max)
                global_min=min(0.9, global_min)
                # global_min=0.8
                # print(graphsList)
                # print(graphsList[0])
                # graphsList[0].Print()
                print("global max", global_max)
                graphsList[0].SetMaximum(global_max*1.1)
                graphsList[0].SetMinimum(global_min*0.9)
                graphsList[0].GetXaxis().SetLimits(scanMin-0.05, windowMaxYList[i])
                graphsList[0].Draw("ALP0")
                for gr in graphsList[1:]:
                    gr.Draw("same LP0")
                canvas.BuildLegend()
                canvas.SetGrid()
                windowString="["+str(windowMinX)+","+str(windowMaxX)+"]x["+str(windowMinY)+","+str(windowMaxYList[i])+"]"
                windowString=windowString.replace(".","p")
                canvas.SaveAs(outputDir+"scan_yDirection_yRestricted_"+windowString+"_xCut_"+str(otherDirectionCutList[i])+"_"+WP_tagName+".png")
                canvas.SaveAs(outputDir+"scan_yDirection_yRestricted_"+windowString+"_xCut_"+str(otherDirectionCutList[i])+"_"+WP_tagName+".pdf")


            # Now restrict x<0.7 and scan whole y-range

            global_max=0.0
            global_min=1.0

            windowMinX=lowerbound
            windowMaxX=0.7
            windowMinY=lowerbound
            windowMaxY=1.0
            scanDirection="y"
            scanMin=lowerbound+0.01

            scanMaxList=[  0.95, 0.95, ]
            otherDirectionCutList=[  0.65, 0.67,  ]

            graphXOffset=0.0000*1

            for i in range(len(scanMaxList)):
                print(i)
                verb=0
                graphXOffset=0.0000*i
                graphsList=[]
                global_max=0.
                global_min=10.
                colorList=[   ROOT.kRed,  ROOT.kMagenta,  ]
                modes_list=["data", "data", ]
                region_list=["SR","SR inv Chi2", ]
                histoStringList=[histoString_SR, histoString_SR_invChi2,]
                if do_VRs:
                    colorList=[   ROOT.kRed,  ROOT.kMagenta, ROOT.kBlue, ROOT.kCyan,  ]
                    modes_list=["data", "data", "data", "data",]
                    region_list=["SR","SR inv Chi2", "VR", "VR inv Chi2" ]                
                    histoStringList=[histoString_SR, histoString_SR_invChi2, histoString_VR, histoString_VR_invChi2]
                
                for theRegion, theHistoString, theColor, theMode in zip(region_list, histoStringList, colorList, modes_list):
                    print("At", mass, lifetime, "loop counter", i, theRegion, theHistoString, theColor, theMode)
                    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                                            windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                                            scanMin=scanMin, scanMax=scanMaxList[i],
                                            histoFile=histoFile,
                                            histoString=theHistoString, 
                                            dataMode=theMode, verbosity=verb )
                    graphsList.append(ROOT.TGraphErrors())
                    for k in range(len(cutValue_list)):
                        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
                        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
                    graphsList[-1].SetLineColor(theColor)
                    graphsList[-1].SetName(theRegion+" "+theMode+" window=["+str(windowMinX)+","+str(windowMaxX)+"]x["+str(windowMinY)+","+str(windowMaxY)+"], xCut="+str(otherDirectionCutList[i]))
                    graphsList[-1].GetXaxis().SetTitle("ABCD y Cut")
                    graphsList[-1].GetYaxis().SetTitle("predicted/observed")
                    if theMode=="simQCD":
                        graphsList[-1].SetLineStyle(2)
                    else:
                        graphsList[-1].SetLineStyle(1)    
                    global_max=max(global_max,max(nonClosure_list))
                    global_min=min(global_min, min(nonClosure_list))

                canvas=ROOT.TCanvas("canvas_M15_1_xWindow","canvas_M15_1_xWindow",1024,768)
                global_max=min(2.5, global_max)
                print("global max", global_max)
                global_min=min(0.9, global_min)
                # global_min=0.8
                # print(graphsList)
                # print(graphsList[0])
                # graphsList[0].Print()
                graphsList[0].SetMaximum(global_max*1.1)
                graphsList[0].SetMinimum(global_min*0.9)
                graphsList[0].GetXaxis().SetLimits(scanMin-0.05, 1.0)
                graphsList[0].Draw("ALP0")
                for gr in graphsList[1:]:
                    gr.Draw("same LP0")
                canvas.BuildLegend()
                canvas.SetGrid()
                windowString="["+str(windowMinX)+","+str(windowMaxX)+"]x["+str(windowMinY)+","+str(windowMaxY)+"]"
                windowString=windowString.replace(".","p")
                canvas.SaveAs(outputDir+"scan_yDirection_xRestricted_"+windowString+"_xCut_"+str(otherDirectionCutList[i])+"_"+WP_tagName+".png")
                canvas.SaveAs(outputDir+"scan_yDirection_xRestricted_"+windowString+"_xCut_"+str(otherDirectionCutList[i])+"_"+WP_tagName+".pdf")




###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
# extend range for non-SR lines

# restrict in x and scan in x

# create toplevel dir
toplevelDirName="outputs_comparison_sideband_closure_allRegions_New_extended"
if not os.path.isdir(toplevelDirName):
    os.mkdir(toplevelDirName)

masses=["15","40","55"]
lifetimes=["0","0p05","0p1","1",]

histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
filename="all_histos_main_v6_invertedChi2.root"

do_VRs=False

for mass in masses:
    for lifetime in lifetimes:
        WP_tagName="M"+mass+"_"+lifetime
        print("")
        print("At "+WP_tagName)
        outputDir=toplevelDirName+"/"+"comparison_sideband_x-scan_"+WP_tagName+"/"
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)

        if mass=="15":
            histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG"
            histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_VR"
            histoString_SR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_invertedChi2"
            histoString_VR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_invertedChi2_VR"
        if mass=="40" or mass=="55":
            histoString_SR="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG"
            histoString_VR="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_VR"
            histoString_SR_invChi2="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_invertedChi2"
            histoString_VR_invChi2="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_invertedChi2_VR"
        histoString_SR=histoString_SR.replace("$WPTAG",WP_tagName)
        histoString_VR=histoString_VR.replace("$WPTAG",WP_tagName)
        histoString_SR_invChi2=histoString_SR_invChi2.replace("$WPTAG",WP_tagName)
        histoString_VR_invChi2=histoString_VR_invChi2.replace("$WPTAG",WP_tagName)
        
        lowerbound_list=[0.0, 0.3, 0.4, 0.5, 0.6,]
        for lowerbound in lowerbound_list:
            global_max=0.0
            global_min=1.0

            windowMinX=lowerbound
            windowMinY=lowerbound
            windowMaxY=1.0
            scanDirection="x"
            scanMin=lowerbound+0.01

            scanMaxList=[ 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69]
            otherDirectionCutList=[ 0.7,  0.83, 0.85, 0.87, 0.89, 0.91, 0.93, 0.95]
            windowMaxXList=[  0.7 for s in scanMaxList]

            for i in range(len(scanMaxList)):
                print(i)
                verb=1
                verb=0
                graphXOffset=0.0000*i
                graphsList=[]
                global_max=0.
                global_min=10.
                colorList=[   ROOT.kRed,  ROOT.kMagenta,  ]
                modes_list=["data", "data", ]
                region_list=["SR","SR inv Chi2", ]
                histoStringList=[histoString_SR, histoString_SR_invChi2,]
                if do_VRs:
                    colorList=[   ROOT.kRed,  ROOT.kMagenta, ROOT.kBlue, ROOT.kCyan,  ]
                    modes_list=["data", "data", "data", "data",]
                    region_list=["SR","SR inv Chi2", "VR", "VR inv Chi2" ]                
                    histoStringList=[histoString_SR, histoString_SR_invChi2, histoString_VR, histoString_VR_invChi2]

                for theRegion, theHistoString, theColor, theMode in zip(region_list, histoStringList, colorList, modes_list):
                    print("At", mass, lifetime, "loop counter", i, theRegion, theHistoString, theColor, theMode)
                    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                                            windowMinX=windowMinX, windowMaxX=windowMaxXList[i], windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                                            scanMin=scanMin, scanMax=scanMaxList[i],
                                            histoFile=histoFile,
                                            histoString=theHistoString, 
                                            dataMode=theMode, verbosity=verb )
                    graphsList.append(ROOT.TGraphErrors())
                    for k in range(len(cutValue_list)):
                        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
                        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
                    graphsList[-1].SetLineColor(theColor)
                    graphsList[-1].SetName(theRegion+" "+theMode+" window=["+str(windowMinX)+","+str(windowMaxXList[i])+"]x["+str(windowMinY)+","+str(windowMaxY)+"], yCut="+str(otherDirectionCutList[i]))
                    graphsList[-1].GetXaxis().SetTitle("ABCD x Cut")
                    graphsList[-1].GetYaxis().SetTitle("prediction/observed")
                    if theMode=="simQCD":
                        graphsList[-1].SetLineStyle(2)
                    else:
                        graphsList[-1].SetLineStyle(1)    
                    global_max=max(global_max,max(nonClosure_list))
                    global_min=min(global_min, min(nonClosure_list))

                # now do non-SR extension
                colorList=[   ROOT.kMagenta,  ]
                modes_list=["data", ]
                region_list=["SR inv Chi2", ]
                histoStringList=[histoString_SR_invChi2,]
                if do_VRs:
                    colorList=[    ROOT.kMagenta, ROOT.kBlue, ROOT.kCyan,  ]
                    modes_list=[ "data", "data", "data",]
                    region_list=["SR inv Chi2", "VR", "VR inv Chi2" ]                
                    histoStringList=[ histoString_SR_invChi2, histoString_VR, histoString_VR_invChi2]
                for theRegion, theHistoString, theColor, theMode in zip(region_list, histoStringList, colorList, modes_list):
                    print("At", mass, lifetime, "loop counter", i, theRegion, theHistoString, theColor, theMode)
                    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                                            windowMinX=windowMinX, windowMaxX=1.0, windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                                            scanMin=scanMin, scanMax=1.0,
                                            histoFile=histoFile,
                                            histoString=theHistoString, 
                                            dataMode=theMode, verbosity=verb )
                    graphsList.append(ROOT.TGraphErrors())
                    for k in range(len(cutValue_list)):
                        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
                        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
                    graphsList[-1].SetLineColor(theColor)
                    graphsList[-1].SetName(theRegion+" "+theMode+" window=["+str(windowMinX)+","+str(1.0)+"]x["+str(windowMinY)+","+str(windowMaxY)+"], yCut="+str(otherDirectionCutList[i]))
                    graphsList[-1].GetXaxis().SetTitle("ABCD x Cut")
                    graphsList[-1].GetYaxis().SetTitle("prediction/observed")
                    if theMode=="simQCD":
                        graphsList[-1].SetLineStyle(2)
                    else:
                        graphsList[-1].SetLineStyle(2)    
                    global_max=max(global_max,max(nonClosure_list))
                    global_min=min(global_min, min(nonClosure_list))


                canvas=ROOT.TCanvas("canvas_M15_1_xWindow","canvas_M15_1_xWindow",1024,768)
                global_max=min(2.5, global_max)
                global_min=min(0.9, global_min)
                # global_min=0.8
                # print(graphsList)
                # print(graphsList[0])
                # graphsList[0].Print()
                print("global max", global_max)
                graphsList[0].SetMaximum(global_max*1.1)
                graphsList[0].SetMinimum(global_min*0.9)
                graphsList[0].GetXaxis().SetLimits(scanMin-0.05, 1.0)
                graphsList[0].Draw("ALP0")
                for gr in graphsList[1:]:
                    gr.Draw("same LP0")
                canvas.BuildLegend()
                canvas.SetGrid()
                windowString="["+str(windowMinX)+","+str(windowMaxXList[i])+"]x["+str(windowMinY)+","+str(windowMaxY)+"]"
                windowString=windowString.replace(".","p")
                canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_"+windowString+"_yCut_"+str(otherDirectionCutList[i])+"_"+WP_tagName+".png")
                canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_"+windowString+"_yCut_"+str(otherDirectionCutList[i])+"_"+WP_tagName+".pdf")


# now restrict in y-direction and scan in y-direction

# create toplevel dir
toplevelDirName="outputs_comparison_sideband_closure_allRegions_New_extended"
if not os.path.isdir(toplevelDirName):
    os.mkdir(toplevelDirName)

masses=["15","40","55"]
lifetimes=["0","0p05","0p1","1",]

histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
filename="all_histos_main_v6_invertedChi2.root"

do_VRs=False

for mass in masses:
    for lifetime in lifetimes:
        WP_tagName="M"+mass+"_"+lifetime
        print("")
        print("At "+WP_tagName)
        outputDir=toplevelDirName+"/"+"comparison_sideband_y-scan_"+WP_tagName+"/"
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)

        if mass=="15":
            histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG"
            histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_VR"
            histoString_SR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_invertedChi2"
            histoString_VR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_invertedChi2_VR"
        if mass=="40" or mass=="55":
            histoString_SR="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG"
            histoString_VR="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_VR"
            histoString_SR_invChi2="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_invertedChi2"
            histoString_VR_invChi2="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_invertedChi2_VR"
        histoString_SR=histoString_SR.replace("$WPTAG",WP_tagName)
        histoString_VR=histoString_VR.replace("$WPTAG",WP_tagName)
        histoString_SR_invChi2=histoString_SR_invChi2.replace("$WPTAG",WP_tagName)
        histoString_VR_invChi2=histoString_VR_invChi2.replace("$WPTAG",WP_tagName)
        
        lowerbound_list=[0.0, 0.3, 0.4, 0.5, 0.6,]
        for lowerbound in lowerbound_list:
            global_max=0.0
            global_min=1.0

            windowMinX=lowerbound
            windowMinY=lowerbound
            windowMaxX=1.0
            scanDirection="y"
            scanMin=lowerbound+0.01

            scanMaxList=[ 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69]
            otherDirectionCutList=[ 0.7,  0.83, 0.85, 0.87, 0.89, 0.91, 0.93, 0.95]
            windowMaxYList=[  0.7 for s in scanMaxList]

            for i in range(len(scanMaxList)):
                print(i)
                verb=0
                graphXOffset=0.0000*i
                graphsList=[]
                global_max=0.
                global_min=10.
                colorList=[   ROOT.kRed,  ROOT.kMagenta,  ]
                modes_list=["data", "data", ]
                region_list=["SR","SR inv Chi2", ]
                histoStringList=[histoString_SR, histoString_SR_invChi2,]
                if do_VRs:
                    colorList=[   ROOT.kRed,  ROOT.kMagenta, ROOT.kBlue, ROOT.kCyan,  ]
                    modes_list=["data", "data", "data", "data",]
                    region_list=["SR","SR inv Chi2", "VR", "VR inv Chi2" ]                
                    histoStringList=[histoString_SR, histoString_SR_invChi2, histoString_VR, histoString_VR_invChi2]

                for theRegion, theHistoString, theColor, theMode in zip(region_list, histoStringList, colorList, modes_list):
                    print("At", mass, lifetime, "loop counter", i, theRegion, theHistoString, theColor, theMode)
                    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                                            windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinY, windowMaxY=windowMaxYList[i], otherDirectionCut=otherDirectionCutList[i],
                                            scanMin=scanMin, scanMax=scanMaxList[i],
                                            histoFile=histoFile,
                                            histoString=theHistoString, 
                                            dataMode=theMode, verbosity=verb )
                    graphsList.append(ROOT.TGraphErrors())
                    for k in range(len(cutValue_list)):
                        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
                        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
                    graphsList[-1].SetLineColor(theColor)
                    graphsList[-1].SetName(theRegion+" "+theMode+" window=["+str(windowMinX)+","+str(windowMaxX)+"]x["+str(windowMinY)+","+str(windowMaxYList[i])+"], xCut="+str(otherDirectionCutList[i]))
                    graphsList[-1].GetXaxis().SetTitle("ABCD y Cut")
                    graphsList[-1].GetYaxis().SetTitle("prediction/observed")
                    if theMode=="simQCD":
                        graphsList[-1].SetLineStyle(2)
                    else:
                        graphsList[-1].SetLineStyle(1)    
                    global_max=max(global_max,max(nonClosure_list))
                    global_min=min(global_min, min(nonClosure_list))

                # now do non-SR extension
                colorList=[   ROOT.kMagenta,  ]
                modes_list=["data", ]
                region_list=["SR inv Chi2", ]
                histoStringList=[histoString_SR_invChi2,]
                if do_VRs:
                    colorList=[    ROOT.kMagenta, ROOT.kBlue, ROOT.kCyan,  ]
                    modes_list=[ "data", "data", "data",]
                    region_list=["SR inv Chi2", "VR", "VR inv Chi2" ]                
                    histoStringList=[ histoString_SR_invChi2, histoString_VR, histoString_VR_invChi2]
                for theRegion, theHistoString, theColor, theMode in zip(region_list, histoStringList, colorList, modes_list):
                    print("At", mass, lifetime, "loop counter", i, theRegion, theHistoString, theColor, theMode)
                    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                                            windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinY, windowMaxY=1.0, otherDirectionCut=otherDirectionCutList[i],
                                            scanMin=scanMin, scanMax=1.0,
                                            histoFile=histoFile,
                                            histoString=theHistoString, 
                                            dataMode=theMode, verbosity=verb )
                    graphsList.append(ROOT.TGraphErrors())
                    for k in range(len(cutValue_list)):
                        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
                        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
                    graphsList[-1].SetLineColor(theColor)
                    graphsList[-1].SetName(theRegion+" "+theMode+" window=["+str(windowMinX)+","+str(windowMaxX)+"]x["+str(windowMinY)+","+str(1.0)+"], xCut="+str(otherDirectionCutList[i]))
                    graphsList[-1].GetXaxis().SetTitle("ABCD y Cut")
                    graphsList[-1].GetYaxis().SetTitle("prediction/observed")
                    if theMode=="simQCD":
                        graphsList[-1].SetLineStyle(2)
                    else:
                        graphsList[-1].SetLineStyle(2)    
                    global_max=max(global_max,max(nonClosure_list))
                    global_min=min(global_min, min(nonClosure_list))


                canvas=ROOT.TCanvas("canvas_M15_1_yWindow","canvas_M15_1_yWindow",1024,768)
                global_max=min(2.5, global_max)
                global_min=min(0.9, global_min)
                # global_min=0.8
                # print(graphsList)
                # print(graphsList[0])
                # graphsList[0].Print()
                print("global max", global_max)
                graphsList[0].SetMaximum(global_max*1.1)
                graphsList[0].SetMinimum(global_min*0.9)
                graphsList[0].GetXaxis().SetLimits(scanMin-0.05, 1.0)
                graphsList[0].Draw("ALP0")
                for gr in graphsList[1:]:
                    gr.Draw("same LP0")
                canvas.BuildLegend()
                canvas.SetGrid()
                windowString="["+str(windowMinX)+","+str(windowMaxX)+"]x["+str(windowMinY)+","+str(windowMaxYList[i])+"]"
                windowString=windowString.replace(".","p")
                canvas.SaveAs(outputDir+"scan_yDirection_yRestricted_"+windowString+"_xCut_"+str(otherDirectionCutList[i])+"_"+WP_tagName+".png")
                canvas.SaveAs(outputDir+"scan_yDirection_yRestricted_"+windowString+"_xCut_"+str(otherDirectionCutList[i])+"_"+WP_tagName+".pdf")



###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

# scan diagonal for SR_invertedChi2 for different upper restrictions
print("scan diagonal")
toplevelDirName="outputs_comparison_sideband_closure_allRegions_New_extended"
if not os.path.isdir(toplevelDirName):
    os.mkdir(toplevelDirName)

masses=["15","40","55"]
lifetimes=["0","0p05","0p1","1",]

histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
filename="all_histos_main_v6_invertedChi2.root"

do_VRs=False

for mass in masses:
    for lifetime in lifetimes:
        WP_tagName="M"+mass+"_"+lifetime
        print("")
        print("At "+WP_tagName)
        outputDir=toplevelDirName+"/"+"SRPrime_sideband_diagonal-scan_"+WP_tagName+"/"
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)

        if mass=="15":
            histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG"
            histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_VR"
            histoString_SR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_invertedChi2"
            histoString_VR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_$WPTAG_invertedChi2_VR"
        if mass=="40" or mass=="55":
            histoString_SR="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG"
            histoString_VR="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_VR"
            histoString_SR_invChi2="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_invertedChi2"
            histoString_VR_invChi2="resolved_Higgs_combi_chi2_A_three_phiSorted_12_$SAMPLE_sigprob_$WPTAG_invertedChi2_VR"
        histoString_SR=histoString_SR.replace("$WPTAG",WP_tagName)
        histoString_VR=histoString_VR.replace("$WPTAG",WP_tagName)
        histoString_SR_invChi2=histoString_SR_invChi2.replace("$WPTAG",WP_tagName)
        histoString_VR_invChi2=histoString_VR_invChi2.replace("$WPTAG",WP_tagName)
        
        lowerbound_list=[0.0, 0.3, 0.4, 0.5, 0.6,]
        for lowerbound in lowerbound_list:
            global_max=0.0
            global_min=1.0

            windowMinX=lowerbound
            windowMinY=lowerbound
            scanDirection="diag"
            scanMin=lowerbound+0.01

            windowMaxList=[0.7, 0.75, 0.8, 0.85, 0.9, 1.0]
            scanMaxList=[value-0.01 for value in windowMaxList]
            otherDirectionCutList=[ -1.0 for value in windowMaxList]
            colorList=[ROOT.kRed, ROOT.kMagenta, ROOT.kBlue, ROOT.kCyan, ROOT.kGreen, ROOT.kGreen+3]
            graphsList=[]
            global_max=0.
            global_min=10.
            theHistoString=histoString_SR_invChi2
            theRegion="SR inv Chi2"
            theMode="data"
            for i in range(len(scanMaxList)):
                print(i)
                verb=0
                graphXOffset=0.0000*i
                print("At", mass, lifetime, "loop counter", i, theRegion, theHistoString, theMode, windowMaxList[i])
                cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                                        windowMinX=windowMinX, windowMaxX=windowMaxList[i], windowMinY=windowMinY, windowMaxY=windowMaxList[i], otherDirectionCut=otherDirectionCutList[i],
                                        scanMin=scanMin, scanMax=scanMaxList[i],
                                        histoFile=histoFile,
                                        histoString=theHistoString, 
                                        dataMode=theMode, verbosity=verb )
                graphsList.append(ROOT.TGraphErrors())
                for k in range(len(cutValue_list)):
                    graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
                    graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
                graphsList[-1].SetLineColor(colorList[i])
                graphsList[-1].SetName(theRegion+" "+theMode+" window=["+str(windowMinX)+","+str(windowMaxList[i])+"]x["+str(windowMinY)+","+str(windowMaxList[i])+"]")
                graphsList[-1].GetXaxis().SetTitle("ABCD x,y Cut")
                graphsList[-1].GetYaxis().SetTitle("prediction/observed")
                if theMode=="simQCD":
                    graphsList[-1].SetLineStyle(2)
                else:
                    graphsList[-1].SetLineStyle(1)    
                global_max=max(global_max,max(nonClosure_list))
                global_min=min(global_min, min(nonClosure_list))

            canvas=ROOT.TCanvas("canvas_M15_1_xyWindow","canvas_M15_1_xyWindow",1024,768)
            global_max=min(2.5, global_max)
            global_min=min(0.9, global_min)
            # global_min=0.8
            # print(graphsList)
            # print(graphsList[0])
            # graphsList[0].Print()
            print("global max", global_max)
            graphsList[0].SetMaximum(global_max*1.1)
            graphsList[0].SetMinimum(global_min*0.9)
            graphsList[0].GetXaxis().SetLimits(scanMin-0.05, 1.0)
            graphsList[0].Draw("ALP0")
            for gr in graphsList[1:]:
                gr.Draw("same LP0")
            canvas.BuildLegend()
            canvas.SetGrid()
            windowString="["+str(windowMinX)+","+str(windowMaxList[i])+"]x["+str(windowMinY)+","+str(windowMaxList[i])+"]"
            windowString=windowString.replace(".","p")
            canvas.SaveAs(outputDir+"scan_diagonal_Restricted_"+windowString+"_"+WP_tagName+".png")
            canvas.SaveAs(outputDir+"scan_diagonal_Restricted_"+windowString+"_"+WP_tagName+".pdf")
