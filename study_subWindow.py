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
            closure, closure_error=theABCD.getNonClosureWrtData()
            if verbosity>0: print("closure: "+str(round(closure,3)))
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
            if verbosity>0: print("non-closure: "+str(round(closure,3)))
            cutValue_list.append(cutValue)
            nonClosure_list.append(closure)  
            non_Closure_error_list.append(closure_error)          

    return cutValue_list, nonClosure_list, non_Closure_error_list


##############################################################################################################################################################################
histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1"
histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR"
histoFile=ROOT.TFile("all_histos_datacards_main_v6.root","r")
outputDir="correct_scanOneDirection_M15_1/"

histoString=histoString_VR
global_max=0.0
global_min=1.0

windowMinX=0.0 
windowMaxX=1.
windowMinY=0.0
windowMaxY=1.0
scanDirection="x"

scanMaxList=[0.9, 0.84, 0.69, 0.9, 0.84, 0.69, 0.69, 0.69, 0.84, 0.9,   ]
otherDirectionCutList=[0.85, 0.85, 0.7, 0.85, 0.85, 0.7, 0.7, 0.7, 0.85, 0.85  ]
colorList=[ROOT.kBlack, ROOT.kBlue, ROOT.kCyan, ROOT.kBlack, ROOT.kBlue, ROOT.kCyan, ROOT.kRed, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1 ]
modes_list=["data", "data", "data", "simQCD", "simQCD", "simQCD", "data", "simQCD", "simQCD", "simQCD",   ]
region_list=["VR","VR", "VR", "VR","VR", "VR", "SR", "SR", "SR", "SR",   ]
windowMaxXList=[1.0, 0.85, 0.7, 1.0, 0.85, 0.7, 0.7, 0.7, 0.85, 1.0,   ]
windowMinXList=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   ]
graphsList=[]
assert(len(scanMaxList)==len(colorList))

graphXOffset=0.005*0

for i in range(len(scanMaxList)):
    print(i)
    if region_list[i]=="VR":
        histoString=histoString_VR
    elif region_list[i]=="SR":
        histoString=histoString_SR

    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                             windowMinX=windowMinXList[i], windowMaxX=windowMaxXList[i], windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                             scanMin=0.5, scanMax=scanMaxList[i],
                             histoFile=histoFile, histoString=histoString,
                             dataMode=modes_list[i], verbosity=0)
    graphsList.append(ROOT.TGraphErrors())
    for k in range(len(cutValue_list)):
        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
    graphsList[-1].SetLineColor(colorList[i])
    graphsList[-1].SetName(region_list[i]+" "+modes_list[i]+" xRange=["+str(windowMinXList[i])+","+str(windowMaxXList[i])+"] yCut="+str(otherDirectionCutList[i])+" xCut=[0.5,"+str(scanMaxList[i])+"]")
    graphsList[-1].GetXaxis().SetTitle("ABCD x Cut")
    graphsList[-1].GetYaxis().SetTitle("closure")
    if modes_list[i]=="simQCD":
        graphsList[-1].SetLineStyle(2)
    else:
        graphsList[-1].SetLineStyle(1)    
    global_max=max(global_max,max(nonClosure_list))
    global_min=min(global_min, min(nonClosure_list))
canvas=ROOT.TCanvas("canvas_M15_ctau1_xWindow","canvas_M15_ctau1_xWindow",1024,768)
global_max=min(2.5, global_max)
global_min=0
graphsList[0].SetMaximum(global_max*1.1)
graphsList[0].SetMinimum(global_min*0.9)
graphsList[0].Draw("ALP0")
for gr in graphsList[1:]:
    gr.Draw("same LP0")

canvas.BuildLegend()
canvas.SetGrid()
canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_M15_1.png")
canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_M15_1.pdf")


# test
# histo_data_obs=histoFile.Get(histoString.replace("$SAMPLE","data_obs"))
# histo_SM_Higgs = histoFile.Get(histoString.replace("$SAMPLE","SM_Higgs"))
# histo_TTbar = histoFile.Get(histoString.replace("$SAMPLE","TTbar"))
# histo_Others = histoFile.Get(histoString.replace("$SAMPLE","Others"))
# histo_nonQCD = histo_SM_Higgs.Clone()
# histo_nonQCD.Add(histo_TTbar)
# histo_nonQCD.Add(histo_Others)
# theABCD=ABCD()
# theABCD.setDataHisto(histo_data_obs)
# theABCD.setNonQCDHisto(histo_nonQCD)
# theABCD.subtractNonQCDFromData()
# windowMaxX=1.0
# windowMaxY=1.0
# windowMinX=0.0
# windowMinY=0.0
# cutValue=0.6
# otherDirectionCut=0.7
# theABCD.setWindow("A", cutValue, windowMaxX, otherDirectionCut, windowMaxY)
# theABCD.setWindow("B", cutValue, windowMaxX, windowMinY, otherDirectionCut-0.0001)
# theABCD.setWindow("C", windowMinX, cutValue-0.0001, otherDirectionCut, windowMaxY)
# theABCD.setWindow("D", windowMinX, cutValue-0.0001, windowMinY, otherDirectionCut-0.0001)

# exit(0)
#########################################################################################################################################################################################################################
#########################################################################################################################################################################################################################
#########################################################################################################################################################################################################################

# Now limit the window to 0.7 (0.85) in y-direction and repeat the scan in x-diretion

histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1"
histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR"
histoFile=ROOT.TFile("all_histos_datacards_main_v6.root","r")
outputDir="correct_scanOneDirection_M15_1/"

histoString=histoString_VR
global_max=0.0
global_min=1.0

windowMinX=0.0 
windowMaxX=1.
windowMinY=0.0
windowMaxY=1.0
scanDirection="x"

scanMaxList=[0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95,    0.95, 0.95, 0.95, 0.95 , 0.95,  ]
otherDirectionCutList=[0.85, 0.7, 0.7, 0.6, 0.85, 0.7, 0.7, 0.6,   0.6, 0.6, 0.7, 0.85, 0.7   ]
colorList=[ROOT.kBlack, ROOT.kSpring, ROOT.kBlue, ROOT.kCyan, ROOT.kBlack, ROOT.kSpring, ROOT.kBlue, ROOT.kCyan,    ROOT.kRed, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1 , ROOT.kMagenta-5  ]
modes_list=["data", "data", "data", "data", "simQCD", "simQCD", "simQCD", "simQCD",    "data", "simQCD", "simQCD","simQCD", "simQCD",    ]
region_list=["VR", "VR", "VR", "VR", "VR", "VR", "VR", "VR",    "SR", "SR", "SR", "SR", "SR"   ]
windowMaxYList=[1.0, 1.0, 0.85, 0.7, 1.0, 1.0, 0.85, 0.7,    0.7, 0.7, 0.85, 1.0, 1.0,   ]
windowMinYList=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   ]
graphsList=[]

graphXOffset=0.005*0

for i in range(len(scanMaxList)):
    if region_list[i]=="VR":
        histoString=histoString_VR
    elif region_list[i]=="SR":
        histoString=histoString_SR

    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                             windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinYList[i], windowMaxY=windowMaxYList[i], otherDirectionCut=otherDirectionCutList[i],
                             scanMin=0.5, scanMax=scanMaxList[i],
                             histoFile=histoFile, histoString=histoString,
                             dataMode=modes_list[i], verbosity=0 )
    graphsList.append(ROOT.TGraphErrors())
    for k in range(len(cutValue_list)):
        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
    graphsList[-1].SetLineColor(colorList[i])
    graphsList[-1].SetName(region_list[i]+" "+modes_list[i]+" yRange=["+str(windowMinYList[i])+","+str(windowMaxYList[i])+"] yCut="+str(otherDirectionCutList[i])+" xCut=[0.5,"+str(scanMaxList[i])+"]")
    graphsList[-1].GetXaxis().SetTitle("ABCD x-Cut")
    graphsList[-1].GetYaxis().SetTitle("closure")
    if modes_list[i]=="simQCD":
        graphsList[-1].SetLineStyle(2)
    else:
        graphsList[-1].SetLineStyle(1)    
    global_max=max(global_max,max(nonClosure_list))
    global_min=min(global_min, min(nonClosure_list))

canvas=ROOT.TCanvas("canvas_M15_ctau1_yWindow","canvas_M15_ctau1_yWindow",1024,768)
global_max=min(2.5, global_max)
global_min=0
graphsList[0].SetMaximum(global_max*1.1)
graphsList[0].SetMinimum(global_min*0.9)
graphsList[0].Draw("ALP0")
for gr in graphsList[1:]:
    gr.Draw("same LP0")

canvas.BuildLegend()
canvas.SetGrid()
canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_M15_1.png")
canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_M15_1.pdf")

#########################################################################################################################################################################################################################
#########################################################################################################################################################################################################################
##########################################################################################################################################################################################################################
# invertedChi2

histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2"
histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR"
histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
outputDir="correct_scanOneDirection_M15_1/"

histoString=histoString_VR
global_max=0.0
global_min=1.0

windowMinX=0.0 
windowMaxX=1.
windowMinY=0.0
windowMaxY=1.0
scanDirection="x"

scanMaxList=[0.9, 0.84, 0.69, 0.9, 0.84, 0.69, 0.69, 0.84, 0.9, 0.69, 0.84, 0.9,   ]
otherDirectionCutList=[0.85, 0.85, 0.7, 0.85, 0.85, 0.7, 0.7, 0.85, 0.85, 0.7, 0.85, 0.85  ]
colorList=[ROOT.kBlack, ROOT.kBlue, ROOT.kCyan, ROOT.kBlack, ROOT.kBlue, ROOT.kCyan, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1 ]
modes_list=["data", "data", "data", "simQCD", "simQCD", "simQCD", "data", "data", "data", "simQCD", "simQCD", "simQCD",   ]
region_list=["VR","VR", "VR", "VR","VR", "VR", "SR", "SR", "SR", "SR", "SR", "SR",   ]
windowMaxXList=[1.0, 0.85, 0.7, 1.0, 0.85, 0.7, 0.7, 0.85, 1.0, 0.7, 0.85, 1.0,   ]
windowMinXList=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   ]
graphsList=[]
assert(len(scanMaxList)==len(colorList))

graphXOffset=0.005*0

for i in range(len(scanMaxList)):
    print(i)
    if region_list[i]=="VR":
        histoString=histoString_VR
    elif region_list[i]=="SR":
        histoString=histoString_SR

    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                             windowMinX=windowMinXList[i], windowMaxX=windowMaxXList[i], windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                             scanMin=0.5, scanMax=scanMaxList[i],
                             histoFile=histoFile, histoString=histoString,
                             dataMode=modes_list[i], verbosity=0 )
    graphsList.append(ROOT.TGraphErrors())
    for k in range(len(cutValue_list)):
        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
    graphsList[-1].SetLineColor(colorList[i])
    graphsList[-1].SetName(region_list[i]+" "+modes_list[i]+" xRange=["+str(windowMinXList[i])+","+str(windowMaxXList[i])+"] yCut="+str(otherDirectionCutList[i])+" xCut=[0.5,"+str(scanMaxList[i])+"]")
    graphsList[-1].GetXaxis().SetTitle("ABCD x Cut")
    graphsList[-1].GetYaxis().SetTitle("closure")
    if modes_list[i]=="simQCD":
        graphsList[-1].SetLineStyle(2)
    else:
        graphsList[-1].SetLineStyle(1)    
    global_max=max(global_max,max(nonClosure_list))
    global_min=min(global_min, min(nonClosure_list))

canvas=ROOT.TCanvas("canvas_M15_ctau1_xWindow_invertedChi2","canvas_M15_ctau1_xWindow_invertedChi2",1024,768)
global_max=min(2.5, global_max)
global_min=0
graphsList[0].SetMaximum(global_max*1.1)
graphsList[0].SetMinimum(global_min*0.9)
graphsList[0].Draw("ALP0")
for gr in graphsList[1:]:
    gr.Draw("same LP0")

canvas.BuildLegend()
canvas.SetGrid()
canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_invertedChi2_M15_1.png")
canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_invertedChi2_M15_1.pdf")


# Now limit the window to 0.7 (0.85) in y-direction and repeat the scan in x-diretion

histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2"
histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR"
histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
outputDir="correct_scanOneDirection_M15_1/"

histoString=histoString_VR
global_max=0.0
global_min=1.0

windowMinX=0.0 
windowMaxX=1.
windowMinY=0.0
windowMaxY=1.0
scanDirection="x"

scanMaxList=[0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95,    0.95, 0.95, 0.95, 0.95 , 0.95, 0.95, 0.95 , 0.95,    ]
otherDirectionCutList=[0.85, 0.7, 0.7, 0.6, 0.85, 0.7, 0.7, 0.6,   0.6, 0.6, 0.7, 0.85, 0.7, 0.7, 0.85, 0.7    ]
colorList=[ROOT.kBlack, ROOT.kSpring, ROOT.kBlue, ROOT.kCyan, ROOT.kBlack, ROOT.kSpring, ROOT.kBlue, ROOT.kCyan,    ROOT.kRed, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1 , ROOT.kMagenta-5, ROOT.kOrange, ROOT.kMagenta+1 , ROOT.kMagenta-5   ]
modes_list=["data", "data", "data", "data", "simQCD", "simQCD", "simQCD", "simQCD",    "data", "simQCD", "simQCD","simQCD", "simQCD",  "data","data","data",  ]
region_list=["VR", "VR", "VR", "VR", "VR", "VR", "VR", "VR",    "SR", "SR", "SR", "SR", "SR" ,"SR", "SR", "SR"   ]
windowMaxYList=[1.0, 1.0, 0.85, 0.7, 1.0, 1.0, 0.85, 0.7,    0.7, 0.7, 0.85, 1.0, 1.0 , 0.85, 1.0, 1.0,  ]
windowMinYList=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  ]
graphsList=[]

graphXOffset=0.005*0

for i in range(len(scanMaxList)):
    if region_list[i]=="VR":
        histoString=histoString_VR
    elif region_list[i]=="SR":
        histoString=histoString_SR

    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                             windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinYList[i], windowMaxY=windowMaxYList[i], otherDirectionCut=otherDirectionCutList[i],
                             scanMin=0.5, scanMax=scanMaxList[i],
                             histoFile=histoFile, histoString=histoString,
                             dataMode=modes_list[i], verbosity=0 )
    graphsList.append(ROOT.TGraphErrors())
    for k in range(len(cutValue_list)):
        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
    graphsList[-1].SetLineColor(colorList[i])
    graphsList[-1].SetName(region_list[i]+" "+modes_list[i]+" yRange=["+str(windowMinYList[i])+","+str(windowMaxYList[i])+"] yCut="+str(otherDirectionCutList[i])+" xCut=[0.5,"+str(scanMaxList[i])+"]")
    graphsList[-1].GetXaxis().SetTitle("ABCD x-Cut")
    graphsList[-1].GetYaxis().SetTitle("closure")
    if modes_list[i]=="simQCD":
        graphsList[-1].SetLineStyle(2)
    else:
        graphsList[-1].SetLineStyle(1)    
    global_max=max(global_max,max(nonClosure_list))
    global_min=min(global_min, min(nonClosure_list))

canvas=ROOT.TCanvas("canvas_M15_ctau1_yWindow_invertedChi2","canvas_M15_ctau1_yWindow_invertedChi2",1024,768)
global_max=min(2.5, global_max)
global_min=0
graphsList[0].SetMaximum(global_max*1.1)
graphsList[0].SetMinimum(global_min*0.9)
graphsList[0].Draw("ALP0")
for gr in graphsList[1:]:
    gr.Draw("same LP0")

canvas.BuildLegend()
canvas.SetGrid()
canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_invertedChi2_M15_1.png")
canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_invertedChi2_M15_1.pdf")

#########################################################################################################################################################################################################################
#########################################################################################################################################################################################################################
##########################################################################################################################################################################################################################
#########################################################################################################################################################################################################################
#########################################################################################################################################################################################################################
##########################################################################################################################################################################################################################

# Now everything again but with lower bound of 0.3 in both directions


histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1"
histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR"
histoFile=ROOT.TFile("all_histos_datacards_main_v6.root","r")
outputDir="correct_scanOneDirection_M15_1/"

histoString=histoString_VR
global_max=0.0
global_min=1.0

windowMinX=0.3
windowMaxX=1.
windowMinY=0.3
windowMaxY=1.0
scanDirection="x"

scanMaxList=[0.9, 0.84, 0.69, 0.9, 0.84, 0.69, 0.69, 0.69, 0.84, 0.9,   ]
otherDirectionCutList=[0.85, 0.85, 0.7, 0.85, 0.85, 0.7, 0.7, 0.7, 0.85, 0.85  ]
colorList=[ROOT.kBlack, ROOT.kBlue, ROOT.kCyan, ROOT.kBlack, ROOT.kBlue, ROOT.kCyan, ROOT.kRed, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1 ]
modes_list=["data", "data", "data", "simQCD", "simQCD", "simQCD", "data", "simQCD", "simQCD", "simQCD",   ]
region_list=["VR","VR", "VR", "VR","VR", "VR", "SR", "SR", "SR", "SR",   ]
windowMaxXList=[1.0, 0.85, 0.7, 1.0, 0.85, 0.7, 0.7, 0.7, 0.85, 1.0,   ]
windowMinXList=[0.3 ]*10
graphsList=[]
assert(len(scanMaxList)==len(colorList))

graphXOffset=0.005*0

for i in range(len(scanMaxList)):
    print(i)
    if region_list[i]=="VR":
        histoString=histoString_VR
    elif region_list[i]=="SR":
        histoString=histoString_SR

    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                             windowMinX=windowMinXList[i], windowMaxX=windowMaxXList[i], windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                             scanMin=0.5, scanMax=scanMaxList[i],
                             histoFile=histoFile, histoString=histoString,
                             dataMode=modes_list[i], verbosity=0 )
    graphsList.append(ROOT.TGraphErrors())
    for k in range(len(cutValue_list)):
        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
    graphsList[-1].SetLineColor(colorList[i])
    graphsList[-1].SetName(region_list[i]+" "+modes_list[i]+" xRange=["+str(windowMinXList[i])+","+str(windowMaxXList[i])+"] yCut="+str(otherDirectionCutList[i])+" xCut=[0.5,"+str(scanMaxList[i])+"]")
    graphsList[-1].GetXaxis().SetTitle("ABCD x Cut")
    graphsList[-1].GetYaxis().SetTitle("closure")
    if modes_list[i]=="simQCD":
        graphsList[-1].SetLineStyle(2)
    else:
        graphsList[-1].SetLineStyle(1)    
    global_max=max(global_max,max(nonClosure_list))
    global_min=min(global_min, min(nonClosure_list))

canvas=ROOT.TCanvas("canvas_M15_ctau1_xWindow","canvas_M15_ctau1_xWindow",1024,768)
global_max=min(2.5, global_max)
global_min=0
graphsList[0].SetMaximum(global_max*1.1)
graphsList[0].SetMinimum(global_min*0.9)
graphsList[0].Draw("ALP0")
for gr in graphsList[1:]:
    gr.Draw("same LP0")

canvas.BuildLegend()
canvas.SetGrid()
canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_lowRestriction_0p3_M15_1.png")
canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_lowRestriction_0p3_M15_1.pdf")


#########################################################################################################################################################################################################################
#########################################################################################################################################################################################################################
#########################################################################################################################################################################################################################

# Now limit the window to 0.7 (0.85) in y-direction and repeat the scan in x-diretion

histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1"
histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR"
histoFile=ROOT.TFile("all_histos_datacards_main_v6.root","r")
outputDir="correct_scanOneDirection_M15_1/"

histoString=histoString_VR
global_max=0.0
global_min=1.0

windowMinX=0.3 
windowMaxX=1.
windowMinY=0.3
windowMaxY=1.0
scanDirection="x"

scanMaxList=[0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95,    0.95, 0.95, 0.95, 0.95 , 0.95,  ]
otherDirectionCutList=[0.85, 0.7, 0.7, 0.6, 0.85, 0.7, 0.7, 0.6,   0.6, 0.6, 0.7, 0.85, 0.7   ]
colorList=[ROOT.kBlack, ROOT.kSpring, ROOT.kBlue, ROOT.kCyan, ROOT.kBlack, ROOT.kSpring, ROOT.kBlue, ROOT.kCyan,    ROOT.kRed, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1 , ROOT.kMagenta-5  ]
modes_list=["data", "data", "data", "data", "simQCD", "simQCD", "simQCD", "simQCD",    "data", "simQCD", "simQCD","simQCD", "simQCD",    ]
region_list=["VR", "VR", "VR", "VR", "VR", "VR", "VR", "VR",    "SR", "SR", "SR", "SR", "SR"   ]
windowMaxYList=[1.0, 1.0, 0.85, 0.7, 1.0, 1.0, 0.85, 0.7,    0.7, 0.7, 0.85, 1.0, 1.0,   ]
windowMinYList=[0.3   ]*13
graphsList=[]

graphXOffset=0.005*0

for i in range(len(scanMaxList)):
    if region_list[i]=="VR":
        histoString=histoString_VR
    elif region_list[i]=="SR":
        histoString=histoString_SR

    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                             windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinYList[i], windowMaxY=windowMaxYList[i], otherDirectionCut=otherDirectionCutList[i],
                             scanMin=0.5, scanMax=scanMaxList[i],
                             histoFile=histoFile, histoString=histoString,
                             dataMode=modes_list[i], verbosity=0 )
    graphsList.append(ROOT.TGraphErrors())
    for k in range(len(cutValue_list)):
        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
    graphsList[-1].SetLineColor(colorList[i])
    graphsList[-1].SetName(region_list[i]+" "+modes_list[i]+" yRange=["+str(windowMinYList[i])+","+str(windowMaxYList[i])+"] yCut="+str(otherDirectionCutList[i])+" xCut=[0.5,"+str(scanMaxList[i])+"]")
    graphsList[-1].GetXaxis().SetTitle("ABCD x-Cut")
    graphsList[-1].GetYaxis().SetTitle("closure")
    if modes_list[i]=="simQCD":
        graphsList[-1].SetLineStyle(2)
    else:
        graphsList[-1].SetLineStyle(1)    
    global_max=max(global_max,max(nonClosure_list))
    global_min=min(global_min, min(nonClosure_list))

canvas=ROOT.TCanvas("canvas_M15_ctau1_yWindow","canvas_M15_ctau1_yWindow",1024,768)
global_max=min(2.5, global_max)
global_min=0
graphsList[0].SetMaximum(global_max*1.1)
graphsList[0].SetMinimum(global_min*0.9)
graphsList[0].Draw("ALP0")
for gr in graphsList[1:]:
    gr.Draw("same LP0")

canvas.BuildLegend()
canvas.SetGrid()
canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_lowRestriction_0p3_M15_1.png")
canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_lowRestriction_0p3_M15_1.pdf")

#########################################################################################################################################################################################################################
#########################################################################################################################################################################################################################
##########################################################################################################################################################################################################################
# invertedChi2

histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2"
histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR"
histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
outputDir="correct_scanOneDirection_M15_1/"

histoString=histoString_VR
global_max=0.0
global_min=1.0

windowMinX=0.3 
windowMaxX=1.
windowMinY=0.3
windowMaxY=1.0
scanDirection="x"

scanMaxList=[0.9, 0.84, 0.69, 0.9, 0.84, 0.69, 0.69, 0.84, 0.9, 0.69, 0.84, 0.9,   ]
otherDirectionCutList=[0.85, 0.85, 0.7, 0.85, 0.85, 0.7, 0.7, 0.85, 0.85, 0.7, 0.85, 0.85  ]
colorList=[ROOT.kBlack, ROOT.kBlue, ROOT.kCyan, ROOT.kBlack, ROOT.kBlue, ROOT.kCyan, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1 ]
modes_list=["data", "data", "data", "simQCD", "simQCD", "simQCD", "data", "data", "data", "simQCD", "simQCD", "simQCD",   ]
region_list=["VR","VR", "VR", "VR","VR", "VR", "SR", "SR", "SR", "SR", "SR", "SR",   ]
windowMaxXList=[1.0, 0.85, 0.7, 1.0, 0.85, 0.7, 0.7, 0.85, 1.0, 0.7, 0.85, 1.0,   ]
windowMinXList=[0.3  ]*12
graphsList=[]
assert(len(scanMaxList)==len(colorList))

graphXOffset=0.005*0

for i in range(len(scanMaxList)):
    print(i)
    if region_list[i]=="VR":
        histoString=histoString_VR
    elif region_list[i]=="SR":
        histoString=histoString_SR

    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                             windowMinX=windowMinXList[i], windowMaxX=windowMaxXList[i], windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                             scanMin=0.5, scanMax=scanMaxList[i],
                             histoFile=histoFile, histoString=histoString,
                             dataMode=modes_list[i], verbosity=0 )
    graphsList.append(ROOT.TGraphErrors())
    for k in range(len(cutValue_list)):
        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
    graphsList[-1].SetLineColor(colorList[i])
    graphsList[-1].SetName(region_list[i]+" "+modes_list[i]+" xRange=["+str(windowMinXList[i])+","+str(windowMaxXList[i])+"] yCut="+str(otherDirectionCutList[i])+" xCut=[0.5,"+str(scanMaxList[i])+"]")
    graphsList[-1].GetXaxis().SetTitle("ABCD x Cut")
    graphsList[-1].GetYaxis().SetTitle("closure")
    if modes_list[i]=="simQCD":
        graphsList[-1].SetLineStyle(2)
    else:
        graphsList[-1].SetLineStyle(1)    
    global_max=max(global_max,max(nonClosure_list))
    global_min=min(global_min, min(nonClosure_list))

canvas=ROOT.TCanvas("canvas_M15_ctau1_xWindow_invertedChi2","canvas_M15_ctau1_xWindow_invertedChi2",1024,768)
global_max=min(2.5, global_max)
global_min=0
graphsList[0].SetMaximum(global_max*1.1)
graphsList[0].SetMinimum(global_min*0.9)
graphsList[0].Draw("ALP0")
for gr in graphsList[1:]:
    gr.Draw("same LP0")

canvas.BuildLegend()
canvas.SetGrid()
canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_lowRestriction_0p3_invertedChi2_M15_1.png")
canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_lowRestriction_0p3_invertedChi2_M15_1.pdf")


# Now limit the window to 0.7 (0.85) in y-direction and repeat the scan in x-diretion

histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2"
histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR"
histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
outputDir="correct_scanOneDirection_M15_1/"

histoString=histoString_VR
global_max=0.0
global_min=1.0

windowMinX=0.3
windowMaxX=1.
windowMinY=0.3
windowMaxY=1.0
scanDirection="x"

scanMaxList=[0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95,    0.95, 0.95, 0.95, 0.95 , 0.95, 0.95, 0.95 , 0.95,    ]
otherDirectionCutList=[0.85, 0.7, 0.7, 0.6, 0.85, 0.7, 0.7, 0.6,   0.6, 0.6, 0.7, 0.85, 0.7, 0.7, 0.85, 0.7    ]
colorList=[ROOT.kBlack, ROOT.kSpring, ROOT.kBlue, ROOT.kCyan, ROOT.kBlack, ROOT.kSpring, ROOT.kBlue, ROOT.kCyan,    ROOT.kRed, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1 , ROOT.kMagenta-5, ROOT.kOrange, ROOT.kMagenta+1 , ROOT.kMagenta-5   ]
modes_list=["data", "data", "data", "data", "simQCD", "simQCD", "simQCD", "simQCD",    "data", "simQCD", "simQCD","simQCD", "simQCD",  "data","data","data",  ]
region_list=["VR", "VR", "VR", "VR", "VR", "VR", "VR", "VR",    "SR", "SR", "SR", "SR", "SR" ,"SR", "SR", "SR"   ]
windowMaxYList=[1.0, 1.0, 0.85, 0.7, 1.0, 1.0, 0.85, 0.7,    0.7, 0.7, 0.85, 1.0, 1.0 , 0.85, 1.0, 1.0,  ]
windowMinYList=[0.3 ]*16
graphsList=[]

graphXOffset=0.005*0

for i in range(len(scanMaxList)):
    if region_list[i]=="VR":
        histoString=histoString_VR
    elif region_list[i]=="SR":
        histoString=histoString_SR

    cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection(scanDirection=scanDirection,
                             windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinYList[i], windowMaxY=windowMaxYList[i], otherDirectionCut=otherDirectionCutList[i],
                             scanMin=0.5, scanMax=scanMaxList[i],
                             histoFile=histoFile, histoString=histoString,
                             dataMode=modes_list[i], verbosity=0 )
    graphsList.append(ROOT.TGraphErrors())
    for k in range(len(cutValue_list)):
        graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
        graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
    graphsList[-1].SetLineColor(colorList[i])
    graphsList[-1].SetName(region_list[i]+" "+modes_list[i]+" yRange=["+str(windowMinYList[i])+","+str(windowMaxYList[i])+"] yCut="+str(otherDirectionCutList[i])+" xCut=[0.5,"+str(scanMaxList[i])+"]")
    graphsList[-1].GetXaxis().SetTitle("ABCD x-Cut")
    graphsList[-1].GetYaxis().SetTitle("closure")
    if modes_list[i]=="simQCD":
        graphsList[-1].SetLineStyle(2)
    else:
        graphsList[-1].SetLineStyle(1)    
    global_max=max(global_max,max(nonClosure_list))
    global_min=min(global_min, min(nonClosure_list))

canvas=ROOT.TCanvas("canvas_M15_ctau1_yWindow_invertedChi2","canvas_M15_ctau1_yWindow_invertedChi2",1024,768)
global_max=min(2.5, global_max)
global_min=0
graphsList[0].SetMaximum(global_max*1.1)
graphsList[0].SetMinimum(global_min*0.9)
graphsList[0].Draw("ALP0")
for gr in graphsList[1:]:
    gr.Draw("same LP0")

canvas.BuildLegend()
canvas.SetGrid()
canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_lowRestriction_0p3_invertedChi2_M15_1.png")
canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_lowRestriction_0p3_invertedChi2_M15_1.pdf")