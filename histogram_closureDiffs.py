import ROOT
import math
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptFit()


infileName="2022_03_14_table_closures_V2.csv"
valueList=[]

valueList_VR_VRbar=[]
valueList_VR_VRbar_X=[]
valueList_VR_VRbar_Y=[]
valueList_SR_SRbar_X=[]
valueList_SR_SRbar_Y=[]


with open(infileName) as infile:
    l=list(infile)
    for line in l:
        tableCounter=0
        if "mm" in line:
            tableCounter=0
        if "GeV" in line or "mm" in line or "," not in line:
            continue
        tableCounter+=1
        print(tableCounter)
        print(line)
        sl=line.split(";")
        print(sl)
        fsl=[]
        for frag in sl:
            if "," in frag:
                fsl.append(float(frag.replace(",",".")))
        print(fsl)
        valueList+=fsl
        if tableCounter==1:
            valueList_VR_VRbar+=fsl
        if tableCounter==2:
            valueList_VR_VRbar_X+=fsl
        if tableCounter==3:
            valueList_VR_VRbar_Y+=fsl
        if tableCounter==4:
            valueList_VR_VRbar_X+=fsl
        if tableCounter==5:
            valueList_VR_VRbar_Y+=fsl
            

print(len(valueList))
minValue=min(valueList)
maxValue=max(valueList)
nBins=int(len(valueList)/2.)
nBins=140
maxValue=10
minValue=-4
histoAll=ROOT.TH1D("closure_diff_sigs","all",int(nBins), minValue, maxValue)
for val in valueList:
    histoAll.Fill(val)
histoAll.GetXaxis().SetTitle("closure difference significance")
print(histoAll.GetBinContent(nBins+1))
histoAll.SetBinContent(nBins, histoAll.GetBinContent(nBins)+histoAll.GetBinContent(nBins+1))
histoAll.SetBinContent(1, histoAll.GetBinContent(0)+histoAll.GetBinContent(1))

canvas=ROOT.TCanvas("c_closure_diff_sigs","c_closure_diff_sigs",1024, 768)
histoAll.Draw("histo")
canvas.SaveAs("closureHistos/histo_closure_diff_sigs_signed_V1.pdf")

nBins=70
maxValue=10
minValue=-4
histoAll=ROOT.TH1D("closure_diff_sigs","all",int(nBins), minValue, maxValue)
for val in valueList:
    histoAll.Fill(val)
histoAll.GetXaxis().SetTitle("closure difference significance")
print(histoAll.GetBinContent(nBins+1))
histoAll.SetBinContent(nBins, histoAll.GetBinContent(nBins)+histoAll.GetBinContent(nBins+1))
histoAll.SetBinContent(1, histoAll.GetBinContent(0)+histoAll.GetBinContent(1))

canvas=ROOT.TCanvas("c_closure_diff_sigs","c_closure_diff_sigs",1024, 768)
histoAll.Draw("histo")
canvas.SaveAs("closureHistos/histo_closure_diff_sigs_signed_V2.pdf")

nBins=170
maxValue=30
minValue=-4
histoAll=ROOT.TH1D("closure_diff_sigs","all",int(nBins), minValue, maxValue)
for val in valueList:
    histoAll.Fill(val)
histoAll.GetXaxis().SetTitle("closure difference significance")
print(histoAll.GetBinContent(nBins+1))
histoAll.SetBinContent(nBins, histoAll.GetBinContent(nBins)+histoAll.GetBinContent(nBins+1))
histoAll.SetBinContent(1, histoAll.GetBinContent(0)+histoAll.GetBinContent(1))

canvas=ROOT.TCanvas("c_closure_diff_sigs","c_closure_diff_sigs",1024, 768)
histoAll.Draw("histo")
canvas.SaveAs("closureHistos/histo_closure_diff_sigs_signed_V3.pdf")

nBins=70
maxValue=10
minValue=-4
histoAll=ROOT.TH1D("closure_diff_sigs","all",int(nBins), minValue, maxValue)
for val in valueList:
    histoAll.Fill(val)
histoAll.GetXaxis().SetTitle("closure difference significance")
print(histoAll.GetBinContent(nBins+1))
histoAll.SetBinContent(nBins, histoAll.GetBinContent(nBins)+histoAll.GetBinContent(nBins+1))
histoAll.SetBinContent(1, histoAll.GetBinContent(0)+histoAll.GetBinContent(1))

canvas=ROOT.TCanvas("c_closure_diff_sigs","c_closure_diff_sigs",1024, 768)
histoAll.Draw("histo")
func=ROOT.TF1("f1","gaus",-4,4)
histoAll.Fit("gaus","","",-4,4)
func=histoAll.GetFunction("gaus")
func.Draw("same")
# func.Draw("same")
canvas.SaveAs("closureHistos/histo_closure_diff_sigs_signed_onlyPeak_V2.pdf")



# labelList=["VR-VRbar","VR_VRbar_X","VR_VRbar_Y","SR_SRbar_X","SR_SRbar_Y"]
# histoList=[histo_VR_VRbar, histo_VR_VRbar_X, histo_VR_VRbar_Y, histo_SR_SRbar_X, histo_SR_SRbar_Y]
# listlist=[valueList_VR_VRbar, valueList_VR_VRbar_X, valueList_VR_VRbar_Y, valueList_SR_SRbar_X, valueList_SR_SRbar_Y]

# histoAll_VR_VRbar=ROOT.TH1D("closure_diff_sigs_VR_VRbar","VR-VR_bar",int(nBins), minValue, maxValue)
# for val in valueList_VR_VRbar:
#     histoAll_VR_VRbar.Fill(val)
# histoAll_VR_VRbar.GetXaxis().SetTitle("closure difference significance")
# print(histoAll_VR_VRbar.GetBinContent(nBins+1))
# histoAll_VR_VRbar.SetBinContent(nBins, histoAll_VR_VRbar.GetBinContent(nBins)+histoAll_VR_VRbar.GetBinContent(nBins+1))

