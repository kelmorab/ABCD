import ROOT
import math
import ctypes
import numpy
ROOT.gROOT.SetBatch(True)


class ABCD:
    # CA
    # DB
    def __init__(self):
        self.binTolerance=0.001
        self.histo=None
        self.histo_data=None
        self.histo_nonQCD=None
        self.histo_QCD=None
        self.regions_list=["A","B","C","D"]
        self.windowsAreSet=[False for region in self.regions_list]
        self.cut_value_dict={}
        for region in self.regions_list:
            for coord in ["X","Y"]:
                for i in ["1","2"]:
                    self.cut_value_dict[region+coord+i]=0.0

    def setDataHisto(self, histo):
        self.histo_data=histo

    def setNonQCDHisto(self, histo):
        self.histo_nonQCD=histo

    def setQCDHisto(self, histo):
        self.histo_QCD=histo
        self.setActiveHisto(self.histo_QCD)

    def setActiveHisto(self, histo):
        self.histo=histo

    def subtractNonQCDFromData(self):
        self.histo_QCD=self.histo_data.Clone()
        self.histo_QCD.Add(self.histo_nonQCD, -1.0)
        # print(self.histo_QCD.Integral(), self.histo_data.Integral(), self.histo_nonQCD.Integral())
        self.setActiveHisto(self.histo_QCD)

    def findBin(self, x,y):
        theBinX=self.histo.GetXaxis().FindBin(x)
        theBinY=self.histo.GetYaxis().FindBin(y)
        return theBinX, theBinY        

    def printWindow(self, region="A"):
        if region not in self.regions_list:
            print("WARNING: Only know ABCD. Doing Nothing!")
            return False
        s=region+": "    
        s2=""
        for coord in ["X","Y"]:
            for i in ["1","2"]:
                s+=coord+i+" = "+str(self.cut_value_dict[region+coord+i])+" , "
        # print("This corresponds to the following bins")
        s2+=" "+str(self.findBin(self.cut_value_dict[region+"X1"],self.cut_value_dict[region+"Y1"]))
        s2+=" "+str(self.findBin(self.cut_value_dict[region+"X2"],self.cut_value_dict[region+"Y1"]))
        s2+=" "+str(self.findBin(self.cut_value_dict[region+"X1"],self.cut_value_dict[region+"Y2"]))
        s2+=" "+str(self.findBin(self.cut_value_dict[region+"X2"],self.cut_value_dict[region+"Y2"]))
        s+=" --> "+s2
        print(s)

    # def checkForOverlap(self):
    # TODO

    def setWindow(self,region="A", X1=0.0, X2=0.0, Y1=0.0, Y2=0.0, verbose=0):
        defaultReturn = None
        if region not in self.regions_list:
            print("WARNING: Only know ABCD. Doing Nothing!")
            return defaultReturn
        vals=[X1, X2, Y1, Y2]
        if verbose>=2:
            print(vals)
        if X1==X2 or Y1==Y2:
            if verbose>=1:
                print("WARNING: Areas can not be 0. Check your window values. Doing Nothing")
                print(vals)
            return defaultReturn    
        if verbose>=2:
            print("Setting up region "+region)
        X1Bin, Y1Bin= self.findBin(X1, Y1)
        X2Bin, Y2Bin= self.findBin(X2, Y2)
        actualX1=self.histo.GetXaxis().GetBinLowEdge(X1Bin)
        self.cut_value_dict[region+"X1"]=actualX1
        if actualX1!=X1 and verbose>=2:
            print("X1 does not fit on bin edge. Setting to lower bin edge ->"+str(actualX1))
            print(X1, actualX1)
        if X1Bin==X2Bin:
            X2Bin=X1Bin+1
            if verbose>=2:
                print("Both X values are in same bin! Moving x2 up to end of next bin")
        actualX2=self.histo.GetXaxis().GetBinUpEdge(X2Bin)-self.histo.GetXaxis().GetBinWidth(X2Bin)*self.binTolerance
        if verbose>=2:
            print(self.histo.GetXaxis().GetBinUpEdge(X2Bin), self.histo.GetXaxis().GetBinWidth(X2Bin), self.histo.GetXaxis().GetBinWidth(X2Bin)*self.binTolerance )
        self.cut_value_dict[region+"X2"]=actualX2
        if actualX2!=X2 and verbose>=2:
            print("X2 does not fit on bin edge. Setting to slightly below upper bin edge ->"+str(actualX2))

        actualY1=self.histo.GetYaxis().GetBinLowEdge(Y1Bin)
        self.cut_value_dict[region+"Y1"]=actualY1
        if actualY1!=Y1 and verbose>=2:
            print("Y1 does not fit on bin edge. Setting to lower bin edge ->"+str(actualY1))
        if Y1Bin==Y2Bin:
            Y2Bin=Y1Bin+1
            if verbose>=2:
                print("Both Y values are in same bin! Moving Y2 up to end of neYt bin")
        actualY2=self.histo.GetYaxis().GetBinUpEdge(Y2Bin)-self.histo.GetYaxis().GetBinWidth(Y2Bin)*self.binTolerance
        self.cut_value_dict[region+"Y2"]=actualY2
        if actualY2!=Y2 and verbose>=2:
            print("Y2 does not fit on bin edge. Setting to slightly below upper bin edge ->"+str(actualY2))
        
        self.windowsAreSet[self.regions_list.index(region)]=True
        if verbose>=2:
            print("Set up region "+region)
            print("The actual window edges are:")
        if verbose>=1:
            self.printWindow(region)
        return True

    def setWindowFromBins(self, region="A", binX1=1, binX2=1, binY1=1, binY2=1, verbose=0):
        if verbose>=2:
            print("setting window from given bins:",region,  binX1, binX2, binY1, binY2)
        X1=self.histo.GetXaxis().GetBinLowEdge(binX1)
        X2=self.histo.GetXaxis().GetBinLowEdge(binX2)
        Y1=self.histo.GetYaxis().GetBinLowEdge(binY1)
        Y2=self.histo.GetYaxis().GetBinLowEdge(binY2)
        self.setWindow(region, X1, X2, Y1, Y2, verbose)
        return True
        

    def getIntegralAndErrorInRegion(self, region="A", histoName="QCD", verbose=0):
        if histoName=="QCD":
            self.setActiveHisto(self.histo_QCD)
        if histoName=="data":
            self.setActiveHisto(self.histo_data)
        if histoName=="nonQCD":
            self.setActiveHisto(self.histo_nonQCD)
            
        integral=0
        e=ctypes.c_double(0)
        x1=self.cut_value_dict[region+"X1"]
        x2=self.cut_value_dict[region+"X2"]
        y1=self.cut_value_dict[region+"Y1"]
        y2=self.cut_value_dict[region+"Y2"]
        x1bin, y1bin=self.findBin(x1,y1)
        x2bin, y2bin=self.findBin(x2,y2)
        
        integral=self.histo.IntegralAndError(x1bin, x2bin, y1bin, y2bin, e)
        if verbose>=2:
            print("integration over bins: ", x1bin, x2bin, y1bin, y2bin, " -> ", integral)
        error=e.value
        return integral, error

    def getIntegralAndErrorInWholeHisto(self, histoName="QCD", verbose=0):
        if histoName=="QCD":
            self.setActiveHisto(self.histo_QCD)
        if histoName=="data":
            self.setActiveHisto(self.histo_data)
        if histoName=="nonQCD":
            self.setActiveHisto(self.histo_nonQCD)
        integral=0
        e=ctypes.c_double(0)
        x1bin=1
        y1bin=1
        x2bin=self.histo.GetXaxis().GetNbins()
        y2bin=self.histo.GetYaxis().GetNbins()
        integral=self.histo.IntegralAndError(x1bin, x2bin, y1bin, y2bin, e)
        if verbose>=2:
            print("integration over bins: ", x1bin, x2bin, y1bin, y2bin, " -> ", integral)
        error=e.value
        return integral, error

    #TODO
    #def removeOverlap
    # move thresolds in direction of smallest error in integrals

    def predictRegionA(self, verbose=0):
        self.setActiveHisto(self.histo_QCD)
        pred_A=0
        pred_error_A=0
        B, BE =self.getIntegralAndErrorInRegion("B")
        C, CE =self.getIntegralAndErrorInRegion("C")
        D, DE =self.getIntegralAndErrorInRegion("D")
        
        if D<=0.0:
            if verbose>=1:
                print("WARNING: Integral in D<=0 = "+ str(D)+ " . Doing Nothing!")
            return pred_A, pred_error_A
        pred_A = B*C/D
        pred_error_A = math.sqrt( (C/D)**2 * BE**2 + (B/D)**2 * CE**2 + (B*C/D/D)**2 * DE**2   )
        return pred_A, pred_error_A

    def getNonClosure(self):
        pred_A, pred_error_A = self.predictRegionA()
        integral, error = self.getIntegralAndErrorInRegion("A")
        if integral<=0 or pred_A<=0:
            return 0.0
        closure=pred_A/integral
        return closure    

    def getNonClosureWrtData(self):
        pred_A, pred_error_A = self.predictRegionA()
        integral_nonQCD, error_nonQCD = self.getIntegralAndErrorInRegion("A", "nonQCD")
        integral_data, error_data = self.getIntegralAndErrorInRegion("A", "data")
        if integral_data<=0 or pred_A<=0:
            return 0.0
        pred_plus_nonQCD=pred_A + integral_nonQCD
        nonClosure = pred_plus_nonQCD / integral_data
        return nonClosure
        


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
            pritn(cutValue)
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


##############################################################################################################################################################################
histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1"
histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR"
histoString_SR_invertedChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2"
histoString_VR_invertedChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR"

histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")

histoString=histoString_VR
global_max=0.0
global_min=1.0

histo_data_obs_Chi2_SR=histoFile.Get(histoString_SR.replace("$SAMPLE","data_obs"))
histo_SM_Higgs_Chi2_SR = histoFile.Get(histoString_SR.replace("$SAMPLE","SM_Higgs"))
histo_TTbar_Chi2_SR = histoFile.Get(histoString_SR.replace("$SAMPLE","TTbar"))
histo_Others_Chi2_SR = histoFile.Get(histoString_SR.replace("$SAMPLE","Others"))
histo_nonQCD_Chi2_SR = histo_SM_Higgs_Chi2_SR.Clone()
histo_nonQCD_Chi2_SR.Add(histo_TTbar_Chi2_SR)
histo_nonQCD_Chi2_SR.Add(histo_Others_Chi2_SR)
theABCD_SR=ABCD()
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
theABCD_VR=ABCD()
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
theABCD_SR_invertedChi2=ABCD()
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
theABCD_VR_invertedChi2=ABCD()
theABCD_VR_invertedChi2.setDataHisto(histo_data_obs_Chi2_VR_invertedChi2)
theABCD_VR_invertedChi2.setNonQCDHisto(histo_nonQCD_Chi2_VR_invertedChi2)
theABCD_VR_invertedChi2.subtractNonQCDFromData()

cutX=0.85
cutY=0.85

theABCD_SR.setWindow("A", cutX, 1, cutY, 1)
theABCD_SR.setWindow("B", cutX, 1, 0, cutY-0.0001)
theABCD_SR.setWindow("C", 0, cutX-0.0001, cutY, 1)
theABCD_SR.setWindow("D", 0, cutX-0.0001, 0, cutY-0.0001)

theABCD_VR.setWindow("A", cutX, 1, cutY, 1)
theABCD_VR.setWindow("B", cutX, 1, 0, cutY-0.0001)
theABCD_VR.setWindow("C", 0, cutX-0.0001, cutY, 1)
theABCD_VR.setWindow("D", 0, cutX-0.0001, 0, cutY-0.0001)

theABCD_SR_invertedChi2.setWindow("A", cutX, 1, cutY, 1)
theABCD_SR_invertedChi2.setWindow("B", cutX, 1, 0, cutY-0.0001)
theABCD_SR_invertedChi2.setWindow("C", 0, cutX-0.0001, cutY, 1)
theABCD_SR_invertedChi2.setWindow("D", 0, cutX-0.0001, 0, cutY-0.0001)

theABCD_VR_invertedChi2.setWindow("A", cutX, 1, cutY, 1)
theABCD_VR_invertedChi2.setWindow("B", cutX, 1, 0, cutY-0.0001)
theABCD_VR_invertedChi2.setWindow("C", 0, cutX-0.0001, cutY, 1)
theABCD_VR_invertedChi2.setWindow("D", 0, cutX-0.0001, 0, cutY-0.0001)


print("theABCD_SR", theABCD_SR.getIntegralAndErrorInWholeHisto()[0])
print("theABCD_VR", theABCD_VR.getIntegralAndErrorInWholeHisto()[0])
print("theABCD_SR_invertedChi2", theABCD_SR_invertedChi2.getIntegralAndErrorInWholeHisto()[0])
print("theABCD_VR_invertedChi2", theABCD_VR_invertedChi2.getIntegralAndErrorInWholeHisto()[0])