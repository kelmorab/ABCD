import ROOT
import math
import ctypes

histoString="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_0"
histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_0_VR"

histoFile=ROOT.TFile("all_histos_datacards_main_v6.root","r")


class ABCD:
    # CA
    # DB
    def __init__(self, histo=None):
        self.binTolerance=0.001
        self.histo=histo
        self.regions_list=["A","B","C","D"]
        self.windowsAreSet=[False for region in self.regions_list]
        self.cut_value_dict={}
        for region in self.regions_list:
            for coord in ["X","Y"]:
                for i in ["1","2"]:
                    self.cut_value_dict[region+coord+i]=0.0

    def findBin(self, x,y):
        theBinX=self.histo.GetXaxis().FindBin(x)
        theBinY=self.histo.GetYaxis().FindBin(y)
        return theBinX, theBinY        

    def printWindow(self, region="A"):
        if region not in self.regions_list:
            print("WARNING: Only know ABCD. Doing Nothing!")
            return False
        s=region+": "    
        for coord in ["X","Y"]:
            for i in ["1","2"]:
                s+=coord+i+" = "+str(self.cut_value_dict[region+coord+i])+" , "
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
        

    def getIntegralAndErrorInRegion(self, region="A", verbose=0):
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

    def getIntegralAndErrorInWholeHisto(self, verbose=0):
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

    def getClosure(self):
        pred_A, pred_error_A = self.predictRegionA()
        integral, error = self.getIntegralAndErrorInRegion("A")
        if integral<=0 or pred_A<=0:
            return 0.0
        closure=integral/pred_A
        return closure    


    # TODO
    # correlations
    # profiles
    # reduced histograms -> correlations there


class ABCD_Analyzer:
    def __init__(self):
        self.histo=None
        self.theABCD=None

    def setHisto(self, histo):
        self.histo=histo
        self.theABCD=ABCD(self.histo)

    def scanClosuresForBinwisePartitionsOfWidth(self, windowWidth=1):
        closures_list=[]
        iComb=0
        nBinsX=self.histo.GetNbinsX()
        nBinsY=self.histo.GetNbinsY()
        for iBinCX1 in range(1,nBinsX+1):
            if iBinCX1+2*windowWidth>nBinsX:
                break
            iBinCX2=iBinCX1+windowWidth
            for iBinAX1 in range(iBinCX2+1, nBinsX+1):
                if iBinAX1+windowWidth>nBinsX:
                    break
                iBinAX2=iBinAX1+windowWidth

                for iBinDY1 in range(1, nBinsY+1):
                    if iBinDY1+2*windowWidth>nBinsY:
                        break
                    iBinDY2=iBinDY1+windowWidth
                    for iBinCY1 in range(iBinDY2+1, nBinsY+1):
                        if iBinCY1+windowWidth>nBinsY:
                            break
                        iBinCY2=iBinCY1+windowWidth

                        iBinDX1=iBinCX1
                        iBinDX2=iBinCX2
                        iBinAY1=iBinCY1
                        iBinAY2=iBinCY2
                        iBinBY1=iBinDY1
                        iBinBY2=iBinDY2
                        iBinBX1=iBinAX1
                        iBinBX2=iBinAX2
                        self.theABCD.setWindowFromBins("A",iBinAX1, iBinAX2, iBinAY1, iBinAY2)
                        self.theABCD.setWindowFromBins("B",iBinBX1, iBinBX2, iBinBY1, iBinBY2)
                        self.theABCD.setWindowFromBins("C",iBinCX1, iBinCX2, iBinCY1, iBinCY2)
                        self.theABCD.setWindowFromBins("D",iBinDX1, iBinDX2, iBinDY1, iBinDY2)
                        closure=self.theABCD.getClosure()
                        closures_list.append(closure)
                        iComb+=1
                        if iComb%100==0:
                            print("------------- at combination"+str(iComb)+" --------------")

        return closures_list         



histo_data_obs_VR=histoFile.Get(histoString_VR.replace("$SAMPLE","data_obs"))
theHisto=histo_data_obs_VR

theABCD=ABCD(theHisto)
theABCD.setWindow("A", 0.83, 1.0, 0.83, 1.0)
theABCD.setWindow("B", 0.83, 1.0, 0.0, 0.829999)
theABCD.setWindow("C", 0, 0.82999, 0.83, 1.0)
theABCD.setWindow("D", 0, 0.82999, 0, 0.82999)
print("")

print(theABCD.getIntegralAndErrorInWholeHisto()[0])
print(theABCD.getIntegralAndErrorInRegion("A")[0]+theABCD.getIntegralAndErrorInRegion("B")[0]+theABCD.getIntegralAndErrorInRegion("C")[0]+theABCD.getIntegralAndErrorInRegion("D")[0])
print("")
print("ABCD in these regions")
obs_A, obs_error_A = theABCD.getIntegralAndErrorInRegion("A")
print("A:", obs_A, obs_error_A)
pred_A, pred_error_A = theABCD.predictRegionA()
print("predicted A: ",pred_A, pred_error_A )
print("non-closure: "+str(obs_A/pred_A))

print("")
theABCDAnalyzer=ABCD_Analyzer()
theABCDAnalyzer.setHisto(theHisto)
closures_list=theABCDAnalyzer.scanClosuresForBinwisePartitionsOfWidth(40)