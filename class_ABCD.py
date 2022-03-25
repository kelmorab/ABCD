import ROOT
import math
import ctypes
import numpy

def calcMinMaxAllowedObs(significance=1, pred=1, pred_error=0.1, verbosity=0):
    # for now assume obs_error = sqrt(obs)
    alpha=-1*(2*pred+significance**2)
    beta=(pred**2 - significance**2 * pred_error**2)
    x1 = -alpha/2. + math.sqrt((alpha/2.0)**2 - beta)
    x2 = -alpha/2. - math.sqrt((alpha/2.0)**2 - beta)
    return x1, x2

def calcMinMaxAllowedClosure(significance=1, pred=1, pred_error=0.1, verbosity=0):
    x1, x2 = calcMinMaxAllowedObs(significance=significance, pred=pred, pred_error=pred_error, verbosity=verbosity)
    closure1=pred/float(x1)
    closure2=pred/float(x2)
    return closure1, closure2
    
def calcNonClosureSig(pred=1, pred_error=0.1, obs=1, obs_error=1):
    res=(pred-obs)/math.sqrt(pred_error**2 + obs_error**2)
    return res


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

    def emptyBinsCrossWise(self, xValue=0.76, yValue=0.76):
        for histo in [self.histo_data, self.histo_QCD, self.histo_nonQCD]:
            if histo==None:
                continue
            xBin=histo.GetXaxis().FindBin(xValue)
            yBin=histo.GetXaxis().FindBin(yValue)
            a,b,c = ctypes.c_int(0),  ctypes.c_int(0),  ctypes.c_int(0)
            nBinsX=histo.GetNbinsX()
            nBinsY=histo.GetNbinsY()
            nBins=(nBinsX+2)*(nBinsY+2)
            for theBin in range(nBins):
                histo.GetBinXYZ(theBin, a,b,c)
                if histo.GetXaxis().FindBin(xValue)==a.value or histo.GetYaxis().FindBin(yValue)==b.value:
                    histo.SetBinContent(theBin, 0.0)

    def printHistos(self, path=""):
        for histo, name in zip([self.histo_data, self.histo_QCD, self.histo_nonQCD], ["data","QCD","nonQCD"]):
            if histo==None:
                continue
            canvas=ROOT.TCanvas("c","c",1024,768)
            histo.Draw("colz")
            canvas.SetLogz()
            canvas.SaveAs(path+"_"+name+".pdf")
            
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
        pred_A = B*C/float(D)
        pred_error_A = math.sqrt( (C/D)**2 * BE**2 + (B/D)**2 * CE**2 + (B*C/D/D)**2 * DE**2   )
        return pred_A, pred_error_A

    def getNonClosure(self):
        pred_A, pred_error_A = self.predictRegionA()
        integral, error = self.getIntegralAndErrorInRegion("A")
        if integral<=0 or pred_A<=0:
            return 0.0, 0.0
        closure=pred_A/float(integral)
        closure_error=math.sqrt(pred_error_A**2 / integral**2 + pred_A**2 * error**2 / integral**4)
        return closure, closure_error

    def getNonClosureWrtData(self):
        pred_A, pred_error_A = self.predictRegionA()
        integral_nonQCD, error_nonQCD = self.getIntegralAndErrorInRegion("A", "nonQCD")
        integral_data, error_data = self.getIntegralAndErrorInRegion("A", "data")
        if integral_data<=0 or pred_A<=0:
            return 0.0, 0.0
        pred_plus_nonQCD=pred_A + integral_nonQCD
        pred_plus_nonQCD_error = math.sqrt(pred_error_A**2 + error_nonQCD**2 )
        closure = pred_plus_nonQCD / integral_data
        closure_error = math.sqrt(pred_plus_nonQCD_error**2 / integral_data**2 + pred_plus_nonQCD**2 * error_data**2 / integral_data**4)
        return closure, closure_error
        
    def getCorrectionFactorWrtData(self, verbose=False):
        pred_A, pred_error_A = self.predictRegionA()
        integral_nonQCD, error_nonQCD = self.getIntegralAndErrorInRegion("A", "nonQCD")
        integral_data, error_data = self.getIntegralAndErrorInRegion("A", "data")
        if verbose:
            print("getCorrectionFactorWrtData: predA: ", pred_A)
            print("getCorrectionFactorWrtData: integral_nonQCD: ", integral_nonQCD)
            print("getCorrectionFactorWrtData: integral_data: ", integral_data)
        if integral_data<=0 or pred_A<=0:
            return 0.0, 0.0
        data_minus_nonQCD = integral_data - integral_nonQCD
        correctionFactor = abs(data_minus_nonQCD / float(pred_A))
        error_correctionFactor = math.sqrt( error_data**2 /(pred_A**2) + error_nonQCD**2/(pred_A**2) + pred_error_A**2 * ((integral_data - integral_nonQCD)/(pred_A**2))**2  )
        return correctionFactor, error_correctionFactor

    def getCorrectionFactorWrtSimQCD(self):
        pred_A, pred_error_A = self.predictRegionA()
        integral, error = self.getIntegralAndErrorInRegion("A")
        if integral<=0 or pred_A<=0:
            return 0.0, 0.0
        correctionFactor = abs(integral / float(pred_A))
        error_correctionFactor = math.sqrt( error**2 /(pred_A**2) + pred_error_A**2 * integral**2/pred_A**4  )
        return correctionFactor, error_correctionFactor

    def predictRegionAWithCorrectionAndAddUncertainty(self, correctionFactor=1., correctionFactor_unc=0.0, add_unc=0.0, verbose=0):
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
        pred_A = B*C/float(D)
        pred_error_A = math.sqrt( (C/D)**2 * BE**2 + (B/D)**2 * CE**2 + (B*C/D/D)**2 * DE**2   )
        corrected_pred_A = correctionFactor * pred_A
        total_error_pred_A = math.sqrt(pred_A**2 * (correctionFactor_unc**2 ) + correctionFactor**2 * (pred_error_A**2) )
        # this is an absolute uncertainty. But the add_unc. should act as relative uncertainty
        total_error_pred_A_relative=0.0
        if corrected_pred_A>0:
            total_error_pred_A_relative=total_error_pred_A/float(corrected_pred_A)
        total_plusAdd_error_pred_A_relative=math.sqrt(total_error_pred_A_relative**2 + add_unc**2)    
        # now revert back to absolute uncertainty
        total_plus_Add_error_pred_A = total_plusAdd_error_pred_A_relative * corrected_pred_A
        # stuff
        correctionFactor_unc_relative=0.0
        if correctionFactor>0:
            correctionFactor_unc_relative=correctionFactor_unc/correctionFactor
        pred_error_A_relative=0.0
        if pred_A>0:
            pred_error_A_relative = pred_error_A/pred_A
        if verbose:
            print("predictRegionAWithCorrectionAndAddUncertainty: corrected_pred_A, total_error_pred_A", corrected_pred_A, total_plus_Add_error_pred_A)
            print("predictRegionAWithCorrectionAndAddUncertainty: pred_A, pred_error_A, relative pred_error_A: "+str(round(pred_A,3))+" +- "+str(round(pred_error_A,3))+" ("+str(round(pred_error_A_relative,6))+")")
            print("predictRegionAWithCorrectionAndAddUncertainty: correctionFactor, correctionFactor_unc, relative correctionFactor_unc: "+str(round(correctionFactor,3))+" +- "+str(round(correctionFactor_unc,3))+" ("+str(round(correctionFactor_unc_relative,6))+")")
            print("predictRegionAWithCorrectionAndAddUncertainty: additional uncertainty:"+str(add_unc))
        return corrected_pred_A, total_plus_Add_error_pred_A    

    def getClosureFactorWrtDataWithCorrectionAndAddUncertainty(self, correctionFactor=1., correctionFactor_unc=0.0, add_unc=0.0, verbose=0):
        pred_A, pred_error_A = self.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor, correctionFactor_unc, add_unc, verbose)
        integral_nonQCD, error_nonQCD = self.getIntegralAndErrorInRegion("A", "nonQCD")
        integral_data, error_data = self.getIntegralAndErrorInRegion("A", "data")
        if integral_data<=0 or pred_A<=0:
            return 0.0, 0.0
        pred_plus_nonQCD=pred_A + integral_nonQCD
        pred_plus_nonQCD_error = math.sqrt(pred_error_A**2 + error_nonQCD**2 )
        closure = integral_data / pred_plus_nonQCD
        closure_error = math.sqrt(error_data**2 / pred_plus_nonQCD**2 + integral_data**2 * pred_plus_nonQCD_error**2 / pred_plus_nonQCD**4)
        return closure, closure_error

    def getClosureFactorWrtSimQCDWithCorrectionAndAddUncertainty(self, correctionFactor=1., correctionFactor_unc=0.0, add_unc=0.0, verbose=0):
        pred_A, pred_error_A = self.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor, correctionFactor_unc, add_unc, verbose)
        integral, error = self.getIntegralAndErrorInRegion("A")
        if integral<=0 or pred_A<=0:
            return 0.0, 0.0
        closure = integral / pred_A
        closure_error = math.sqrt(error**2 / pred_A**2 + integral**2 * pred_error_A**2 / pred_A**4)
        if verbose:
            print("getClosureFactorWrtSimQCDWithCorrectionAndAddUncertainty: pred_A, pred_error_A: ",pred_A, pred_error_A)
            print("getClosureFactorWrtSimQCDWithCorrectionAndAddUncertainty: closure, closure_error: ",closure, closure_error)
        return closure, closure_error

    def getNonClosureWrtSimQCDWithCorrectionAndAddUncertainty(self, correctionFactor=1., correctionFactor_unc=0.0, add_unc=0.0, verbose=0):
        pred_A, pred_error_A = self.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor, correctionFactor_unc, add_unc, verbose)
        integral, error = self.getIntegralAndErrorInRegion("A")
        if integral<=0 or pred_A<=0:
            return 0.0, 0.0
        closure=pred_A/float(integral)
        closure_error=math.sqrt(pred_error_A**2 / integral**2 + pred_A**2 * error**2 / integral**4)
        return closure, closure_error

    def getNonClosureWrtDataWithCorrectionAndAddUncertainty(self, correctionFactor=1., correctionFactor_unc=0.0, add_unc=0.0, verbose=0):
        pred_A, pred_error_A = self.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor, correctionFactor_unc, add_unc, verbose)
        integral_nonQCD, error_nonQCD = self.getIntegralAndErrorInRegion("A", "nonQCD")
        integral_data, error_data = self.getIntegralAndErrorInRegion("A", "data")
        if integral_data<=0 or pred_A<=0:
            return 0.0, 0.0
        pred_plus_nonQCD=pred_A + integral_nonQCD
        pred_plus_nonQCD_error = math.sqrt(pred_error_A**2 + error_nonQCD**2 )
        closure = pred_plus_nonQCD / integral_data
        closure_error = math.sqrt(pred_plus_nonQCD_error**2 / integral_data**2 + pred_plus_nonQCD**2 * error_data**2 / integral_data**4)
        if verbose:
            print("getClosureFactorWrtSimQCDWithCorrectionAndAddUncertainty: pred_A, pred_error_A: ",pred_A, pred_error_A)
            print("getClosureFactorWrtSimQCDWithCorrectionAndAddUncertainty: closure, closure_error: ",closure, closure_error)
        return closure, closure_error

    def getSimpleClosureSignificance(self, dataMode="data", verbosity=0, doUnsigned=True):
        if dataMode=="data":
            pred_A, pred_error_A = self.predictRegionA()
            integral_nonQCD, error_nonQCD = self.getIntegralAndErrorInRegion("A", "nonQCD")
            integral_data, error_data = self.getIntegralAndErrorInRegion("A", "data")
            if integral_data<=0 or pred_A<=0:
                return 0.0
            pred_plus_nonQCD=pred_A + integral_nonQCD
            pred_plus_nonQCD_error = math.sqrt(pred_error_A**2 + error_nonQCD**2 )
            diff=pred_plus_nonQCD-integral_data
            total_error=math.sqrt(pred_plus_nonQCD_error**2 + error_data**2 )
            if total_error<=0:
                return 0.0
            significance=diff/float(total_error)
            if doUnsigned:
                significance=abs(significance)
            if verbosity:
                print("diff, total_error, significance:", diff, total_error, significance)
            return significance
        elif dataMode=="simQCD":
            pred_A, pred_error_A = self.predictRegionA()
            integral_data, error_data = self.getIntegralAndErrorInRegion("A")
            if integral_data<=0 or pred_A<=0:
                return 0.0
            diff=pred_A-integral_data
            total_error=math.sqrt(pred_error_A**2 + error_data**2 )
            if total_error<=0:
                return 0.0
            significance=diff/float(total_error)
            if doUnsigned:
                significance=abs(significance)
            if verbosity:
                print("diff, total_error, significance:", diff, total_error, significance)
            return significance

    def getMaxAllowedNonClosures(self, sig=1, dataMode="data", verbosity=0):
        if dataMode=="data":
            pred_A, pred_error_A = self.predictRegionA()
            integral_nonQCD, error_nonQCD = self.getIntegralAndErrorInRegion("A", "nonQCD")
            integral_data, error_data = self.getIntegralAndErrorInRegion("A", "data")
            if verbosity>0:
                print("data, data_error:"+str(integral_data)+" "+str(error_data))
            if integral_data<=0 or pred_A<=0:
                return 0.0, 0.0
            pred_plus_nonQCD=pred_A + integral_nonQCD
            pred_plus_nonQCD_error = math.sqrt(pred_error_A**2 + error_nonQCD**2 )
            if pred_plus_nonQCD_error<=0:
                return 0.0, 0.0
            closure1, closure2 = calcMinMaxAllowedClosure(significance=sig, pred=pred_plus_nonQCD, pred_error=pred_plus_nonQCD_error, verbosity=verbosity)
            if verbosity:
                print("maximally allowed nonclosures: ",closure1, closure2)
            return abs(1.0 - closure1), abs(1.0 - closure2)
        elif dataMode=="simQCD":
            pred_A, pred_error_A = self.predictRegionA()
            integral_data, error_data = self.getIntegralAndErrorInRegion("A")
            if integral_data<=0 or pred_A<=0:
                return 0.0, 0.0
            if pred_error_A<=0:
                return 0.0, 0.0
            closure1, closure2 = calcMinMaxAllowedClosure(significance=sig, pred=pred_A, pred_error=pred_error_A, verbosity=verbosity)
            if verbosity:
                print("maximally allowed nonclosures: ",closure1, closure2)
            return abs(1.0 - closure1), abs(1.0 - closure2)
class FullABCD:
    def __init__(self):
        self.theABCD_SR=None
        self.theABCD_VR=None
        self.theABCD_SR_prime=None
        self.theABCD_VR_prime=None
        self.theABCD_SR_file=None
        self.theABCD_VR_file=None
        self.theABCD_SR_prime_file=None
        self.theABCD_VR_prime_file=None
        

    def setup_ABCD(self, region="", filename="", histoString="", dataMode="data"):
        if region=="SR":
            self.theABCD_SR_file=ROOT.TFile(filename,"r")
            self.theABCD_SR=None
            if dataMode=="data":
                histo_data_obs=self.theABCD_SR_file.Get(histoString.replace("$SAMPLE","data_obs"))
                histo_SM_Higgs = self.theABCD_SR_file.Get(histoString.replace("$SAMPLE","SM_Higgs"))
                histo_TTbar = self.theABCD_SR_file.Get(histoString.replace("$SAMPLE","TTbar"))
                histo_Others = self.theABCD_SR_file.Get(histoString.replace("$SAMPLE","Others"))
                histo_nonQCD = histo_SM_Higgs.Clone()
                histo_nonQCD.Add(histo_TTbar)
                histo_nonQCD.Add(histo_Others)
                self.theABCD_SR=ABCD()
                self.theABCD_SR.setDataHisto(histo_data_obs)
                self.theABCD_SR.setNonQCDHisto(histo_nonQCD)
                self.theABCD_SR.subtractNonQCDFromData()
            elif dataMode=="simQCD":
                histo_QCD=self.theABCD_SR_file.Get(histoString.replace("$SAMPLE","QCD"))
                self.theABCD_SR=ABCD()
                self.theABCD_SR.setQCDHisto(histo_QCD)
        if region=="VR":
            self.theABCD_VR_file=ROOT.TFile(filename,"r")
            self.theABCD_VR=None
            if dataMode=="data":
                histo_data_obs=self.theABCD_VR_file.Get(histoString.replace("$SAMPLE","data_obs"))
                histo_SM_Higgs = self.theABCD_VR_file.Get(histoString.replace("$SAMPLE","SM_Higgs"))
                histo_TTbar = self.theABCD_VR_file.Get(histoString.replace("$SAMPLE","TTbar"))
                histo_Others = self.theABCD_VR_file.Get(histoString.replace("$SAMPLE","Others"))
                histo_nonQCD = histo_SM_Higgs.Clone()
                histo_nonQCD.Add(histo_TTbar)
                histo_nonQCD.Add(histo_Others)
                self.theABCD_VR=ABCD()
                self.theABCD_VR.setDataHisto(histo_data_obs)
                self.theABCD_VR.setNonQCDHisto(histo_nonQCD)
                self.theABCD_VR.subtractNonQCDFromData()
            elif dataMode=="simQCD":
                histo_QCD=self.theABCD_VR_file.Get(histoString.replace("$SAMPLE","QCD"))
                self.theABCD_VR=ABCD()
                self.theABCD_VR.setQCDHisto(histo_QCD)
        if region=="SR_prime":
            self.theABCD_SR_prime_file=ROOT.TFile(filename,"r")
            self.theABCD_SR_prime=None
            if dataMode=="data":
                histo_data_obs=self.theABCD_SR_prime_file.Get(histoString.replace("$SAMPLE","data_obs"))
                histo_SM_Higgs = self.theABCD_SR_prime_file.Get(histoString.replace("$SAMPLE","SM_Higgs"))
                histo_TTbar = self.theABCD_SR_prime_file.Get(histoString.replace("$SAMPLE","TTbar"))
                histo_Others = self.theABCD_SR_prime_file.Get(histoString.replace("$SAMPLE","Others"))
                histo_nonQCD = histo_SM_Higgs.Clone()
                histo_nonQCD.Add(histo_TTbar)
                histo_nonQCD.Add(histo_Others)
                self.theABCD_SR_prime=ABCD()
                self.theABCD_SR_prime.setDataHisto(histo_data_obs)
                self.theABCD_SR_prime.setNonQCDHisto(histo_nonQCD)
                self.theABCD_SR_prime.subtractNonQCDFromData()
            elif dataMode=="simQCD":
                histo_QCD=self.theABCD_SR_prime_file.Get(histoString.replace("$SAMPLE","QCD"))
                self.theABCD_SR_prime=ABCD()
                self.theABCD_SR_prime.setQCDHisto(histo_QCD)
        if region=="VR_prime":
            self.theABCD_VR_prime_file=ROOT.TFile(filename,"r")
            self.theABCD_VR_prime=None
            if dataMode=="data":
                histo_data_obs=self.theABCD_VR_prime_file.Get(histoString.replace("$SAMPLE","data_obs"))
                histo_SM_Higgs = self.theABCD_VR_prime_file.Get(histoString.replace("$SAMPLE","SM_Higgs"))
                histo_TTbar = self.theABCD_VR_prime_file.Get(histoString.replace("$SAMPLE","TTbar"))
                histo_Others = self.theABCD_VR_prime_file.Get(histoString.replace("$SAMPLE","Others"))
                histo_nonQCD = histo_SM_Higgs.Clone()
                histo_nonQCD.Add(histo_TTbar)
                histo_nonQCD.Add(histo_Others)
                self.theABCD_VR_prime=ABCD()
                self.theABCD_VR_prime.setDataHisto(histo_data_obs)
                self.theABCD_VR_prime.setNonQCDHisto(histo_nonQCD)
                self.theABCD_VR_prime.subtractNonQCDFromData()
            elif dataMode=="simQCD":
                histo_QCD=self.theABCD_VR_prime_file.Get(histoString.replace("$SAMPLE","QCD"))
                self.theABCD_VR_prime=ABCD()
                self.theABCD_VR_prime.setQCDHisto(histo_QCD)

    def setWindow(self,region="A", X1=0.0, X2=0.0, Y1=0.0, Y2=0.0, verbose=0):
        if self.theABCD_SR!=None:
            self.theABCD_SR.setWindow(region, X1, X2, Y1, Y2, verbose)
        if self.theABCD_VR!=None:
            self.theABCD_VR.setWindow(region, X1, X2, Y1, Y2, verbose)
        if self.theABCD_SR_prime!=None:
            self.theABCD_SR_prime.setWindow(region, X1, X2, Y1, Y2, verbose)
        if self.theABCD_VR_prime!=None:
            self.theABCD_VR_prime.setWindow(region, X1, X2, Y1, Y2, verbose)

    def getFullClosureAndUncertainty(self, dataMode="data", verbose=0, applyAdditionalUncertainty=True):
        if self.theABCD_SR==None or self.theABCD_VR==None or self.theABCD_SR_prime==None or self.theABCD_VR_prime==None:
            print(self.theABCD_SR, self.theABCD_VR, self.theABCD_SR_prime, self.theABCD_VR_prime)
            return -1.0, 0.0
        if dataMode=="data":
            # get correction factor in VR_prime
            closure_VR_noCF, closure_VR_error_noCF = self.theABCD_VR.getClosureFactorWrtDataWithCorrectionAndAddUncertainty(1.0, 0.0)
            if verbose:
                print("VR WITHOUT CF closure_VR, closure_VR_error:", closure_VR_noCF, closure_VR_error_noCF)
            correctionFactor_VR, correctionFactor_VR_unc = self.theABCD_VR_prime.getCorrectionFactorWrtData()
            if verbose:
                print("VR_prime correctionFactor_VR, correctionFactor_VR_unc:", correctionFactor_VR, correctionFactor_VR_unc)
            # get closure in VR using the correction factor
            pred_A_VR, pred_A_VR_error = self.theABCD_VR.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc)
            if verbose:
                print("VR pred_A_VR, pred_A_VR_error:", pred_A_VR, pred_A_VR_error)
            closure_VR, closure_VR_error = self.theABCD_VR.getClosureFactorWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc)
            if verbose:
                print("VR WITH CF closure_VR, closure_VR_error:", closure_VR, closure_VR_error)
            
            add_unc=0.0
            if applyAdditionalUncertainty:
                add_unc=abs(1.0-closure_VR)
                if verbose:
                    print("additional uncertainty from VR: "+str(add_unc))
            # get correction factor in SR_prime
            correctionFactor_SR, correctionFactor_SR_unc = self.theABCD_SR_prime.getCorrectionFactorWrtData()
            if verbose:
                print("SR_prime correctionFactor_SR, correctionFactor_SR_unc:", correctionFactor_SR, correctionFactor_SR_unc)
            nonClosure_SR, nonClosure_SR_uncertainty = self.theABCD_SR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor_SR, correctionFactor_SR_unc, add_unc, verbose)
            if verbose:
                print("nonClosure_SR, nonClosure_SR_uncertainty:", nonClosure_SR, nonClosure_SR_uncertainty)
            return nonClosure_SR, nonClosure_SR_uncertainty

        if dataMode=="simQCD":
            # get correction factor in VR_prime
            correctionFactor_VR, correctionFactor_VR_unc = self.theABCD_VR_prime.getCorrectionFactorWrtSimQCD()
            if verbose:
                print("VR_prime correctionFactor_VR, correctionFactor_VR_unc:", correctionFactor_VR, correctionFactor_VR_unc)
            # get closure in VR using the correction factor
            pred_A_VR, pred_A_VR_error = self.theABCD_VR.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc)
            if verbose:
                print("VR pred_A_VR, pred_A_VR_error:", pred_A_VR, pred_A_VR_error)
            closure_VR, closure_VR_error = self.theABCD_VR.getClosureFactorWrtSimQCDWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc)
            if verbose:
                print("VR closure_VR, closure_VR_error:", closure_VR, closure_VR_error)
            add_unc=0.0
            if applyAdditionalUncertainty:
                add_unc=abs(1.0-closure_VR)
            # get correction factor in SR_prime
            correctionFactor_SR, correctionFactor_SR_unc = self.theABCD_SR_prime.getCorrectionFactorWrtSimQCD()
            if verbose:
                print("SR_prime correctionFactor_SR, correctionFactor_SR_unc:", correctionFactor_SR, correctionFactor_SR_unc)
            nonClosure_SR, nonClosure_SR_uncertainty = self.theABCD_SR.getNonClosureWrtSimQCDWithCorrectionAndAddUncertainty(correctionFactor_SR, correctionFactor_SR_unc, add_unc, verbose)
            if verbose:
                print("nonClosure_SR, nonClosure_SR_uncertainty:", nonClosure_SR, nonClosure_SR_uncertainty)
            return nonClosure_SR, nonClosure_SR_uncertainty

