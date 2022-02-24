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


##############################
doTests=False
if doTests:
    # test
    print("test")
    print("cut = 0.85")
    histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
    filename="all_histos_main_v6_invertedChi2.root"
    outputDir="NewScheme_scanOneDirection_M15_1/"
    dataMode="data"
    histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1"
    histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR"
    histoString_SR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2"
    histoString_VR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR"

    theFullABCD = FullABCD()
    theFullABCD.setup_ABCD("SR",filename,histoString_SR, dataMode )
    theFullABCD.setup_ABCD("VR",filename,histoString_VR, dataMode )
    theFullABCD.setup_ABCD("SR_prime",filename,histoString_SR_invChi2, dataMode )
    theFullABCD.setup_ABCD("VR_prime",filename,histoString_VR_invChi2, dataMode )

    cutValue=0.85
    windowMaxX=1.0
    windowMinX=0.3
    windowMaxY=1.0
    windowMinY=0.3
    otherDirectionCut=cutValue
    theFullABCD.setWindow("A", cutValue, windowMaxX, otherDirectionCut, windowMaxY)
    theFullABCD.setWindow("B", cutValue, windowMaxX, windowMinY, otherDirectionCut-0.0001)
    theFullABCD.setWindow("C", windowMinX, cutValue-0.0001, otherDirectionCut, windowMaxY)
    theFullABCD.setWindow("D", windowMinX, cutValue-0.0001, windowMinY, otherDirectionCut-0.0001)


    # get nonClosure in VR without the correction Factor
    closure_VR, closure_VR_error = theFullABCD.theABCD_VR.getClosureFactorWrtDataWithCorrectionAndAddUncertainty(1.0, 0.0, 0.0)
    print("VR No Correcton ->  closure_VR, closure_VR_error:", round(closure_VR,9), round(closure_VR_error,9))
    nonClosure_VR, nonClosure_error_VR = theFullABCD.theABCD_VR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(1.0, 0.0, 0.0)
    print("VR No Correcton ->  nonClosure_VR, nonClosure_error_VR:", round(nonClosure_VR,9), round(nonClosure_error_VR,9))

    # get correction factor in VR_prime
    correctionFactor_VR, correctionFactor_VR_unc = theFullABCD.theABCD_VR_prime.getCorrectionFactorWrtData()
    print("VR_prime correctionFactor_VR, correctionFactor_VR_unc:", round(correctionFactor_VR,9), round(correctionFactor_VR_unc,9))
    nonClosure_VR, nonClosure_error_VR = theFullABCD.theABCD_VR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc, 0.0)
    print("VR_prime nonClosure_VR, nonClosure_error_VR:", round(nonClosure_VR,9), round(nonClosure_error_VR,9))

    # get closure in VR using the correction factor
    # pred_A_VR, pred_A_VR_error = theFullABCD.theABCD_VR.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc)
    # print("VR pred_A_VR, pred_A_VR_error:", pred_A_VR, pred_A_VR_error)
    closure_VR, closure_VR_error = theFullABCD.theABCD_VR.getClosureFactorWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc)
    print("VR WITH Correcton -> closure_VR, closure_VR_error:", round(closure_VR,9), round(closure_VR_error,9))
    # get nonClosure in VR without the correction Factor
    nonClosure_VR, nonClosure_error_VR = theFullABCD.theABCD_VR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc, 0.0)
    print("VR WITH Correcton ->  nonClosure_VR, nonClosure_error_VR:", round(nonClosure_VR,9), round(nonClosure_error_VR,9))

    # add unc.
    add_unc = abs(1. - closure_VR)
    print("additional uncertainty: ", round(add_unc, 3))
    # get correction factor in SR
    correctionFactor_SR, correctionFactor_SR_unc = theFullABCD.theABCD_SR_prime.getCorrectionFactorWrtData()
    print("SR_prime correctionFactor_SR, correctionFactor_SR_unc:", round(correctionFactor_SR,9), round(correctionFactor_SR_unc,9))
    # total_uncertainty=math.sqrt(correctionFactor_SR_unc**2 + add_unc**2)
    # print("")
    # print("correction factor: "+ str(round(correctionFactor_SR)))
    # print("total uncertainty: "+ str(round(total_uncertainty,3)))
    # predict 
    pred_A, pred_A_error = theFullABCD.theABCD_SR.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor_SR,correctionFactor_SR_unc, add_unc, verbose=True)
    relative_pred_A_error=pred_A_error/float(pred_A)
    print("corrected pred_A, pred_A_error, relative_pred_A_error: "+str(round(pred_A,3))+" +- "+str(round(pred_A_error,3))+" ("+str(round(relative_pred_A_error,9))+")")

    histoFile.Close()


    print("")
    print("test")
    print("cut = 0.89")
    histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
    filename="all_histos_main_v6_invertedChi2.root"
    outputDir="NewScheme_scanOneDirection_M15_1/"
    dataMode="data"
    histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1"
    histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR"
    histoString_SR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2"
    histoString_VR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR"

    theFullABCD = FullABCD()
    theFullABCD.setup_ABCD("SR",filename,histoString_SR, dataMode )
    theFullABCD.setup_ABCD("VR",filename,histoString_VR, dataMode )
    theFullABCD.setup_ABCD("SR_prime",filename,histoString_SR_invChi2, dataMode )
    theFullABCD.setup_ABCD("VR_prime",filename,histoString_VR_invChi2, dataMode )

    cutValue=0.89
    windowMaxX=1.0
    windowMinX=0.3
    windowMaxY=1.0
    windowMinY=0.3
    otherDirectionCut=cutValue
    theFullABCD.setWindow("A", cutValue, windowMaxX, otherDirectionCut, windowMaxY)
    theFullABCD.setWindow("B", cutValue, windowMaxX, windowMinY, otherDirectionCut-0.0001)
    theFullABCD.setWindow("C", windowMinX, cutValue-0.0001, otherDirectionCut, windowMaxY)
    theFullABCD.setWindow("D", windowMinX, cutValue-0.0001, windowMinY, otherDirectionCut-0.0001)


    # get nonClosure in VR without the correction Factor
    closure_VR, closure_VR_error = theFullABCD.theABCD_VR.getClosureFactorWrtDataWithCorrectionAndAddUncertainty(1.0, 0.0, 0.0)
    print("VR No Correcton ->  closure_VR, closure_VR_error:", round(closure_VR,9), round(closure_VR_error,9))
    nonClosure_VR, nonClosure_error_VR = theFullABCD.theABCD_VR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(1.0, 0.0, 0.0)
    print("VR No Correcton ->  nonClosure_VR, nonClosure_error_VR:", round(nonClosure_VR,9), round(nonClosure_error_VR,9))

    # get correction factor in VR_prime
    correctionFactor_VR, correctionFactor_VR_unc = theFullABCD.theABCD_VR_prime.getCorrectionFactorWrtData()
    print("VR_prime correctionFactor_VR, correctionFactor_VR_unc:", round(correctionFactor_VR,9), round(correctionFactor_VR_unc,9))
    nonClosure_VR, nonClosure_error_VR = theFullABCD.theABCD_VR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc, 0.0)
    print("VR_prime nonClosure_VR, nonClosure_error_VR:", round(nonClosure_VR,4), round(nonClosure_error_VR,4))

    # get closure in VR using the correction factor
    # pred_A_VR, pred_A_VR_error = theFullABCD.theABCD_VR.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc)
    # print("VR pred_A_VR, pred_A_VR_error:", pred_A_VR, pred_A_VR_error)
    closure_VR, closure_VR_error = theFullABCD.theABCD_VR.getClosureFactorWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc)
    print("VR WITH Correcton -> closure_VR, closure_VR_error:", round(closure_VR,9), round(closure_VR_error,9))
    # get nonClosure in VR without the correction Factor
    nonClosure_VR, nonClosure_error_VR = theFullABCD.theABCD_VR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc, 0.0)
    print("VR WITH Correcton ->  nonClosure_VR, nonClosure_error_VR:", round(nonClosure_VR,9), round(nonClosure_error_VR,9))

    # add unc.
    add_unc = abs(1. - closure_VR)
    print("additional uncertainty: ", round(add_unc, 3))
    # get correction factor in SR
    correctionFactor_SR, correctionFactor_SR_unc = theFullABCD.theABCD_SR_prime.getCorrectionFactorWrtData(verbose=True)
    print("SR_prime correctionFactor_SR, correctionFactor_SR_unc:", round(correctionFactor_SR,8), round(correctionFactor_SR_unc,8))
    # total_uncertainty=math.sqrt(correctionFactor_SR_unc**2 + add_unc**2)
    # print("")
    # print("correction factor: "+ str(round(correctionFactor_SR,3)))
    # print("total uncertainty: "+ str(round(total_uncertainty,6)))
    # predict 
    pred_A, pred_A_error = theFullABCD.theABCD_SR.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor_SR,correctionFactor_SR_unc, add_unc, verbose=True)
    relative_pred_A_error=pred_A_error/float(pred_A)
    print("corrected pred_A, pred_A_error, relative_pred_A_error: "+str(round(pred_A,3))+" +- "+str(round(pred_A_error,3))+" ("+str(round(relative_pred_A_error,9))+")")

    histoFile.Close()



    print("")
    print("test restricted")
    print("cut = 0.6")
    histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
    filename="all_histos_main_v6_invertedChi2.root"
    outputDir="NewScheme_scanOneDirection_M15_1/"
    dataMode="data"
    histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1"
    histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR"
    histoString_SR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2"
    histoString_VR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR"

    theFullABCD = FullABCD()
    theFullABCD.setup_ABCD("SR",filename,histoString_SR, dataMode )
    theFullABCD.setup_ABCD("VR",filename,histoString_VR, dataMode )
    theFullABCD.setup_ABCD("SR_prime",filename,histoString_SR_invChi2, dataMode )
    theFullABCD.setup_ABCD("VR_prime",filename,histoString_VR_invChi2, dataMode )

    cutValue=0.6
    windowMaxX=0.7
    windowMinX=0.3
    windowMaxY=1.0
    windowMinY=0.3
    otherDirectionCut=0.7
    theFullABCD.setWindow("A", cutValue, windowMaxX, otherDirectionCut, windowMaxY)
    theFullABCD.setWindow("B", cutValue, windowMaxX, windowMinY, otherDirectionCut-0.0001)
    theFullABCD.setWindow("C", windowMinX, cutValue-0.0001, otherDirectionCut, windowMaxY)
    theFullABCD.setWindow("D", windowMinX, cutValue-0.0001, windowMinY, otherDirectionCut-0.0001)

    theFullABCD.theABCD_SR.printWindow("A")
    theFullABCD.theABCD_SR.printWindow("B")
    theFullABCD.theABCD_SR.printWindow("C")
    theFullABCD.theABCD_SR.printWindow("D")

    # get nonClosure in VR without the correction Factor
    closure_VR, closure_VR_error = theFullABCD.theABCD_VR.getClosureFactorWrtDataWithCorrectionAndAddUncertainty(1.0, 0.0, 0.0)
    print("VR No Correcton ->  closure_VR, closure_VR_error:", round(closure_VR,3), round(closure_VR_error,3))
    nonClosure_VR, nonClosure_error_VR = theFullABCD.theABCD_VR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(1.0, 0.0, 0.0)
    print("VR No Correcton ->  nonClosure_VR, nonClosure_error_VR:", round(nonClosure_VR,3), round(nonClosure_error_VR,3))

    # get correction factor in VR_prime
    correctionFactor_VR, correctionFactor_VR_unc = theFullABCD.theABCD_VR_prime.getCorrectionFactorWrtData()
    print("VR_prime correctionFactor_VR, correctionFactor_VR_unc:", round(correctionFactor_VR,3), round(correctionFactor_VR_unc,3))
    nonClosure_VR, nonClosure_error_VR = theFullABCD.theABCD_VR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc, 0.0)
    print("VR_prime nonClosure_VR, nonClosure_error_VR:", round(nonClosure_VR,3), round(nonClosure_error_VR,3))

    # get closure in VR using the correction factor
    # pred_A_VR, pred_A_VR_error = theFullABCD.theABCD_VR.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc)
    # print("VR pred_A_VR, pred_A_VR_error:", pred_A_VR, pred_A_VR_error)
    closure_VR, closure_VR_error = theFullABCD.theABCD_VR.getClosureFactorWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc)
    print("VR WITH Correcton -> closure_VR, closure_VR_error:", round(closure_VR,3), round(closure_VR_error,3))
    # get nonClosure in VR without the correction Factor
    nonClosure_VR, nonClosure_error_VR = theFullABCD.theABCD_VR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc, 0.0)
    print("VR WITH Correcton ->  nonClosure_VR, nonClosure_error_VR:", round(nonClosure_VR,3), round(nonClosure_error_VR,3))

    # add unc.
    add_unc = abs(1. - closure_VR)
    print("additional uncertainty: ", round(add_unc, 3))
    # get correction factor in SR
    correctionFactor_SR, correctionFactor_SR_unc = theFullABCD.theABCD_SR_prime.getCorrectionFactorWrtData(verbose=True)
    print("SR_prime correctionFactor_SR, correctionFactor_SR_unc:", round(correctionFactor_SR,8), round(correctionFactor_SR_unc,8))
    # total_uncertainty=math.sqrt(correctionFactor_SR_unc**2 + add_unc**2)
    # print("")
    # print("correction factor: "+ str(round(correctionFactor_SR,3)))
    # print("total uncertainty: "+ str(round(total_uncertainty,6)))
    # predict 
    pred_A, pred_A_error = theFullABCD.theABCD_SR.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor_SR,correctionFactor_SR_unc, add_unc, verbose=True)
    relative_pred_A_error=pred_A_error/float(pred_A)
    print("corrected pred_A, pred_A_error, relative_pred_A_error: "+str(round(pred_A,3))+" +- "+str(round(pred_A_error,3))+" ("+str(round(relative_pred_A_error,9))+")")

    histoFile.Close()



    # print("test unrestricted")
    # print("cut = 0.6")
    # histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
    # filename="all_histos_main_v6_invertedChi2.root"
    # outputDir="NewScheme_scanOneDirection_M15_1/"
    # dataMode="data"
    # histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1"
    # histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR"
    # histoString_SR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2"
    # histoString_VR_invChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR"

    # theFullABCD = FullABCD()
    # theFullABCD.setup_ABCD("SR",filename,histoString_SR, dataMode )
    # theFullABCD.setup_ABCD("VR",filename,histoString_VR, dataMode )
    # theFullABCD.setup_ABCD("SR_prime",filename,histoString_SR_invChi2, dataMode )
    # theFullABCD.setup_ABCD("VR_prime",filename,histoString_VR_invChi2, dataMode )

    # cutValue=0.6
    # windowMaxX=0.7
    # windowMinX=0.0
    # windowMaxY=1.0
    # windowMinY=0.0
    # otherDirectionCut=0.7
    # theFullABCD.setWindow("A", cutValue, windowMaxX, otherDirectionCut, windowMaxY)
    # theFullABCD.setWindow("B", cutValue, windowMaxX, windowMinY, otherDirectionCut-0.0001)
    # theFullABCD.setWindow("C", windowMinX, cutValue-0.0001, otherDirectionCut, windowMaxY)
    # theFullABCD.setWindow("D", windowMinX, cutValue-0.0001, windowMinY, otherDirectionCut-0.0001)

    # theFullABCD.theABCD_SR.printWindow("A")
    # theFullABCD.theABCD_SR.printWindow("B")
    # theFullABCD.theABCD_SR.printWindow("C")
    # theFullABCD.theABCD_SR.printWindow("D")

    # # get nonClosure in VR without the correction Factor
    # closure_VR, closure_VR_error = theFullABCD.theABCD_VR.getClosureFactorWrtDataWithCorrectionAndAddUncertainty(1.0, 0.0, 0.0)
    # print("VR No Correcton ->  closure_VR, closure_VR_error:", round(closure_VR,3), round(closure_VR_error,3))
    # nonClosure_VR, nonClosure_error_VR = theFullABCD.theABCD_VR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(1.0, 0.0, 0.0)
    # print("VR No Correcton ->  nonClosure_VR, nonClosure_error_VR:", round(nonClosure_VR,3), round(nonClosure_error_VR,3))

    # # get correction factor in VR_prime
    # correctionFactor_VR, correctionFactor_VR_unc = theFullABCD.theABCD_VR_prime.getCorrectionFactorWrtData()
    # print("VR_prime correctionFactor_VR, correctionFactor_VR_unc:", round(correctionFactor_VR,3), round(correctionFactor_VR_unc,3))
    # nonClosure_VR, nonClosure_error_VR = theFullABCD.theABCD_VR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc, 0.0)
    # print("VR_prime nonClosure_VR, nonClosure_error_VR:", round(nonClosure_VR,3), round(nonClosure_error_VR,3))

    # # get closure in VR using the correction factor
    # # pred_A_VR, pred_A_VR_error = theFullABCD.theABCD_VR.predictRegionAWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc)
    # # print("VR pred_A_VR, pred_A_VR_error:", pred_A_VR, pred_A_VR_error)
    # closure_VR, closure_VR_error = theFullABCD.theABCD_VR.getClosureFactorWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc)
    # print("VR WITH Correcton -> closure_VR, closure_VR_error:", round(closure_VR,3), round(closure_VR_error,3))
    # # get nonClosure in VR without the correction Factor
    # nonClosure_VR, nonClosure_error_VR = theFullABCD.theABCD_VR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor_VR, correctionFactor_VR_unc, 0.0)
    # print("VR WITH Correcton ->  nonClosure_VR, nonClosure_error_VR:", round(nonClosure_VR,3), round(nonClosure_error_VR,3))

    # # add unc.
    # add_unc = abs(1. - closure_VR)
    # print("additional uncertainty: ", round(add_unc, 3))
    # # get correction factor in SR
    # correctionFactor_SR, correctionFactor_SR_unc = theFullABCD.theABCD_SR_prime.getCorrectionFactorWrtData()
    # print("SR_prime correctionFactor_SR, correctionFactor_SR_unc:", round(correctionFactor_SR,3), round(correctionFactor_SR_unc,3))
    # total_uncertainty=math.sqrt(correctionFactor_SR_unc**2 + add_unc**2)
    # print("")
    # print("correction factor: "+ str(round(correctionFactor_SR)))
    # print("total uncertainty: "+ str(round(total_uncertainty,3)))
    # nonClosure_SR, nonClosure_SR_uncertainty = theFullABCD.theABCD_SR.getNonClosureWrtDataWithCorrectionAndAddUncertainty(correctionFactor_SR, correctionFactor_SR_unc, add_unc)
    # print("SR WITH COrrection: nonClosure_SR, nonClosure_SR_uncertainty ", round(nonClosure_SR,7), round(nonClosure_SR_uncertainty,7))
    # histoFile.Close()


# ##############################################################################################################################################################################

# create toplevel dir
toplevelDirName="outputs_sideband_scans"
if not os.path.isdir(toplevelDirName):
    os.mkdir(toplevelDirName)

masses=["15","40","55"]
lifetimes=["0","0p05","0p1","1",]

histoFile=ROOT.TFile("all_histos_main_v6_invertedChi2.root","r")
filename="all_histos_main_v6_invertedChi2.root"

for mass in masses:
    for lifetime in lifetimes:
        WP_tagName="M"+mass+"_"+lifetime
        print("")
        print("At "+WP_tagName)
        outputDir=toplevelDirName+"/"+"sideband_scan_"+WP_tagName+"/"
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

        # # limit window in x-direction and scan in x-direction
        # global_max=0.0
        # global_min=1.0

        # windowMinX=0.0 
        # windowMaxX=1.
        # windowMinY=0.0
        # windowMaxY=1.0
        # scanDirection="x"

        # scanMaxList=[ 0.69, 0.69, 0.84, 0.9,   ]
        # otherDirectionCutList=[ 0.7, 0.7, 0.85, 0.85  ]
        # colorList=[ ROOT.kRed, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1 ]
        # modes_list=[ "data", "simQCD", "simQCD", "simQCD",   ]
        # region_list=[ "SR", "SR", "SR", "SR",   ]
        # windowMaxXList=[ 0.7, 0.7, 0.85, 1.0,   ]
        # windowMinXList=[ 0.0, 0.0, 0.0, 0.0,   ]
        # graphsList=[]
        # assert(len(scanMaxList)==len(colorList))

        # graphXOffset=0.005*0

        # for i in range(len(scanMaxList)):
        #     print(i)

        #     cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection_NewScheme(scanDirection=scanDirection,
        #                             windowMinX=windowMinXList[i], windowMaxX=windowMaxXList[i], windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
        #                             scanMin=0.4, scanMax=scanMaxList[i],
        #                             filename=filename,
        #                             histoString_SR=histoString_SR, histoString_VR=histoString_VR, histoString_SR_invChi2=histoString_SR_invChi2, histoString_VR_invChi2=histoString_VR_invChi2,
        #                             dataMode=modes_list[i], verbosity=0 )
        #     graphsList.append(ROOT.TGraphErrors())
        #     # graphsList[-1].SetDirectory(0)
        #     # print(graphsList)
        #     for k in range(len(cutValue_list)):
        #         # print(cutValue_list[k], nonClosure_list[k])
        #         graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
        #         graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
        #     graphsList[-1].SetLineColor(colorList[i])
        #     graphsList[-1].SetName(region_list[i]+" "+modes_list[i]+" xRange=["+str(windowMinXList[i])+","+str(windowMaxXList[i])+"] yCut="+str(otherDirectionCutList[i])+" xCut=[0.5,"+str(scanMaxList[i])+"]")
        #     graphsList[-1].GetXaxis().SetTitle("ABCD x Cut")
        #     graphsList[-1].GetYaxis().SetTitle("closure")
        #     if modes_list[i]=="simQCD":
        #         graphsList[-1].SetLineStyle(2)
        #     else:
        #         graphsList[-1].SetLineStyle(1)    
        #     global_max=max(global_max,max(nonClosure_list))
        #     global_min=min(global_min, min(nonClosure_list))

        # canvas=ROOT.TCanvas("canvas_M15_1_xWindow","canvas_M15_1_xWindow",1024,768)
        # global_max=min(2.5, global_max)
        # global_min=0
        # # print(graphsList)
        # # print(graphsList[0])
        # graphsList[0].Print()
        # graphsList[0].SetMaximum(global_max*1.1)
        # graphsList[0].SetMinimum(global_min*0.9)
        # graphsList[0].GetXaxis().SetLimits(0.38, 0.92)
        # graphsList[0].Draw("ALP0")
        # for gr in graphsList[1:]:
        #     gr.Draw("same LP0")

        # canvas.BuildLegend()
        # canvas.SetGrid()
        # canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_"+WP_tagName+".png")
        # canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_"+WP_tagName+".pdf")

        # #########################################################################################################################################################################################################################

        # # Now limit the window to 0.7 (0.85) in y-direction and repeat the scan in x-diretion

        # global_max=0.0
        # global_min=1.0

        # windowMinX=0.0 
        # windowMaxX=1.
        # windowMinY=0.0
        # windowMaxY=1.0
        # scanDirection="x"

        # scanMaxList=[   0.95, 0.95, 0.95, 0.95 , 0.95,  ]
        # otherDirectionCutList=[  0.6, 0.6, 0.7, 0.85, 0.7   ]
        # colorList=[    ROOT.kRed, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1 , ROOT.kMagenta-5  ]
        # modes_list=[  "data", "simQCD", "simQCD","simQCD", "simQCD",    ]
        # region_list=[  "SR", "SR", "SR", "SR", "SR"   ]
        # windowMaxYList=[   0.7, 0.7, 0.85, 1.0, 1.0,   ]
        # windowMinYList=[ 0.0, 0.0, 0.0, 0.0, 0.0,   ]
        # graphsList=[]

        # graphXOffset=0.005*0

        # for i in range(len(scanMaxList)):
        #     cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection_NewScheme(scanDirection=scanDirection,
        #                             windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinYList[i], windowMaxY=windowMaxYList[i], otherDirectionCut=otherDirectionCutList[i],
        #                             scanMin=0.5, scanMax=scanMaxList[i],
        #                             filename=filename,
        #                             histoString_SR=histoString_SR, histoString_VR=histoString_VR, histoString_SR_invChi2=histoString_SR_invChi2, histoString_VR_invChi2=histoString_VR_invChi2,
        #                             dataMode=modes_list[i], verbosity=0 )
        #     graphsList.append(ROOT.TGraphErrors())
        #     for k in range(len(cutValue_list)):
        #         graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
        #         graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
        #     graphsList[-1].SetLineColor(colorList[i])
        #     graphsList[-1].SetName(region_list[i]+" "+modes_list[i]+" yRange=["+str(windowMinYList[i])+","+str(windowMaxYList[i])+"] yCut="+str(otherDirectionCutList[i])+" xCut=[0.5,"+str(scanMaxList[i])+"]")
        #     graphsList[-1].GetXaxis().SetTitle("ABCD x-Cut")
        #     graphsList[-1].GetYaxis().SetTitle("closure")
        #     if modes_list[i]=="simQCD":
        #         graphsList[-1].SetLineStyle(2)
        #     else:
        #         graphsList[-1].SetLineStyle(1)    
        #     global_max=max(global_max,max(nonClosure_list))
        #     global_min=min(global_min, min(nonClosure_list))

        # canvas=ROOT.TCanvas("canvas_M15_1_yWindow","canvas_M15_1_yWindow",1024,768)
        # global_max=min(2.5, global_max)
        # global_min=0
        # graphsList[0].SetMaximum(global_max*1.1)
        # graphsList[0].SetMinimum(global_min*0.9)
        # graphsList[0].Draw("ALP0")
        # for gr in graphsList[1:]:
        #     gr.Draw("same LP0")

        # canvas.BuildLegend()
        # canvas.SetGrid()
        # canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_"+WP_tagName+".png")
        # canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_"+WP_tagName+".pdf")

        #########################################################################################################################################################################################################################
        #########################################################################################################################################################################################################################
        ##########################################################################################################################################################################################################################
        #########################################################################################################################################################################################################################
        #########################################################################################################################################################################################################################
        ##########################################################################################################################################################################################################################

        # Now everything again but with lower bound of 0.3 in both directions

        # global_max=0.0
        # global_min=1.0

        # windowMinX=0.3 
        # windowMaxX=1.
        # windowMinY=0.3
        # windowMaxY=1.0
        # scanDirection="x"

        # scanMaxList=[ 0.69, 0.69, 0.84, 0.9,   ]
        # otherDirectionCutList=[ 0.7, 0.7, 0.85, 0.85  ]
        # colorList=[ ROOT.kRed, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1 ]
        # modes_list=[ "data", "simQCD", "simQCD", "simQCD",   ]
        # region_list=[ "SR", "SR", "SR", "SR",   ]
        # windowMaxXList=[ 0.7, 0.7, 0.85, 1.0,   ]
        # windowMinXList=[ 0.3   ]*4
        # graphsList=[]
        # assert(len(scanMaxList)==len(colorList))

        # graphXOffset=0.005*0

        # for i in range(len(scanMaxList)):
        #     print(i)
        #     verb=0
        #     if i==0:
        #         verb=1

        #     cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection_NewScheme(scanDirection=scanDirection,
        #                             windowMinX=windowMinXList[i], windowMaxX=windowMaxXList[i], windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
        #                             scanMin=0.5, scanMax=scanMaxList[i],
        #                             filename=filename,
        #                             histoString_SR=histoString_SR, histoString_VR=histoString_VR, histoString_SR_invChi2=histoString_SR_invChi2, histoString_VR_invChi2=histoString_VR_invChi2,
        #                             dataMode=modes_list[i], verbosity=verb )
        #     graphsList.append(ROOT.TGraphErrors())
        #     for k in range(len(cutValue_list)):
        #         graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
        #         graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
        #     graphsList[-1].SetLineColor(colorList[i])
        #     graphsList[-1].SetName(region_list[i]+" "+modes_list[i]+" xRange=["+str(windowMinXList[i])+","+str(windowMaxXList[i])+"] yCut="+str(otherDirectionCutList[i])+" xCut=[0.5,"+str(scanMaxList[i])+"]")
        #     graphsList[-1].GetXaxis().SetTitle("ABCD x Cut")
        #     graphsList[-1].GetYaxis().SetTitle("closure")
        #     if modes_list[i]=="simQCD":
        #         graphsList[-1].SetLineStyle(2)
        #     else:
        #         graphsList[-1].SetLineStyle(1)    
        #     global_max=max(global_max,max(nonClosure_list))
        #     global_min=min(global_min, min(nonClosure_list))

        # canvas=ROOT.TCanvas("canvas_M15_1_xWindow","canvas_M15_1_xWindow",1024,768)
        # global_max=min(2.5, global_max)
        # global_min=0
        # # print(graphsList)
        # # print(graphsList[0])
        # graphsList[0].Print()
        # graphsList[0].SetMaximum(global_max*1.1)
        # graphsList[0].SetMinimum(global_min*0.9)
        # graphsList[0].GetXaxis().SetLimits(0.38, 0.92)
        # graphsList[0].Draw("ALP0")
        # for gr in graphsList[1:]:
        #     gr.Draw("same LP0")

        # canvas.BuildLegend()
        # canvas.SetGrid()
        # canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_lowRestriction_0p3_"+WP_tagName+".png")
        # canvas.SaveAs(outputDir+"scan_xDirection_xRestricted_lowRestriction_0p3_"+WP_tagName+".pdf")


        # #########################################################################################################################################################################################################################
        # #########################################################################################################################################################################################################################
        # #########################################################################################################################################################################################################################

        # # Now limit the window to 0.7 (0.85) in y-direction and repeat the scan in x-diretion

        # global_max=0.0
        # global_min=1.0

        # windowMinX=0.3
        # windowMaxX=1.
        # windowMinY=0.3
        # windowMaxY=1.0
        # scanDirection="x"

        # scanMaxList=[   0.95, 0.95, 0.95, 0.95 , 0.95,  ]
        # otherDirectionCutList=[  0.6, 0.6, 0.7, 0.85, 0.7   ]
        # colorList=[    ROOT.kRed, ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1 , ROOT.kMagenta-5  ]
        # modes_list=[  "data", "simQCD", "simQCD","simQCD", "simQCD",    ]
        # region_list=[  "SR", "SR", "SR", "SR", "SR"   ]
        # windowMaxYList=[   0.7, 0.7, 0.85, 1.0, 1.0,   ]
        # windowMinYList=[ 0.3   ]*5
        # graphsList=[]

        # graphXOffset=0.005*0

        # for i in range(len(scanMaxList)):
        #     cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection_NewScheme(scanDirection=scanDirection,
        #                             windowMinX=windowMinX, windowMaxX=windowMaxX, windowMinY=windowMinYList[i], windowMaxY=windowMaxYList[i], otherDirectionCut=otherDirectionCutList[i],
        #                             scanMin=0.5, scanMax=scanMaxList[i],
        #                             filename=filename,
        #                             histoString_SR=histoString_SR, histoString_VR=histoString_VR, histoString_SR_invChi2=histoString_SR_invChi2, histoString_VR_invChi2=histoString_VR_invChi2,
        #                             dataMode=modes_list[i], verbosity=0 )
        #     graphsList.append(ROOT.TGraphErrors())
        #     for k in range(len(cutValue_list)):
        #         graphsList[-1].SetPoint(k,cutValue_list[k]+graphXOffset,nonClosure_list[k])
        #         graphsList[-1].SetPointError(k,0,non_Closure_error_list[k])
        #     graphsList[-1].SetLineColor(colorList[i])
        #     graphsList[-1].SetName(region_list[i]+" "+modes_list[i]+" yRange=["+str(windowMinYList[i])+","+str(windowMaxYList[i])+"] yCut="+str(otherDirectionCutList[i])+" xCut=[0.5,"+str(scanMaxList[i])+"]")
        #     graphsList[-1].GetXaxis().SetTitle("ABCD x-Cut")
        #     graphsList[-1].GetYaxis().SetTitle("closure")
        #     if modes_list[i]=="simQCD":
        #         graphsList[-1].SetLineStyle(2)
        #     else:
        #         graphsList[-1].SetLineStyle(1)    
        #     global_max=max(global_max,max(nonClosure_list))
        #     global_min=min(global_min, min(nonClosure_list))

        # canvas=ROOT.TCanvas("canvas_M15_1_yWindow","canvas_M15_1_yWindow",1024,768)
        # global_max=min(2.5, global_max)
        # global_min=0
        # graphsList[0].SetMaximum(global_max*1.1)
        # graphsList[0].SetMinimum(global_min*0.9)
        # graphsList[0].Draw("ALP0")
        # for gr in graphsList[1:]:
        #     gr.Draw("same LP0")

        # canvas.BuildLegend()
        # canvas.SetGrid()
        # canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_lowRestriction_0p3_"+WP_tagName+".png")
        # canvas.SaveAs(outputDir+"scan_xDirection_yRestricted_lowRestriction_0p3_"+WP_tagName+".pdf")




        #########################################################################################################################################################################################################################
        #########################################################################################################################################################################################################################
        ##########################################################################################################################################################################################################################
        #########################################################################################################################################################################################################################
        #########################################################################################################################################################################################################################
        ##########################################################################################################################################################################################################################

        # Now scan x-range for different y-cuts

        global_max=0.0
        global_min=1.0

        windowMinX=0.3 
        windowMaxX=1.
        windowMinY=0.3
        windowMaxY=1.0
        scanDirection="x"

        scanMaxList=[ 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69  ]
        otherDirectionCutList=[ 0.6, 0.65, 0.7, 0.75, 0.8, 0.83, 0.85, 0.87, 0.89, 0.91 ]
        colorList=[  ROOT.kOrange, ROOT.kMagenta+1,  ROOT.kRed, ROOT.kMagenta-5, ROOT.kCyan, ROOT.kOrange-1, ROOT.kGreen, ROOT.kGreen+3, ROOT.kCyan-8,ROOT.kBlue, ]
        modes_list=[ "data", "data", "data", "data", "data","data", "data", "data", "data", "data"  ]
        region_list=[ "SR", "SR", "SR", "SR", "SR","SR", "SR", "SR", "SR", "SR"  ]
        windowMaxXList=[ 0.7, 0.7, 0.7, 0.7, 0.7,0.7,0.7, 0.7,0.7,0.7,]
        windowMinXList=[ 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,0.3,0.3,    ]
        graphsList=[]
        assert(len(scanMaxList)==len(colorList))


        for i in range(len(scanMaxList)):
            print(i)
            verb=0
            if i==9:
                verb=1
            graphXOffset=0.0000*i

            cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection_NewScheme(scanDirection=scanDirection,
                                    windowMinX=windowMinXList[i], windowMaxX=windowMaxXList[i], windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                                    scanMin=0.4, scanMax=scanMaxList[i],
                                    filename=filename,
                                    histoString_SR=histoString_SR, histoString_VR=histoString_VR, histoString_SR_invChi2=histoString_SR_invChi2, histoString_VR_invChi2=histoString_VR_invChi2,
                                    dataMode=modes_list[i], verbosity=verb )
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

        canvas=ROOT.TCanvas("canvas_M15_1_xWindow","canvas_M15_1_xWindow",1024,768)
        global_max=min(2.5, global_max)
        global_min=0.8
        # print(graphsList)
        # print(graphsList[0])
        graphsList[0].Print()
        graphsList[0].SetMaximum(global_max*1.1)
        graphsList[0].SetMinimum(global_min*0.9)
        graphsList[0].GetXaxis().SetLimits(0.38, 0.72)
        graphsList[0].Draw("ALP0")
        for gr in graphsList[1:]:
            gr.Draw("same LP0")

        canvas.BuildLegend()
        canvas.SetGrid()
        canvas.SaveAs(outputDir+"scan_xDirection_yDirection_lowerRestriction_03_"+WP_tagName+".png")
        canvas.SaveAs(outputDir+"scan_xDirection_yDirection_lowerRestriction_03_"+WP_tagName+".pdf")


        # Now restrict y<0.7 and scan whole x-range

        global_max=0.0
        global_min=1.0

        windowMinX=0.3 
        windowMaxX=1.
        windowMinY=0.3
        windowMaxY=0.7
        scanDirection="x"

        scanMaxList=[ 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, ]
        otherDirectionCutList=[ 0.4,  0.5, 0.6, 0.65, 0.67, 0.69, ]
        colorList=[ ROOT.kCyan, ROOT.kBlue, ROOT.kMagenta+1, ROOT.kMagenta-5, ROOT.kRed, ROOT.kGreen, ]
        modes_list=[ "data",  "data", "data", "data", "data", "data",  ]
        region_list=[ "SR",  "SR", "SR", "SR", "SR", "SR",   ]
        windowMaxXList=[ 1.0, 1.0,1.0,1.0,1.0,1.0,]
        windowMinXList=[ 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,  ]
        graphsList=[]
        assert(len(scanMaxList)==len(colorList))

        graphXOffset=0.0005*1

        for i in range(len(scanMaxList)):
            print(i)
            verb=0
            # if i==9:
            #     verb=1

            cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection_NewScheme(scanDirection=scanDirection,
                                    windowMinX=windowMinXList[i], windowMaxX=windowMaxXList[i], windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                                    scanMin=0.4, scanMax=scanMaxList[i],
                                    filename=filename,
                                    histoString_SR=histoString_SR, histoString_VR=histoString_VR, histoString_SR_invChi2=histoString_SR_invChi2, histoString_VR_invChi2=histoString_VR_invChi2,
                                    dataMode=modes_list[i], verbosity=verb )
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

        canvas=ROOT.TCanvas("canvas_M15_1_xWindow","canvas_M15_1_xWindow",1024,768)
        global_max=min(2.5, global_max)
        global_min=0.8
        # print(graphsList)
        # print(graphsList[0])
        # graphsList[0].Print()
        graphsList[0].SetMaximum(global_max*1.1)
        graphsList[0].SetMinimum(global_min*0.9)
        graphsList[0].GetXaxis().SetLimits(0.38, 1.0)
        graphsList[0].Draw("ALP0")
        for gr in graphsList[1:]:
            gr.Draw("same LP0")

        canvas.BuildLegend()
        canvas.SetGrid()
        canvas.SaveAs(outputDir+"scan_xDirection_yDirection_yRestricted_lowerRestriction_03_"+WP_tagName+".png")
        canvas.SaveAs(outputDir+"scan_xDirection_yDirection_yRestricted_lowerRestriction_03_"+WP_tagName+".pdf")




        #########################################################################################################################################################################################################################
        #########################################################################################################################################################################################################################
        ##########################################################################################################################################################################################################################
        #########################################################################################################################################################################################################################
        #########################################################################################################################################################################################################################
        ##########################################################################################################################################################################################################################

        # Now without adding the additional uncertainty of the validaton regions
        # Now scan x-range for different y-cuts

        global_max=0.0
        global_min=1.0

        windowMinX=0.3 
        windowMaxX=1.
        windowMinY=0.3
        windowMaxY=1.0
        scanDirection="x"

        scanMaxList=[ 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69  ]
        otherDirectionCutList=[ 0.6, 0.65, 0.7, 0.75, 0.8, 0.83, 0.85, 0.87, 0.89, 0.91  ]
        colorList=[ ROOT.kOrange, ROOT.kMagenta+1,  ROOT.kRed, ROOT.kMagenta-5, ROOT.kCyan, ROOT.kOrange-1, ROOT.kGreen, ROOT.kGreen+3, ROOT.kCyan-8,  ROOT.kBlue, ]
        modes_list=[ "data", "data", "data", "data", "data","data", "data", "data", "data", "data"  ]
        region_list=[ "SR", "SR", "SR", "SR", "SR","SR", "SR", "SR", "SR", "SR"  ]
        windowMaxXList=[ 0.7, 0.7, 0.7, 0.7, 0.7,0.7,0.7, 0.7,0.7,0.7,]
        windowMinXList=[ 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,0.3,0.3,    ]
        graphsList=[]
        assert(len(scanMaxList)==len(colorList))

        graphXOffset=0.005*0

        for i in range(len(scanMaxList)):
            print(i)
            verb=0
            if i==9:
                verb=1

            cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection_NewScheme(scanDirection=scanDirection,
                                    windowMinX=windowMinXList[i], windowMaxX=windowMaxXList[i], windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                                    scanMin=0.4, scanMax=scanMaxList[i],
                                    filename=filename,
                                    histoString_SR=histoString_SR, histoString_VR=histoString_VR, histoString_SR_invChi2=histoString_SR_invChi2, histoString_VR_invChi2=histoString_VR_invChi2,
                                    dataMode=modes_list[i], verbosity=verb, applyAdditionalUncertainty=False )
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

        canvas=ROOT.TCanvas("canvas_M15_1_xWindow","canvas_M15_1_xWindow",1024,768)
        global_max=min(2.5, global_max)
        global_min=0.8
        # print(graphsList)
        # print(graphsList[0])
        graphsList[0].Print()
        graphsList[0].SetMaximum(global_max*1.1)
        graphsList[0].SetMinimum(global_min*0.9)
        graphsList[0].GetXaxis().SetLimits(0.38, 0.72)
        graphsList[0].Draw("ALP0")
        for gr in graphsList[1:]:
            gr.Draw("same LP0")

        canvas.BuildLegend()
        canvas.SetGrid()
        canvas.SaveAs(outputDir+"scan_xDirection_yDirection_lowerRestriction_03_noAddVRUncertainty_"+WP_tagName+".png")
        canvas.SaveAs(outputDir+"scan_xDirection_yDirection_lowerRestriction_03_noAddVRUncertainty_"+WP_tagName+".pdf")


        # Now restrict y<0.7 and scan whole x-range

        global_max=0.0
        global_min=1.0

        windowMinX=0.3 
        windowMaxX=1.
        windowMinY=0.3
        windowMaxY=0.7
        scanDirection="x"

        scanMaxList=[ 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, ]
        otherDirectionCutList=[ 0.4,  0.5, 0.6, 0.65, 0.67, 0.69, ]
        colorList=[ ROOT.kCyan, ROOT.kBlue, ROOT.kMagenta+1, ROOT.kMagenta-5, ROOT.kRed, ROOT.kGreen, ]
        modes_list=[ "data",  "data", "data", "data", "data", "data",  ]
        region_list=[ "SR",  "SR", "SR", "SR", "SR", "SR",   ]
        windowMaxXList=[ 1.0, 1.0,1.0,1.0,1.0,1.0,]
        windowMinXList=[ 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,  ]
        graphsList=[]
        assert(len(scanMaxList)==len(colorList))

        graphXOffset=0.005*0

        for i in range(len(scanMaxList)):
            print(i)
            verb=0
            # if i==0:
            #     verb=1

            cutValue_list, nonClosure_list, non_Closure_error_list = scanClosureOneDirection_NewScheme(scanDirection=scanDirection,
                                    windowMinX=windowMinXList[i], windowMaxX=windowMaxXList[i], windowMinY=windowMinY, windowMaxY=windowMaxY, otherDirectionCut=otherDirectionCutList[i],
                                    scanMin=0.4, scanMax=scanMaxList[i],
                                    filename=filename,
                                    histoString_SR=histoString_SR, histoString_VR=histoString_VR, histoString_SR_invChi2=histoString_SR_invChi2, histoString_VR_invChi2=histoString_VR_invChi2,
                                    dataMode=modes_list[i], verbosity=verb, applyAdditionalUncertainty=False )
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

        canvas=ROOT.TCanvas("canvas_M15_1_xWindow","canvas_M15_1_xWindow",1024,768)
        global_max=min(2.5, global_max)
        global_min=0.8
        # print(graphsList)
        # print(graphsList[0])
        # graphsList[0].Print()
        graphsList[0].SetMaximum(global_max*1.1)
        graphsList[0].SetMinimum(global_min*0.9)
        graphsList[0].GetXaxis().SetLimits(0.38, 1.0)
        graphsList[0].Draw("ALP0")
        for gr in graphsList[1:]:
            gr.Draw("same LP0")

        canvas.BuildLegend()
        canvas.SetGrid()
        canvas.SaveAs(outputDir+"scan_xDirection_yDirection_yRestricted_lowerRestriction_03_noAddVRUncertainty_"+WP_tagName+".png")
        canvas.SaveAs(outputDir+"scan_xDirection_yDirection_yRestricted_lowerRestriction_03_noAddVRUncertainty_"+WP_tagName+".pdf")

