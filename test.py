import ROOT
import math
import ctypes
import numpy
ROOT.gROOT.SetBatch(True)
from class_ABCD import ABCD

histoFileName="all_histos_main_v6_invertedChi2.root"
histoFile=ROOT.TFile(histoFileName,"r")

theABCD_SR=ABCD()


histoString_SR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1"
histoString_VR="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_VR"
histoString_SR_invertedChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2"
histoString_VR_invertedChi2="boosted_Higgs_combi_chi2_A_phiSorted_01_$SAMPLE_sigprob_M15_1_invertedChi2_VR"

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


