import math
import time
import os
import sys
import ROOT
import random
import warnings

import rootlogon

from array import array

ROOT.gSystem.Load('exroot/lib/libExRootAnalysis.so')
ROOT.gSystem.Load('delphes/libDelphes.so')
ROOT.gSystem.Load('libpdf.so')

ROOT.gSystem.CompileMacro('higgs_drv.cc', 'O', 'higgs_drv')
ROOT.gSystem.Load('higgs_drv.so')

ROOT.gROOT.SetBatch(True)

import numpy as np

nevts     = 50000
verbose   = True
GeV       = 1000.0
nexternal = 6
ninitial  = 2
afloat_t  = np.float64
aint_t    = np.int32

def read_evgn(path, nevts, doboost=False):
    f = ROOT.TFile.Open(path, 'read')
    t = f.Get('physics')
    nevts = min(nevts, t.GetEntries())
    event_data = np.zeros((nevts*nexternal*4,), dtype=afloat_t)
    event_info = np.zeros((nevts*3,), dtype=afloat_t)
    evtlist = sorted(random.sample(range(t.GetEntries()), nevts))
    for ievt, evtnum in enumerate(evtlist):
        t.GetEntry(evtnum)
        lm, vbar, lp, v = None, None, None, None
        ga, gb = ROOT.TLorentzVector(), ROOT.TLorentzVector()
        ga.SetPtEtaPhiM(t.mc_pt[0]/GeV, t.mc_eta[0], t.mc_phi[0], t.mc_m[0]/GeV)
        gb.SetPtEtaPhiM(t.mc_pt[1]/GeV, t.mc_eta[1], t.mc_phi[1], t.mc_m[1]/GeV)
        for ip in xrange(t.mc_n):
            if t.mc_pdgId[ip] in [11, 13] and t.mc_status[ip] == 3:
                lm = ROOT.TLorentzVector()
                lm.SetPtEtaPhiM(t.mc_pt[ip]/GeV, t.mc_eta[ip], t.mc_phi[ip], t.mc_m[ip]/GeV)
            elif t.mc_pdgId[ip] in [-11, -13] and t.mc_status[ip] == 3:
                lp = ROOT.TLorentzVector()
                lp.SetPtEtaPhiM(t.mc_pt[ip]/GeV, t.mc_eta[ip], t.mc_phi[ip], t.mc_m[ip]/GeV)
            elif t.mc_pdgId[ip] in [12, 14] and t.mc_status[ip] == 3:
                v = ROOT.TLorentzVector()
                v.SetPtEtaPhiM(t.mc_pt[ip]/GeV, t.mc_eta[ip], t.mc_phi[ip], t.mc_m[ip]/GeV)
            elif t.mc_pdgId[ip] in [-12, -14] and t.mc_status[ip] == 3:
                vbar = ROOT.TLorentzVector()
                vbar.SetPtEtaPhiM(t.mc_pt[ip]/GeV, t.mc_eta[ip], t.mc_phi[ip], t.mc_m[ip]/GeV)
            if None not in [lm,lp,v,vbar]:
                break
        h = lm + vbar + lp + v
        if doboost:
            boost = h.BoostVector()
            for p in [lm, lp, v, vbar, h, ga, gb]:
                p.Boost( -boost.X(), -boost.Y(), 0)
            T = np.sqrt( (h.Pz()*h.Pz()) + (h.M()*h.M()) )
            EA = (T + h.Pz()) / 2
            EB = (T - h.Pz()) / 2
            ga.SetPxPyPzE(0,0, EA,EA)
            gb.SetPxPyPzE(0,0,-EB,EB)
        if ievt < 10:
            print 'INFO :: EVGN sanity check mH =', h.M()
            print 'INFO :: EVGN sanity check momentum conserved /', h.Pt() - (ga+gb).Pt()
        particles = [ga, gb, lp, v, lm, vbar]
        for ip,p in enumerate(particles):
            event_data[ievt*nexternal*4+ip*4+0] = p.E()
            event_data[ievt*nexternal*4+ip*4+1] = p.Px()
            event_data[ievt*nexternal*4+ip*4+2] = p.Py()
            event_data[ievt*nexternal*4+ip*4+3] = p.Pz()
        event_info[ievt*3+0] = h.M()
        event_info[ievt*3+1] = h.Rapidity()
        event_info[ievt*3+2] = h.E()

    f.Close()
    return event_data, nevts, event_info

def read_lhef(path, nevts):    
    f = ROOT.TFile.Open(path, 'read')
    t = f.Get('Delphes')
    nevts = min(nevts, t.GetEntries())
    event_data = np.zeros((nevts*nexternal*4,), dtype=afloat_t)
    event_info = np.zeros((nevts*3,), dtype=afloat_t)
    evtlist = sorted(random.sample(range(t.GetEntries()), nevts))
    for ievt, evtnum in enumerate(evtlist):
        t.GetEntry(evtnum)
        ga = t.Particle[0]
        gb = t.Particle[1]
        lm, vbar, lp, v = None, None, None, None
        n = t.Particle.GetEntries()
        for ip in xrange(0, n):
            if t.Particle[ip].PID   in [-11, -13]:
                lp = t.Particle[ip]
            elif t.Particle[ip].PID in [ 11,  13]:
                lm = t.Particle[ip]
            elif t.Particle[ip].PID in [-12, -14]:
                vbar = t.Particle[ip]
            elif t.Particle[ip].PID in [ 12,  14]:
                v = t.Particle[ip]
        h = ROOT.TLorentzVector()
        for p in [lp, v, lm, vbar]:
            h += ROOT.TLorentzVector(p.Px, p.Py, p.Pz, p.E)
        if ievt < 10:
            print 'INFO :: LHEF sanity check mH =', h.M()
            print 'INFO :: LHEF sanity check momentum conserved /', h.Pt() - np.sqrt(pow(ga.Px+gb.Px, 2) + pow(ga.Px+gb.Px, 2))
        ## g g > e+ ve mu- vm~ 
        particles = [ga, gb, lp, v, lm, vbar]
        
        # print h.M()
        for ip,p in enumerate(particles):
            event_data[ievt*nexternal*4+ip*4+0] = p.E
            event_data[ievt*nexternal*4+ip*4+1] = p.Px
            event_data[ievt*nexternal*4+ip*4+2] = p.Py
            event_data[ievt*nexternal*4+ip*4+3] = p.Pz
        event_info[ievt*3+0] = h.M()
        event_info[ievt*3+1] = h.Rapidity()
        event_info[ievt*3+2] = h.E()
    f.Close()
    return event_data, nevts, event_info

print
print

from ROOT import PDFGrid
from ROOT import quadme

pdf_obj = PDFGrid.load('CTEQ6L', PDFGrid.kgluon)
pdf_dat = np.ndarray(shape=(pdf_obj.get_nx()*pdf_obj.get_nQ(),), dtype=afloat_t, buffer=pdf_obj.get_grid_data(PDFGrid.kgluon))
        
ca, kSM, kHWW, kAWW, Lambda = 0.5, 1.0, 0.0, 2.0, 125.0
ca, kSM, kHWW, kAWW, Lambda = afloat_t(ca), afloat_t(kSM), afloat_t(kHWW), afloat_t(kAWW), afloat_t(Lambda)

dummy_float_arry = np.array([0],dtype=afloat_t)
dummy_int_arry = np.array([0],dtype=aint_t)

event_data, nevts, info_lhef = read_lhef('/global/dschouten/inputs/spin/ReweightingTest_MG5/emmupDF.root', nevts)

output_lhef_num = np.zeros((nevts,), dtype=afloat_t)
output_lhef_den = np.zeros((nevts,), dtype=afloat_t)

ca = afloat_t(0.5)

quadme(output_lhef_num, event_data, nevts, dummy_float_arry, dummy_int_arry, 0, 
       pdf_dat, pdf_obj.get_nx(), pdf_obj.get_nQ(), np.log(pdf_obj.get_xmin()), np.log(pdf_obj.get_xmax()), np.log(pdf_obj.get_Qmin()), np.log(pdf_obj.get_Qmax()),
       ca, kSM, kHWW, kAWW, Lambda, -1)

ca = afloat_t(0.75)

quadme(output_lhef_den, event_data, nevts, dummy_float_arry, dummy_int_arry, 0, 
       pdf_dat, pdf_obj.get_nx(), pdf_obj.get_nQ(), np.log(pdf_obj.get_xmin()), np.log(pdf_obj.get_xmax()), np.log(pdf_obj.get_Qmin()), np.log(pdf_obj.get_Qmax()),
       ca, kSM, kHWW, kAWW, Lambda, -1)

event_data, nevts, info_evgn = read_evgn('/global/dschouten/inputs/spin/ReweightingTest_MG5/leplepDF_200k.root', nevts, False)

output_evgn_num = np.zeros((nevts,), dtype=afloat_t)
output_evgn_den = np.zeros((nevts,), dtype=afloat_t)

ca = afloat_t(0.5)

quadme(output_evgn_num, event_data, nevts, dummy_float_arry, dummy_int_arry, 0, 
       pdf_dat, pdf_obj.get_nx(), pdf_obj.get_nQ(), np.log(pdf_obj.get_xmin()), np.log(pdf_obj.get_xmax()), np.log(pdf_obj.get_Qmin()), np.log(pdf_obj.get_Qmax()),
       ca, kSM, kHWW, kAWW, Lambda, -1)

ca = afloat_t(0.75)

quadme(output_evgn_den, event_data, nevts, dummy_float_arry, dummy_int_arry, 0, 
       pdf_dat, pdf_obj.get_nx(), pdf_obj.get_nQ(), np.log(pdf_obj.get_xmin()), np.log(pdf_obj.get_xmax()), np.log(pdf_obj.get_Qmin()), np.log(pdf_obj.get_Qmax()),
       ca, kSM, kHWW, kAWW, Lambda, -1)


ofile = ROOT.TFile.Open('/tmp/reweight_cmp.root', 'recreate')

evgn_ntup = ROOT.TNtuple('evgn', '', 'me:mh')
lhef_ntup = ROOT.TNtuple('lhef', '', 'me:mh')

res = [evgn_ntup.Fill(x, y) for x,y in zip(output_evgn_num/output_evgn_den, info_evgn[::3])]
res = [lhef_ntup.Fill(x, y) for x,y in zip(output_lhef_num/output_lhef_den, info_lhef[::3])]

evgn_ntup.Write()
lhef_ntup.Write()

ofile.Close()

ifile = ROOT.TFile.Open('/tmp/reweight_cmp.root', 'read')

ROOT.gROOT.cd()

can = ROOT.TCanvas()

ifile.Get('evgn').Draw('log(1+me)>>EVGN(100,0,1)', '', 'enorm')
ifile.Get('lhef').Draw('log(1+me)>>LHEF(100,0,1)', '', 'normsame')

ROOT.gROOT.FindObject('EVGN').SetLineColor(ROOT.kRed)
ROOT.gROOT.FindObject('EVGN').SetMarkerColor(ROOT.kRed)

leg = ROOT.TLegend(0.25, 0.7, 0.5, 0.9)
leg.AddEntry(ROOT.gROOT.FindObject('EVGN'), 'D3PD (after Pythia)', 'lp')
leg.AddEntry(ROOT.gROOT.FindObject('LHEF'), 'LHEF', 'l')
leg.SetFillColor(0)
leg.Draw('SAME')

can.Print('/tmp/reweight_cmp.pdf')
