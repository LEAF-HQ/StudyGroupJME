import os, math
from root_numpy import root2array, rec2array
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from array import array
from collections import OrderedDict

from tdrstyle_all import *
import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText = 'Private work'
TDR.cms_lumi_TeV = "Run 2, 138 fb^{-1}"


outdir=os.environ['ANALYSISPATH']+'/python/PDFs/'
os.system('mkdir -p '+outdir)

years = ['UL2018','UL2017','UL2016','EOY2018D','EOY2018ABC','EOY2017','EOY2016']
srs = ['ER1']
cats = ['tautau','mutau','etau']
event_vars = ['run','lumi','event', 'rho']
jet_vars = ['pt','eta','phi', 'pt_raw', 'jet_id', 'area']
lep_vars = ['pt','eta','phi', 'selector_bits']
n_jets = 7
n_leptons = 5
phi_max = 3.2
dR_min = 0.4
pt_min = 30
bins = {
    'HT':        (   30,     2500,     0,   range(30,2500,20)),
    'ST':        (   30,     2500,     0,   range(30,2500,20)),
    'ST_MET':    (  400,     2500,     0,   range(30,2500,20)),
    'tau_pt':    (   30,     1500,     0,   range(30,1500,20)),
    'leading_eta': ( -5.2,      5.2, 104,   []),
    'eta':       (   -5.2,      5.2,   0,   list(sorted([0.000, 0.522, 1.044, 1.566, 2.043, 2.500, 2.650, 2.853, 2.964, 3.139, 3.839, 5.191, -0.522, -1.044, -1.566, -2.043, -2.500, -2.650, -2.853, -2.964, -3.139, -3.839, -5.191]))),
    'pt':        (   30,     2500,     0,   range(30,2500,20)),
    'phi':       (-phi_max,  phi_max,  50,  []),
    'leading_phi':(-phi_max, phi_max,  50,  []),
    'd_eta':     (   -0.3,      0.3,   100, []),
    'd_phi':     (   -0.2,      0.2,   100, []),
    'd_pt':      ( -200,      200,     100, []),
    'd_R':       (    0.0,      0.3,   100, []),
    'R_eta':     (    0.5,      1.5,   200, []),
    'R_phi':     (    0.5,      1.5,   200, []),
    'R_pt':      (    0.5,      2.0,   250, []),
    'met_R_pt':  ( -100,      100,     50,  []),
    'met_R_phi': ( -100,      100,     50,  []),
    'met_d_pt':  (-1000,     1000,    100,  []),
    'met_d_phi': (    0.0,      3.5,   50,  []),
    'efrac':     (    0.0,      1.0,   20,  []),
}

colors = {
  'tautau': rt.kOrange+1,
  'mutau':  rt.kRed+1,
  'etau':   rt.kGreen+2,
  'eoy':    rt.kOrange+1,
  'ul':     rt.kGreen+2,
  'SG':     rt.kRed+1,
}
styles = {
    'eoy':     rt.kSolid,
    'ul':      rt.kSolid,
    'missing': rt.kDashed,
}



def deltaPhi(p1,p2):
    result = math.fabs(p1 - p2)
    if (result > math.pi): result = 2 * math.pi - result
    return result

def deltaEta(e1,e2):
  return math.fabs(e1 - e2)

def deltaR(e1,e2,p1,p2):
    if e1==0.0 and p1==0.0: return 10.0
    if e2==0.0 and p2==0.0: return 10.0
    de = deltaEta(e1=e1,e2=e2)
    dp = deltaPhi(p1=p1, p2=p2)
    return math.sqrt(de * de + dp * dp)

def deltaPhi_NP(p1,p2):
    result = np.fabs(p1 - p2)
    result[result> np.pi] = 2 * np.pi - result[result> np.pi]
    return result

def deltaEta_NP(e1,e2):
  return np.fabs(e1 - e2)

def deltaR_NP(e1,e2,p1,p2):
    de = deltaEta_NP(e1=e1,e2=e2)
    dp = deltaPhi_NP(p1=p1, p2=p2)
    dR = np.sqrt(de * de + dp * dp)
    dR.where(~((e1==0.0)&(p1==0.0)|((e2==0.0)&(p2==0.0))),10.0, inplace=True)
    return dR

def HistRatio(h1,h2,hr):
    for bin in range(h1.GetNbinsX()):
        val_n = h1.GetBinContent(bin)
        err_n = h1.GetBinError(bin)
        val_d = h2.GetBinContent(bin)
        err_d = h2.GetBinError(bin)
        hr.SetBinContent(bin, val_n/val_d if val_d!=0.0 else 0.0)
        hr.SetBinError(bin, (np.sqrt(err_n/val_n)**2+(err_d/val_d)**2) if val_d!=0.0 and val_n!=0.0 else 0)

def ResetCanvas(x_min,x_max,y_min,y_max, var, ylabel='Events', isRatio=False):
    if isRatio:
        canv = tdrDiCanvas(var, x_min,x_max,y_min,y_max, 0,2, var, ylabel, 'EOY/UL')
    else:
        canv = tdrCanvas(var, x_min,x_max,y_min,y_max, var, ylabel)
    leg = tdrLeg(0.65,0.7,0.9,0.9, textSize=0.035)
    lines = {}
    if isRatio:
        lines['ref'] = rt.TLine(x_min,1,x_max,1)
        lines['ref'].SetLineStyle(rt.kDotted)
        canv.cd(2)
        lines['ref'].Draw('same')
        canv.cd(1)
    if 'eta' in var:
        for x in [-2.964, -2.650, +2.650, +2.964]:
            lines[x] = rt.TLine(x,y_min,x,y_max)
            # lines[x].SetLineStyle(rt.kDotted)
            lines[x].Draw('same')
    return (canv,leg,lines)

def FilterEvents(df1,df2):
    to_remove=[]
    for i1, ev1 in df1.iterrows():
        ev2 = df2[(df2['run']==ev1['run'])*(df2['lumi']==ev1['lumi'])*(df2['cat']==ev1['cat'])]
        if len(ev2)!=1: to_remove.append(i1)
    df1_missing = df1[df1.index.isin(to_remove)]
    df1 = df1[~df1.index.isin(to_remove)]
    return df1, df1_missing


def SkipJet(ev, jet, n_leptons=3):
    if jet['pt']<pt_min: return True
    if jet['ID']<4: return True
    for t in range(n_leptons):
        for flav in ['ele','muo','tau']:
            if not ev[flav+'_ID_'+str(t)]: continue
            if ev[flav+'_pt_'+str(t)] < pt_min: continue
            lep_eta = ev[flav+'_eta_'+str(t)]
            lep_phi = ev[flav+'_phi_'+str(t)]
            dR = deltaR(jet['eta'],lep_eta,jet['phi'],lep_phi)
            if dR<dR_min:
                return True
    return False

def Save2DHist(canv,hist2D,name,var_x,var_y,line=None):
    SetAlternative2DColor(hist2D)
    hist2D.GetXaxis().SetNdivisions(505)
    hist2D.GetYaxis().SetNdivisions(505)
    hist2D.Draw('same colz')
    if line: line.Draw("same")
    rt.gPad.RedrawAxis()
    canv.SaveAs(outdir+name+'.pdf')
    vars_log = ['pt','ST','HT','MET']
    # if any([x in var_x for x in vars_log]):
    #     canv.SetLogx(True)
    # if any([x in var_y for x in vars_log]):
    #     canv.SetLogy(True)
    #     canv.SaveAs(outdir+name+'_log.pdf')
    canv.Close()

def Save1DHists(canv,leg,hists,pdfName,var_x):
    for name,h in hists.items():
        if h.GetEntries() == 0: continue
        lcolor = colors['ul'] if 'ul' in name else (colors['eoy'] if not 'SG' in name else colors['SG'])
        lstyle = styles['missing'] if 'missing' in name else styles['eoy']
        canv.cd(2 if 'ratio' in name else 1)
        tdrDraw(h, 'hist E', lcolor=lcolor, lstyle=lstyle, marker=0, fstyle=0)
        if not 'ratio' in name:
            leg.AddEntry(h, name.replace('missing', ' unmatch')+' tot:'+str(int(h.GetEntries())), 'l')
    rt.gPad.RedrawAxis()
    canv.SaveAs(outdir+pdfName+'.pdf')
    # vars_log = ['pt','ST','HT','MET']
    # if any([x in var_x for x in vars_log]):
    #     canv.cd(1).SetLogx(True)
    #     canv.cd(2).SetLogx(True)
    #     canv.SaveAs(outdir+pdfName+'_log.pdf')
    canv.Close()

def Plot2D(df1,df2,var1, var2, extraName1='', extraName2='', extraName='', nbinx=None, x_min=None, x_max=None, nbiny=None, y_min=None, y_max=None, binsx_=[],binsy_=[], nobj=-1):
    if extraName1!='': extraName1='_'+extraName1
    if extraName2!='': extraName2='_'+extraName2
    if extraName!='': extraName='_'+extraName
    pdfName = '2D_'+var1+'_'+var2+extraName1+extraName2+extraName
    if any([ x==None for x in [nbinx,x_min,x_max,nbiny,y_min,y_max]]):
        x_min, x_max, nbinx, binsx_ = bins[var1.replace('jet_','').replace('raw','').replace('met_','').replace('leading_pt','pt').replace('leading_','').replace('_SG','').replace('MET','pt').replace('tau_phi','phi').replace('tau_eta','eta').replace('ST_pt','ST')]
        y_min, y_max, nbiny, binsy_ = bins[var2.replace('jet_','').replace('raw','').replace('met_','').replace('leading_pt','pt').replace('leading_','').replace('_SG','').replace('MET','pt').replace('tau_phi','phi').replace('tau_eta','eta').replace('ST_pt','ST')]
    if len(binsx_)==0:
        hist2D   = rt.TH2F(pdfName, pdfName, nbinx, x_min, x_max, nbiny, y_min, y_max)
    else:
        hist2D   = rt.TH2F(pdfName, pdfName, len(binsx_)-1, array('d', binsx_), len(binsy_)-1, array('d', binsy_))
    if nobj<0:
        for index, ev1 in df1.iterrows():
            ev2 = df2[(df2['run']==ev1['run'])*(df2['lumi']==ev1['lumi'])*(df2['cat']==ev1['cat'])]
            if len(ev2)!=1: raise Exception('Event not found: unexpected behavior')
            ev2 = ev2.iloc[0]
            if nobj == -1:
                _ = hist2D.Fill(ev1[var1],ev2[var2])
            else:
                for i in range(int(np.fabs(nobj))):
                    jet1 = dict([(v,ev1['jet_'+v+'_'+str(i)]) for v in ['pt','eta','phi','ID']])
                    if SkipJet(ev1, jet1): continue
                    jetvar1 = ev1[var1+'_'+str(i)]
                    jetvar2 = None
                    for j in range(int(np.fabs(nobj))):
                        jet2 = dict([(v,ev2['jet_'+v+'_'+str(j)]) for v in ['pt','eta','phi', 'ID']])
                        if SkipJet(ev2, jet2): continue
                        dR = deltaR(jet1['eta'],jet2['eta'],jet1['phi'],jet2['phi'])
                        if dR<dR_min:
                            jetvar2 = ev2[var2+'_'+str(j)]
                            break
                    if jetvar2!= None:
                        _ = hist2D.Fill(jetvar1,jetvar2)
    else:
        for index, ev1 in df1.iterrows():
            for i in range(nobj):
                jet1 = dict([(v,ev1['jet_'+v+'_'+str(i)]) for v in ['pt','eta','phi','ID']])
                if SkipJet(ev1, jet1): continue
                jetvar1 = ev1[var1+'_'+str(i)]
                jetvar2 = ev1[var2+'_'+str(i)]
                _ = hist2D.Fill(jetvar1,jetvar2)
    canv = tdrCanvas(pdfName, x_min,x_max,y_min,y_max, var1+extraName1, var2+extraName2, square=kSquare, is2D=True, isExtraSpace=True)
    line = rt.TLine(x_min,x_min,y_max,y_max) if var1.replace('_SG','') == var2.replace('_SG','') else None
    Save2DHist(canv,hist2D,pdfName, var1, var2, line=line)


def Plot1D(dfs, var, varName=None, extraName='', nbinx=None, x_min=None, x_max=None, refName='ul', nobj=-1):
    if varName==None: varName=var
    if extraName!='': extraName='_'+extraName
    pdfName = '1D_'+var+extraName
    if any([ x==None for x in [nbinx,x_min,x_max]]):
        x_min, x_max, nbinx, binsx_ = bins[varName.replace('jet_','').replace('_raw','').replace('met_','').replace('ST_pt','ST') if not '_efrac' in varName else 'efrac']
    hist1D = OrderedDict()
    for name in dfs.keys():
        if len(binsx_)==0:
            for _name in ['','missing','ratio','SG', 'SGratio']:
                if 'SG' in _name and not (name=='eoy' and var=='ST_MET'): continue
                hist1D[name+_name] = rt.TH1F(pdfName+name+_name, pdfName+name+_name, nbinx, x_min, x_max)
        else:
            for _name in ['','missing','ratio','SG', 'SGratio']:
                if 'SG' in _name and not (name=='eoy' and var=='ST_MET'): continue
                hist1D[name+_name] = rt.TH1F(pdfName+name+_name, pdfName+name+_name, len(binsx_)-1, array('d', binsx_))
    for name1, df1 in dfs.items():
        for name2, df2 in dfs.items():
            if name1 == name2: continue
            for i1, ev1 in df1.iterrows():
                ev2 = df2[(df2['run']==ev1['run'])*(df2['lumi']==ev1['lumi'])*(df2['cat']==ev1['cat'])]
                if len(ev2)!=1: raise Exception('Event not found: unexpected behavior')
                ev2 = ev2.iloc[0]
                if nobj == -1:
                    _ = hist1D[name1].Fill(ev1[var])
                    if name1=='eoy' and var=='ST_MET':
                        _ = hist1D[name1+'SG'].Fill(ev1[var+'_SG'])
                else:
                    for i in range(nobj):
                        jet1 = dict([(v,ev1['jet_'+v+'_'+str(i)]) for v in ['pt','eta','phi','ID']])
                        if SkipJet(ev1, jet1): continue
                        var1 = ev1[var+'_'+str(i)]
                        var2 = None
                        for j in range(nobj):
                            jet2 = dict([(v,ev2['jet_'+v+'_'+str(j)]) for v in ['pt','eta','phi', 'ID']])
                            if SkipJet(ev2, jet2): continue
                            dR = deltaR(jet1['eta'],jet2['eta'],jet1['phi'],jet2['phi'])
                            if dR<dR_min:
                                var2 = ev2[var+'_'+str(j)]
                                break
                        if var2!= None:
                            R = (var1/var2) if (var2!=0) else (1 if var1==0 else 0)
                            d = (var1-var2) if (var!='phi')  else deltaPhi(var1,var2)
                            _ = hist1D[name1].Fill(var1)
                        else: _ = hist1D[name1+'missing'].Fill(var1)
    y_min,y_max = (8e-01,9e-01)
    for hist in hist1D.values():
        y_max = max(y_max,hist.GetMaximum())
    for name in dfs.keys():
        if name==refName: continue
        HistRatio(hist1D[name],hist1D[refName],hist1D[name+'ratio'])
        if name+'SG' in hist1D:
            HistRatio(hist1D[name+'SG'],hist1D[refName],hist1D[name+'SGratio'])
    canv,leg,lines = ResetCanvas(x_min,x_max,y_min,y_max*10, var, isRatio=True)
    canv.cd(1).SetLogy(True)
    Save1DHists(canv,leg,hist1D,pdfName,var)

def GetLepPtMin(year,cat, flav):
    if 'tau' == flav: pt = 40 if cat=='tautau' else 20
    if 'muo' == flav: pt = 23 if '16' in year else (25 if '17' in year else 28)
    if 'ele' == flav:
        if cat=='tautau':
            pt = 28 if '16' in year else (36 if '17' in year else 36)
        else: pt = 15
    return pt

def get_ref_file_path(cat, sr):
    # ref_file_path=os.environ['ANALYSISPATH']+'/python/PFDs/'
    ref_file_path=os.environ['ANALYSISPATH']+'/python/'
    fname = ref_file_path+'data__'+cat.replace('etau','eletau')+'_0b_'+sr+'_extended.csv'
    return fname

def LoadDF():
    branches = []
    columns  = []
    for flav in ['run','lumiblock', 'number', 'rho']:
        branches.append(flav)
        columns.append(flav.replace('block','').replace('number','event'))
    for flav in ['met','rawmet']:
        for var in ['pt','phi']:
            branches.append(flav+'.m_'+var)
            columns.append(flav+'_'+var)
    for flav in ['ak4chs']:
        for i_jet in range(n_jets):
            for var in jet_vars:
                branches.append(('jets_'+flav+'.m_'+var.replace('pt_raw','raw_factor')+'['+str(i_jet)+']', 0., 1))
                columns.append(flav.replace('ak4','').replace('chs','')+'jet_'+var.replace('jet_id','ID').replace('pt_raw','raw_factor')+'_'+str(i_jet))
    for flav in ['taus','electrons','muons']:
      for i_lep in range(n_leptons):
        for var in lep_vars:
            branches.append((flav+'.m_'+var+'['+str(i_lep)+']', 0., 1))
            columns.append(flav.replace('taus','tau').replace('electrons','ele').replace('muons','muo')+'_'+var+'_'+str(i_lep))
    df_all = []
    for sr in srs:
        for cat in cats:
            for year in years:
                era = 'eoy' if 'EOY' in year else 'ul'
                mymatrix = rec2array(root2array(filenames='../Inputs/NTuples_'+cat+'_0b_'+sr+'_'+year+'.root', treename='AnalysisTree', branches=branches))
                df = pd.DataFrame(mymatrix, columns=columns)
                df['year'] = year
                df['sr'] = sr
                df['cat'] = cat
                df['ST_MET_SG'] = 0.
                df['MET_SG'] = 0.
                df['leading_jet_pt_SG'] = 0.
                df['leading_tau_pt_SG'] = 0.
                df['subleading_tau_pt_SG'] = 0.
                ref_df = pd.read_csv(get_ref_file_path(cat, sr))
                for i1, ev_ref in ref_df.iterrows():
                    run, lumi = (ev_ref['run'],ev_ref[' lumi'])
                    ev = df[(df['run']==run)*(df['lumi']==lumi)]
                    if len(ev)!=1: continue
                    df.at[ev.index[0],'ST_MET_SG'] = ev_ref[' stmet']
                    df.at[ev.index[0],'MET_SG'] = ev_ref[' met']
                    df.at[ev.index[0],'leading_jet_pt_SG'] = ev_ref[' jet1_pt']
                    df.at[ev.index[0],'leading_tau_pt_SG'] = ev_ref[' e-mu-tau1_pt']
                    df.at[ev.index[0],'subleading_tau_pt_SG'] = ev_ref[' tau1-tau1-tau2_pt']
                for i_jet in range(n_jets):
                    df['rawjet_pt_'+str(i_jet)] = df['jet_pt_'+str(i_jet)] * df['jet_raw_factor_'+str(i_jet)]
                    for v in ['pt','eta','phi']:
                        df['original_jet_'+v+'_'+str(i_jet)] = df['jet_'+v+'_'+str(i_jet)]
                        df['jet_'+v+'_'+str(i_jet)].mask(df['jet_pt_'+str(i_jet)]<pt_min, 0.0, inplace=True)
                        df['jet_'+v+'_'+str(i_jet)].mask(df['jet_ID_'+str(i_jet)]<4, 0.0, inplace=True)
                for i_lep in range(n_leptons):
                    df['tau_ID_'+str(i_lep)] = df['tau_selector_bits_'+str(i_lep)].map(lambda x: bool(int(x) & 1<<4) )
                    df['ele_ID_'+str(i_lep)] = df['ele_selector_bits_'+str(i_lep)].map(lambda x: bool(int(x) & 1<<9) )
                    df['muo_ID_'+str(i_lep)] = df['muo_selector_bits_'+str(i_lep)].map(lambda x: bool(int(x) & 1<<11 & 1<<26) )
                    for flav in ['ele','muo','tau']:
                        for v in ['pt','eta','phi']:
                            df[flav+'_'+v+'_'+str(i_lep)].mask(df[flav+'_pt_'+str(i_lep)]<GetLepPtMin(year,cat,flav), 0.0, inplace=True)
                            df[flav+'_'+v+'_'+str(i_lep)].mask(~df[flav+'_ID_'+str(i_lep)], 0.0, inplace=True)
                        for i_jet in range(n_jets):
                            df['dR_jet_'+str(i_jet)+'_'+flav+'_'+str(i_lep)] = deltaR_NP(df['jet_eta_'+str(i_jet)], df[flav+'_eta_'+str(i_lep)], df['jet_phi_'+str(i_jet)], df[flav+'_phi_'+str(i_lep)])
                for i_jet in range(n_jets):
                    for v in ['pt','eta','phi']:
                        for i_lep in range(n_leptons):
                            for flav in ['ele','muo','tau']:
                                df['jet_'+v+'_'+str(i_jet)].mask(df['dR_jet_'+str(i_jet)+'_'+flav+'_'+str(i_lep)]<dR_min, 0.0, inplace=True)
                for flav in ['ele','muo','tau']:
                    df['leading_'+flav+'_pt'] = df[[flav+'_pt_'+str(i_lep) for i_lep in range(n_leptons)]].max(axis=1)
                    df['leading_'+flav]       = df[[flav+'_pt_'+str(i_lep) for i_lep in range(n_leptons)]].idxmax(axis=1)
                    for v in lep_vars:
                        df['leading_'+flav+'_'+v] = df.apply(lambda x: x[str(x['leading_'+flav]).replace('pt',v)], axis=1)
                for flav in ['jet']:
                    df['leading_'+flav+'_pt_ref'] = df[[flav+'_pt_'+str(i_jet) for i_jet in range(n_jets)]].max(axis=1)
                    df['leading_'+flav] = df[[flav+'_pt_'+str(i_jet) for i_jet in range(n_jets)]].idxmax(axis=1)
                    for v in jet_vars:
                        v = v.replace('jet_id','ID').replace('pt_raw','raw_factor')
                        df['leading_'+flav+'_'+v] = df.apply(lambda x: x[str(x['leading_'+flav]).replace('pt',v)], axis=1)
                for i_lep in range(n_leptons):
                    df['dtau_pt_'+str(i_lep)] = df.apply(lambda x: x['tau_pt_'+str(i_lep)] - x['leading_tau_pt'], axis=1)
                    df['dtau_pt_'+str(i_lep)].mask(df['dtau_pt_'+str(i_lep)]==0., -10000, inplace=True)
                df['subleading_tau_pt'] = df[['dtau_pt_'+str(i_lep) for i_lep in range(n_leptons)]].max(axis=1)
                df['subleading_tau_pt'] += df['leading_tau_pt']
                for i_lep in range(n_leptons):
                    df.drop('dtau_pt_'+str(i_lep), inplace=True, axis=1)
                df['HT'] = 0.
                for i_jet in range(n_jets):
                    df['HT'] = df['HT'] + df['jet_pt_'+str(i_jet)]
                df['ST'] = df['leading_jet_pt']+df['leading_tau_pt']
                df['ST'] = df['ST'] + df['subleading_tau_pt'].mask(df['cat']!='tautau', 0.)
                df['ST'] = df['ST'] + df['leading_ele_pt'].mask(df['cat']!='etau', 0.)
                df['ST'] = df['ST'] + df['leading_muo_pt'].mask(df['cat']!='mutau', 0.)
                df['ST_MET'] = df['ST'] + df['met_pt']
                df_all.append(df)
    return df_all

def PrintEvent(df,run,lumi,vars=['lead','ST','met'], isRef=False):
    if isRef:
        print df[(df['run']==run)*(df[' lumi'])==lumi][list(filter(lambda x: any([v in x for v in vars]), df.keys()))]
    else:
        print df[(df['run']==run)*(df['lumi'])==lumi][list(filter(lambda x: any([v in x for v in vars]), df.keys()))]

def main():
    df_all = pd.concat(LoadDF(), ignore_index=True)
    eoy = df_all[df_all['year'].str.contains('EOY')]
    ul = df_all[df_all['year'].str.contains('UL')]
    eoy, eoy_missing = FilterEvents(eoy,ul)
    ul, ul_missing = FilterEvents(ul,eoy)

    print df_all[(df_all['ST_MET']>1800)*(df_all['met_pt']>800)][['run', 'lumi', 'event','year','cat', 'ST_MET','met_pt']]

    # eoy-vs-SG-eoy
    Plot2D(df1=eoy, df2=eoy, var1='ST_MET_SG',         var2='ST_MET',          extraName1='eoy', extraName2='eoy', nbinx=100, x_min=400, x_max=6000, nbiny=100, y_min=400, y_max=6000)
    Plot2D(df1=eoy, df2=ul,  var1='ST_MET_SG',         var2='ST_MET',          extraName1='eoy', extraName2='ul',  nbinx=100, x_min=400, x_max=6000, nbiny=100, y_min=400, y_max=6000)
    Plot2D(df1=eoy, df2=ul,  var1='MET_SG',            var2='met_pt',          extraName1='eoy', extraName2='ul')
    Plot2D(df1=eoy, df2=ul,  var1='leading_jet_pt_SG', var2='leading_jet_pt',  extraName1='eoy', extraName2='ul')
    Plot2D(df1=eoy, df2=ul,  var1='leading_tau_pt_SG', var2='leading_tau_pt',  extraName1='eoy', extraName2='ul')

    # eoy-vs-ul
    for var in ['ST_MET', 'leading_jet_pt', 'leading_tau_pt', 'met_phi', 'rawmet_phi']:
        Plot2D(df1=eoy, df2=ul, var1=var, var2=var, extraName1='eoy', extraName2='ul')
        if 'leading' in var:
            # split in cat
            for cat in cats:
                Plot2D(df1=eoy[eoy['cat']==cat], df2=ul[ul['cat']==cat],  var1=var, var2=var,  extraName1='eoy', extraName2='ul', extraName=cat)

    for var in ['met_pt','leading_jet_pt','leading_tau_pt']:
        #STMET-vs-obj-pt in eoy and ul
        Plot2D(df1=eoy, df2=eoy, var1='ST_MET', var2=var, extraName1='eoy', extraName2='eoy', nbinx=100, x_min=400, x_max=6000, nbiny=100, y_min=10, y_max=2500)
        Plot2D(df1=eoy, df2=ul,  var1='ST_MET', var2=var, extraName1='eoy', extraName2='ul',  nbinx=100, x_min=400, x_max=6000, nbiny=100, y_min=10, y_max=2500)
        if var=='met_pt': continue
        #MET-vs-obj-pt in eoy and ul
        Plot2D(df1=eoy, df2=eoy, var1='met_pt', var2=var, extraName1='eoy', extraName2='eoy')
        Plot2D(df1=eoy, df2=ul,  var1='met_pt', var2=var, extraName1='eoy', extraName2='ul')
        if var=='leading_jet_pt': continue
        #lead-jet-lead-vs-lead-tau-pt in eoy and ul
        Plot2D(df1=eoy, df2=eoy, var1='leading_jet_pt', var2=var, extraName1='eoy', extraName2='eoy')
        Plot2D(df1=eoy, df2=ul,  var1='leading_jet_pt', var2=var, extraName1='eoy', extraName2='ul')

    #met-phi-vs-lead-(jet/tau)-phi in eoy and ul
    for flav in ['jet', 'tau']:
        Plot2D(df1=eoy, df2=eoy, var1='met_phi', var2='leading_'+flav+'_phi', extraName1='eoy', extraName2='eoy')
        Plot2D(df1=ul,  df2=ul,  var1='met_phi', var2='leading_'+flav+'_phi', extraName1='ul',  extraName2='ul')

    #var eoy-vs-ul
    for var in ['met_pt', 'rawmet_pt', 'HT', 'ST']:
        Plot2D(df1=eoy, df2=ul, var1=var, var2=var, extraName1='eoy', extraName2='ul')
        Plot2D(df1=eoy, df2=ul, var1=var, var2=var, extraName1='eoy', extraName2='ul', extraName='zoomed', nbinx=100, x_min=10, x_max=1500, nbiny=100, y_min=10, y_max=1500)
        Plot2D(df1=eoy, df2=ul, var1=var, var2=var, extraName1='eoy', extraName2='ul', extraName='extend', nbinx=100, x_min=10, x_max=6000, nbiny=100, y_min=10, y_max=6000)

    #lead jet eta-vs-(phi/pt) in eoy
    Plot2D(df1=eoy, df2=eoy,  var1='leading_jet_eta', var2='leading_jet_pt',  extraName1='eoy', extraName2='eoy')
    Plot2D(df1=eoy, df2=eoy,  var1='leading_jet_eta', var2='leading_jet_phi', extraName1='eoy', extraName2='eoy', nbinx=104, x_min=-5.2, x_max=5.2, nbiny=50, y_min=-phi_max, y_max=phi_max)

    #jet eta eoy-vs-ul
    Plot2D(df1=eoy, df2=ul,   var1='jet_eta',         var2='jet_eta',         extraName1='eoy', extraName2='ul', nobj=-n_jets)

    Plot2D(df1=eoy, df2=ul,   var1='jet_pt',          var2='jet_pt',          extraName1='eoy', extraName2='ul', nobj=-n_jets)
    #jet eta-vs-pt in eoy
    Plot2D(df1=eoy, df2=None, var1='jet_eta',         var2='jet_pt',          extraName1='eoy', extraName2='eoy', nobj=n_jets)

    dfs = OrderedDict([('eoy',eoy),('ul',ul)])
    for var in ['met_pt', 'rawmet_pt', 'HT', 'ST','ST_MET','leading_tau_pt']:
        Plot1D(dfs, var,var.replace('leading_','').replace('raw',''))
    for var in jet_vars:
        if 'id' in var: continue
        if 'raw' in var: continue
        if 'n_' in var: continue
        if 'area' in var: continue
        Plot1D(dfs, 'leading_jet_'+var, 'jet_'+var)
        Plot1D(dfs, 'jet_'+var,nobj=n_jets)
        if 'efrac' in var: continue
        if 'pt' in var: continue #case covered above
        Plot2D(df1=eoy, df2=ul, var1='leading_jet_'+var, var2='leading_jet_'+var, extraName1='eoy', extraName2='ul')

if __name__ == '__main__':
    main()

#
#
#
# PrintEvent(eoy,321434,454)
# PrintEvent(ref,321434,454, vars=[''], isRef=True)
