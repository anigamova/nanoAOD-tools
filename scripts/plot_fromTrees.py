from ROOT import *; import glob, numpy as n; from array import array
import CMS_tdrStyle_lumi
CMS_tdrStyle_lumi.extraText       = "Work in progress"
isMC = False
#CMS_tdrStyle_lumi.setTDRStyle()
gStyle.SetOptStat(0);
gStyle.SetPalette(56)

MyFileNames = glob.glob("../output_gen_0p1/*.root")
ch_gen = TChain('Events');
for fName in  MyFileNames:
    ii = ch_gen.Add(fName);
titles = ['FWD',"0-75 GeV", "75-150 GeV", "150-250 GeV 0J", "150-250 GeV GE1J", ">250 GeV"]

#plot_var = 'HTXS_stage1_1_cat_pTjet30GeV';varName = 'STXS 1.2 bin'
#bins=[400+x for x in range(0,7)];binsArr = array("f", bins);

#plot_var = 'LHE_Vpt';varName = 'p_{T}(V) [GeV]'
#bins=[x*20 for x in range(0,20)]
#binsArr = array("f", bins+[500,600]);

plot_var = 'GenBJJ_mass';varName = 'm_{jj} [GeV]'
bins=[124 + x*0.1 for x in range(0,21)]
binsArr = array("f", bins);

HTXS_var_sm = TH1F('HTXS_var_sm',';'+varName+'; #sigma(a.u.)',len(binsArr)-1,binsArr)
HTXS_var_sm.Sumw2()
HTXS_var_sm.SetLineColor(kBlack)

HTXS_var_chq3 = TH1F('HTXS_var_chq3',';'+varName+'; #sigma(a.u.)',len(binsArr)-1,binsArr)
HTXS_var_chq3.Sumw2()
HTXS_var_chq3.SetLineColor(kBlue+2)
HTXS_var_chq3.SetMarkerColor(kBlue+2)
HTXS_var_chq3.SetLineWidth(2)

HTXS_var_chq1 = TH1F('HTXS_var_chq1',';'+varName+'; #sigma(a.u.)',len(binsArr)-1,binsArr)
HTXS_var_chq1.Sumw2()
HTXS_var_chq1.SetLineColor(kMagenta+3)
HTXS_var_chq1.SetMarkerColor(kMagenta+3)
HTXS_var_chq1.SetLineWidth(2)

HTXS_var_chu = TH1F('HTXS_var_chu',';'+varName+'; #sigma(a.u.)',len(binsArr)-1,binsArr)
HTXS_var_chu.Sumw2()
HTXS_var_chu.SetLineColor(kOrange+7)
HTXS_var_chu.SetMarkerColor(kOrange+3)
HTXS_var_chu.SetLineWidth(2)

HTXS_var_chd = TH1F('HTXS_var_chd',';'+varName+'; #sigma(a.u.)',len(binsArr)-1,binsArr)
HTXS_var_chd.Sumw2()
HTXS_var_chd.SetLineColor(kGreen+3)
HTXS_var_chd.SetMarkerColor(kGreen+3)
HTXS_var_chd.SetLineWidth(2)


for i in range(1, 7):
    HTXS_var_sm.GetXaxis().SetBinLabel(i,titles[i-1])
    HTXS_var_chq3.GetXaxis().SetBinLabel(i,titles[i-1])
    HTXS_var_chq1.GetXaxis().SetBinLabel(i,titles[i-1])
    HTXS_var_chu.GetXaxis().SetBinLabel(i,titles[i-1])
    HTXS_var_chd.GetXaxis().SetBinLabel(i,titles[i-1])

ch_gen.Draw(plot_var+">>HTXS_var_sm","(genWeight)")
ch_gen.Draw(plot_var+">>HTXS_var_chq3","(Reweights[4])")
ch_gen.Draw(plot_var+">>HTXS_var_chq1","(Reweights[2])")
ch_gen.Draw(plot_var+">>HTXS_var_chu","(Reweights[6])")
ch_gen.Draw(plot_var+">>HTXS_var_chd","(Reweights[8])")


_MYW = 800; _MYH = 600
_MYT = 0.05*_MYH; _MYB = 0.3*_MYH;
_MYL = 0.1*_MYW; _MYR = 0.05*_MYW


rp2 = TRatioPlot(HTXS_var_chq1,HTXS_var_sm,"divsym");
rp2.SetH2DrawOpt("hist");rp2.SetH1DrawOpt("hist");
rp2.Draw("nogrid");g2 = rp2.GetLowerRefGraph();

rp3 = TRatioPlot(HTXS_var_chu,HTXS_var_sm,"divsym");
rp3.SetH2DrawOpt("hist");rp3.SetH1DrawOpt("hist");
rp3.Draw("nogrid")
g3 = rp3.GetLowerRefGraph();

rp4 = TRatioPlot(HTXS_var_chd,HTXS_var_sm,"divsym");
rp4.SetH2DrawOpt("hist");rp4.SetH1DrawOpt("hist");
rp4.Draw("nogrid")
g4 = rp4.GetLowerRefGraph();


cB=0;
cB=TCanvas("cB","cB",_MYW,_MYH);
rp1 = TRatioPlot(HTXS_var_chq3,HTXS_var_sm,"divsym");
rp1.SetH2DrawOpt("hist");rp1.SetH1DrawOpt("hist");
rp1.GetLowYaxis().SetNdivisions(505)
rp1.SetLeftMargin( _MYL/_MYW );  rp1.SetRightMargin( _MYR/_MYW );
rp1.SetUpTopMargin( _MYT/_MYH );   rp1.SetLowBottomMargin( _MYB/_MYH );
rp1.Draw("nogrid");
rp1.GetUpperPad().cd()
HTXS_var_chq1.Draw("same")
HTXS_var_chu.Draw("same")
HTXS_var_chd.Draw("same")
rp1.GetLowerPad().cd()
g2.Draw("same P")
g3.Draw("same P")
g4.Draw("same P")
rp1.GetLowerRefYaxis().SetTitle("Ratio to SM");
rp1.GetUpperRefYaxis().SetTitle("#sigma(a.u.)");
cB.Update()
rp1.GetUpperPad().cd()
leg = TLegend(0.7, 0.55, 0.9, 0.8);
leg.SetTextFont(42)
leg.SetTextSize(0.04)
leg.SetBorderSize(0)
leg.AddEntry(HTXS_var_sm,"SM", "l");
leg.AddEntry(HTXS_var_chq3,"c_{hq3}= 0.1", "l");
leg.AddEntry(HTXS_var_chu,"c_{hu}= 0.1", "l");
leg.AddEntry(HTXS_var_chd,"c_{hd}= 0.1", "l");
leg.AddEntry(HTXS_var_chq1,"c_{hq1}= 0.1", "l");
leg.Draw('same');
cB.SaveAs(plot_var+'_EFT_4operators.pdf');

