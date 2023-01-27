#include <vector>
#include <TMath.h>

Double_t bwconv(Double_t *x, Double_t *par)
{
   
   //par[0]=Mass
   //par[1]=Gamma
   //par[2]=Area
   //par[3]=Width (sigma) of convoluted Gaussian function
   // ax*x+bx+c
   //par[4]=a
   //par[5]=b
   //par[6]=c

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fbw;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;
      Double_t background;
      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];
      Double_t xb;
      xb = x[0];
      step = (xupp-xlow) / np;

      // Convolution integral of BW and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fbw = .5*par[1]/(3.14159*((xx-par[0])*(xx-par[0]) + .25*par[1]*par[1]));
         sum += fbw * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fbw = .5*par[1]/(3.14159*((xx-par[0])*(xx-par[0]) + .25*par[1]*par[1]));
         sum += fbw * TMath::Gaus(x[0],xx,par[3]);
      }
      background = par[4]*xb*xb + par[5]*xb + par[6];
      return (par[2]*step*sum*invsq2pi/par[3]); }

Double_t Combine(Double_t *x, Double_t *par)
{
   
   //par[0]=Mass
   //par[1]=Gamma
   //par[2]=Area
   //par[3]=Width (sigma) of convoluted Gaussian function
   // ax*x+bx+c
   //par[4]=a
   //par[5]=b
   //par[6]=c

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fbw;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;
      Double_t background;
      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];
      Double_t xb;
      xb = x[0];
      step = (xupp-xlow) / np;

      // Convolution integral of BW and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fbw = .5*par[1]/(3.14159*((xx-par[0])*(xx-par[0]) + .25*par[1]*par[1]));
         sum += fbw * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fbw = .5*par[1]/(3.14159*((xx-par[0])*(xx-par[0]) + .25*par[1]*par[1]));
         sum += fbw * TMath::Gaus(x[0],xx,par[3]);
      }
      background = exp(par[4] + par[5]*x[0]);
      return (background + par[2]*step*sum*invsq2pi/par[3]); }

// define fit BW
Double_t fitBreitWigner(Double_t *x, Double_t *par) {
    Double_t Mean=par[1];
    Double_t Width=par[2];
    Double_t gamma=sqrt(Mean*Mean*(Mean*Mean+Width*Width));
    Double_t con=par[0];
    return con / (pow(x[0]*x[0]-Mean*Mean,2)+Mean*Mean*Width*Width); }

// define fit gauss
Double_t fitGaus(Double_t *x, Double_t *par) {
    Double_t arg=(x[0]-par[1])/par[2];
    return par[0]*exp(-0.5*arg*arg); }

// combine fit
Double_t fitBoth(double *x, double *par) {
    Double_t arg = (x[0]-par[1])/par[2];
    Double_t Gaus =  par[0]*exp(-0.5*arg*arg);
    Double_t Mean=par[4];
    Double_t Width=par[5];
    Double_t gamma=sqrt(Mean*Mean*(Mean*Mean+Width*Width));
    Double_t con=par[3];
    Double_t BW =  con / (pow(x[0]*x[0]-Mean*Mean,2)+Mean*Mean*Width*Width);
    return Gaus + BW; }

Double_t fitExp(double *x, double *par) {
    Double_t arg = par[0]+par[1]*x[0];
    return exp(arg); }

Double_t CombFit(double *x, double *par) {
    Double_t arg = (x[0]-par[1])/par[2];
    Double_t Gaus =  par[0]*exp(-0.5*arg*arg);
    Double_t Mean=par[4];
    Double_t Width=par[5];
    Double_t gamma=sqrt(Mean*Mean*(Mean*Mean+Width*Width));
    Double_t con=par[3];
    Double_t BW =  con / (pow(x[0]*x[0]-Mean*Mean,2)+Mean*Mean*Width*Width);
    Double_t arg1 = par[6]+par[7]*x[0];
    Double_t expo = exp(arg1);
    return Gaus + BW + expo; }
    

void PlotTest() {
    gROOT->LoadMacro("AtlasUtils.C");
    gROOT->SetStyle("ATLAS");
    gStyle->SetOptStat();
    gROOT->ForceStyle();
    gROOT->SetBatch(kTRUE);

    ROOT::EnableImplicitMT();

    // We prepare an input tree to run on
    auto fileName = "MyxAODAnalysis.outputs.root"; // Zmumu without selector
    auto fileName1 = "MyxAODAnalysis1.outputs.root"; // Ztautau
    auto fileName2 = "MyxAODAnalysis2.outputs.root"; // Diboson WZlvvv
    auto fileName3 = "MyxAODAnalysis3.outputs.root"; // Single top
    auto fileName4 = "MyxAODAnalysis4.outputs.root"; // ttbar
    auto fileName5 = "MyxAODAnalysis5.outputs.root"; // Drell-Yan 120-180
    auto fileName6 = "MyxAODAnalysis6.outputs.root"; // Drell-Yan 180-250
    auto fileName7 = "MyxAODAnalysis7.outputs.root"; // Drell-Yan 250-400
    auto fileName8 = "MyxAODAnalysis8.outputs.root"; // Drell-Yan 400-600
    auto fileName9 = "MyxAODAnalysis9.outputs.root"; // Drell-Yan 800-1000
    auto fileName10 = "MyxAODAnalysis10.outputs.root"; // Drell-Yan 1000-1250
    auto fileName11 = "MyxAODAnalysis11.outputs.root"; // Drell-Yan 600-800
    auto fileName12 = "MyxAODAnalysis12.outputs.root"; // Diboson WWlvlv
    auto fileName13 = "MyxAODAnalysis13.outputs.root"; // Diboson WZlvll
    auto fileName14 = "MyxAODAnalysis14.outputs.root"; // Diboson ZZllll
    auto fileName15 = "MyxAODAnalysis15.outputs.root"; // Diboson ZZvvvv
    auto fileName16 = "MyxAODAnalysis16.outputs.root"; // Dijet JF17
    auto fileName17 = "MyxAODAnalysis17.outputs.root"; // Dijet JF23
    auto fileName18 = "MyxAODAnalysis18.outputs.root"; // Dijet JF35
    auto fileName19 = "MyxAODAnalysis19.outputs.root"; // Dijet JF50 (not in use)
    auto fileName20 = "Data1.outputs.root"; // Data loose
    auto fileName21 = "ZtoMuons.outputs.root"; // Zmumu with selector medium
    auto fileName22 = "Data2.outputs.root"; // Data loose
    auto fileName23 = "Data3.outputs.root"; // Data loose
    auto fileName24 = "Higgs1.outputs.root"; // ggH125mumu loose
    auto fileName25 = "Data4.outputs.root"; // Data loose
    auto fileName26 = "Data5.outputs.root"; // Data loose
    auto fileName27 = "Data6.outputs.root"; // Data loose
    auto fileName28 = "Data7.outputs.root"; // Data loose
    auto fileName29 = "Data8.outputs.root"; // Data loose
    auto fileName30 = "Data9.outputs.root"; // Data loose
    auto fileName31 = "Data10.outputs.root"; // Data loose
    auto fileName32 = "Data11.outputs.root"; // Data loose
    auto fileName33 = "Data12.outputs.root"; // Data loose
    auto fileName34 = "Data13.outputs.root"; // Data loose
    auto fileName35 = "Data14.outputs.root"; // Data loose
    auto fileName36 = "Data15.outputs.root"; // Data loose

    auto treeName = "analysis";

    // We read the tree from the file and create a RDataFrame.
    ROOT::RDataFrame d(treeName, fileName, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass", "mChildren", "Z_match", "MinBias", "truthZ", "children", "truthEta", "truthPhi", "matchEta", "matchPhi",
    "ptRes", "etaRes" });
    // Plot MuonEta
    auto MuonEta = d.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    // Plot MuonPhi
    auto MuonPhi = d.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    // Plot MuonPt
    auto MuonPt = d.Histo1D({"MuonPt","",40,0.0,150},{"MuonPt"});
    // Plot MuonE
    auto MuonE = d.Histo1D({"MuonE","",40,0.0,400},{"MuonE"});
    // Plot leading muon Pt
    auto leadingPt = d.Histo1D({"Leading muon Pt","",40,0.0,150},{"leadingPt"});
    // Plot sub-leading muon Pt
    auto subPt = d.Histo1D({"sub-leading muon Pt","",40,0.0,150},{"subPt"});
    // plot the amount of muons in each event
    auto MuonSize = d.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    // scatter plot of Pt
    auto MuonScatter = d.Histo2D({"MuonScattering","",40,0.0,150,40,0.0,150},{"leadingPt"},{"subPt"});
    // plot the Z mass
    auto Z_mass = d.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});
    // plot true Z mass 
    auto truthZ = d.Histo1D({"truthZ","",150,0.0,150},{"truthZ"});
    // plot children Pt
    auto children = d.Histo1D({"children","",40,0.0,150},{"children"});
    // plot muon truth matches
    auto mChildren = d.Histo1D({"mChildren","",40,0.0,150},{"mChildren"});
    // plot matched Z mass
    auto Z_match = d.Histo1D({"Z_match","",150,0.0,150},{"Z_match"});
    // Plot true Eta
    auto truthEta = d.Histo1D({"truthEta","",40,-5,5},{"truthEta"});
    // Plot true Phi
    auto truthPhi = d.Histo1D({"truthPhi","",20,-3.14,3.14},{"truthPhi"});
    // Plot matched Eta
    auto matchEta = d.Histo1D({"matchEta","",40,-5,5},{"matchEta"});
    // Plot matched Phi
    auto matchPhi = d.Histo1D({"matchPhi","",20,-3.14,3.14},{"matchPhi"});
    // test cuts
    //auto pt_cut = d.Define("pt_cut","MuonPt[MuonPt > 20 && MuonEta > -2.7 && MuonEta < 2.7]").Histo1D({"Hist_MuonPt","",40,0.0,150},{"pt_cut"});
    // test cuts
    auto eta_cut = d.Define("eta_cut","MuonEta[MuonPt > 20]").Histo1D({"Hist_MuonEta","",40,-2.7,2.7},{"eta_cut"});
    // match phi cuts
    auto mphi_cut = d.Define("mphi_cut","matchPhi[matchEta > -2.7 && matchEta < 2.7]").Histo1D({"Hist_matchPhi","",40,-3.14,3.14},{"mphi_cut"});
    // truth phi cuts
    auto tphi_cut = d.Define("tphi_cut","truthPhi[truthEta > -2.7 && truthEta < 2.7]").Histo1D({"Hist_truthPhi","",40,-3.14,3.14},{"tphi_cut"});
    // match pt cuts
    auto mpt_cut = d.Define("mpt_cut","mChildren[matchEta > -2.7 && matchEta < 2.7 && matchPhi > -3.14 && matchPhi < 3.14]").Histo1D({"Hist_mChildren","",40,0.0,100},{"mpt_cut"});
    // truth pt cuts
    auto tpt_cut = d.Define("tpt_cut","children[truthEta > -2.7 && truthEta < 2.7 && truthPhi > -3.14 && truthPhi < 3.14]").Histo1D({"Hist_Children","",40,0.0,100},{"tpt_cut"});
    // plot pt resolution
    auto ptRes = d.Histo1D({"ptRes","",40,-0.02,0.02},{"ptRes"});
    // plot eta resolution
    auto etaRes = d.Histo1D({"etaRes","",40,-2.7,2.7},{"etaRes"});
    // scatter plot of truth Pt vs eta
    //auto scatterTruth = d.Histo2D({"scatterTruth","",40,0.0,150,40,-2.7,2.7},{"children"},{"truthEta"});
    // scatter plot of match Pt vs eta
    //auto scatterMatch = d.Histo2D({"scatterMatch","",40,0.0,150,40,-2.7,2.7},{"mChildren"},{"matchEta"});
    // min bias
    auto MinBias = d.Histo1D({"MinBias","",150,0.0,150},{"MinBias"});

    ROOT::RDataFrame d1(treeName, fileName1, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize1 = d1.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass1 = d1.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d2(treeName, fileName2, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize2 = d2.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass2 = d2.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d3(treeName, fileName3, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize3 = d3.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass3 = d3.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d4(treeName, fileName4, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize4 = d4.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass4 = d4.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d5(treeName, fileName5, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize", 
    "Z_mass"});
    auto MuonSize5 = d5.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass5 = d5.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d6(treeName, fileName6, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize", 
    "Z_mass"});
    auto MuonSize6 = d6.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass6 = d6.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d7(treeName, fileName7, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize", 
    "Z_mass"});
    auto MuonSize7 = d7.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass7 = d7.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d8(treeName, fileName8, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize", 
    "Z_mass"});
    auto MuonSize8 = d8.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass8 = d8.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d9(treeName, fileName9, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize9 = d9.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass9 = d9.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d10(treeName, fileName10, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize10 = d10.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass10 = d10.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d11(treeName, fileName11, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize11 = d11.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass11 = d11.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d12(treeName, fileName12, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize12 = d12.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass12 = d12.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d13(treeName, fileName13, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize13 = d13.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass13 = d13.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d14(treeName, fileName14, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize14 = d14.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass14 = d14.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d15(treeName, fileName15, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize15 = d15.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass15 = d15.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d16(treeName, fileName16, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize16 = d16.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass16 = d16.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d17(treeName, fileName17, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSiz",
    "Z_mass"});
    auto MuonSize17 = d17.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass17 = d17.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d18(treeName, fileName18, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize18 = d18.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass18 = d18.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d19(treeName, fileName19, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonSize19 = d19.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass19 = d19.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d20(treeName, fileName20, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta20 = d20.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi20 = d20.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonPt20 = d20.Histo1D({"MuonPt","",40,0.0,150},{"MuonPt"});
    auto leadingPt20 = d20.Histo1D({"Leading muon Pt","",40,0.0,150},{"leadingPt"});
    auto subPt20 = d20.Histo1D({"sub-leading muon Pt","",40,0.0,150},{"subPt"});
    auto MuonSize20 = d20.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass20 = d20.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});
    auto Z_mass20c = d20.Define("Z_massCut","Z_mass[leadingPt > 15 && subPt > 15]").Histo1D({"Z_massHist","",150,0.0,150},{"Z_massCut"});

    ROOT::RDataFrame d21(treeName, fileName21, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass", "mChildren", "Z_match", "MinBias"});
    auto MuonEta21 = d21.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi21 = d21.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonPt21 = d21.Histo1D({"MuonPt","",40,0.0,150},{"MuonPt"});
    auto leadingPt21 = d21.Histo1D({"Leading muon Pt","",40,0.0,150},{"leadingPt"});
    auto subPt21 = d21.Histo1D({"sub-leading muon Pt","",40,0.0,150},{"subPt"});
    auto MuonSize21 = d21.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass21 = d21.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});
    auto Z_match21 = d21.Histo1D({"Z_match","",150,0.0,150},{"Z_match"});
    auto MinBias21 = d21.Histo1D({"MinBias","",150,0.0,150},{"MinBias"});

    ROOT::RDataFrame d22(treeName, fileName22, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta22 = d22.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi22 = d22.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize22 = d22.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass22 = d22.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d23(treeName, fileName23, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta23 = d23.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi23 = d23.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize23 = d23.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass23 = d23.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d24(treeName, fileName24, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta24 = d24.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi24 = d24.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonPt24 = d24.Histo1D({"MuonPt","",40,0.0,150},{"MuonPt"});
    auto leadingPt24 = d24.Histo1D({"Leading muon Pt","",40,0.0,150},{"leadingPt"});
    auto subPt24 = d24.Histo1D({"sub-leading muon Pt","",40,0.0,150},{"subPt"});
    auto MuonSize24 = d24.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto HiggsMass = d24.Histo1D({"Z_mass","",150,0.0,150},{"Z_mass"});

    ROOT::RDataFrame d25(treeName, fileName25, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta25 = d25.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi25 = d25.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize25 = d25.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass25 = d25.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d26(treeName, fileName26, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta26 = d26.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi26 = d26.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize26 = d26.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass26 = d26.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d27(treeName, fileName27, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta27 = d27.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi27 = d27.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize27 = d27.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass27 = d27.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d28(treeName, fileName28, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta28 = d28.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi28 = d28.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize28 = d28.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass28 = d28.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d29(treeName, fileName29, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta29 = d29.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi29 = d29.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize29 = d29.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass29 = d29.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d30(treeName, fileName30, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta30 = d30.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi30 = d30.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize30 = d30.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass30 = d30.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d31(treeName, fileName31, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta31 = d31.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi31 = d31.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize31 = d31.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass31 = d31.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d32(treeName, fileName32, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta32 = d32.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi32 = d32.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize32 = d32.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass32 = d32.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d33(treeName, fileName33, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta33 = d33.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi33 = d33.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize33 = d33.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass33 = d33.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d34(treeName, fileName34, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta34 = d34.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi34 = d34.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize34 = d34.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass34 = d34.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d35(treeName, fileName35, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta35 = d35.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi35 = d35.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize35 = d35.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass35 = d35.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    ROOT::RDataFrame d36(treeName, fileName36, {"MuonEta", "MuonPhi", "MuonPt", "MuonE", "leadingPt", "subPt", "muonSize",
    "Z_mass"});
    auto MuonEta36 = d36.Histo1D({"MuonEta","",20,-2.7,2.7},{"MuonEta"});
    auto MuonPhi36 = d36.Histo1D({"MuonPhi","",20,-3.14,3.14},{"MuonPhi"});
    auto MuonSize36 = d36.Histo1D({"# of muons","",10,-0.5,9.5},{"muonSize"});
    auto Z_mass36 = d36.Histo1D({"Z_mass","",150,0,150},{"Z_mass"});

    // Draw MuonEta
    auto c1 = new TCanvas("c1", "", 200, 10, 800, 600);
    MuonEta->GetXaxis()->SetTitle("#eta");
    MuonEta->GetYaxis()->SetTitle("N/bin");
    MuonEta->GetYaxis()->SetRangeUser(0.0, 2500);
    MuonEta->SetMarkerSize(1.);
    MuonEta->SetMarkerStyle(22.);
    MuonEta->SetMarkerColor(kBlack);
    MuonEta->SetLineColor(kBlack);
    MuonEta->Draw("E0");
    //MuonEta21->Draw("same");
    gPad->RedrawAxis();
    c1->Print("MuonEta.pdf");

    // Draw MuonPhi
    auto c2 = new TCanvas("c2", "", 200, 10, 800, 600);
    MuonPhi->GetXaxis()->SetTitle("#phi");
    MuonPhi->GetYaxis()->SetTitle("N/bin");
    MuonPhi->GetYaxis()->SetRangeUser(0.0, 2000);
    MuonPhi->SetMarkerSize(1.);
    MuonPhi->SetMarkerStyle(22.);
    MuonPhi->SetMarkerColor(kBlack);
    MuonPhi->SetLineColor(kBlack);
    MuonPhi->Draw("E0");
    gPad->RedrawAxis();
    c2->Print("MuonPhi.pdf");

    // Draw MuonPt
    auto c3 = new TCanvas("c3", "", 200, 10, 800, 600);
    MuonPt->GetXaxis()->SetTitle("Pt [GeV]");
    MuonPt->GetYaxis()->SetTitle("N/bin");
    MuonPt->SetMarkerSize(1.);
    MuonPt->SetMarkerStyle(22.);
    MuonPt->SetMarkerColor(kBlack);
    MuonPt->SetLineColor(kBlack);
    MuonPt->Draw("E0");
    //MuonPt21->Draw("same");
    gPad->RedrawAxis();
    c3->Print("MuonPt.pdf");

    // Draw MuonE
    auto c4 = new TCanvas("c4", "", 200, 10, 800, 600);
    MuonE->GetXaxis()->SetTitle("E [GeV]");
    MuonE->GetYaxis()->SetTitle("N/bin");
    MuonE->SetMarkerSize(1.);
    MuonE->SetMarkerStyle(22.);
    MuonE->SetMarkerColor(kBlack);
    MuonE->SetLineColor(kBlack);
    MuonE->Draw("E0");
    gPad->RedrawAxis();
    c4->Print("MuonE.pdf");


    // Draw leading muon Pt
    auto c5 = new TCanvas("c5", "", 200, 10, 800, 600);
    leadingPt->GetXaxis()->SetTitle("leading Pt [GeV]");
    leadingPt->GetYaxis()->SetTitle("N/bin");
    leadingPt->SetMarkerSize(1.);
    leadingPt->SetMarkerStyle(22.);
    leadingPt->SetMarkerColor(kBlack);
    leadingPt->SetLineColor(kBlack);
    leadingPt->Draw("E0");
    gPad->RedrawAxis();
    c5->Print("leadingPt.pdf");

    // Draw sub-leading Pt
    auto c6 = new TCanvas("c6", "", 200, 10, 800, 600);
    subPt->GetXaxis()->SetTitle("sub Pt [GeV]");
    subPt->GetYaxis()->SetTitle("N/bin");
    subPt->SetMarkerSize(1.);
    subPt->SetMarkerStyle(22.);
    subPt->SetMarkerColor(kBlack);
    subPt->SetLineColor(kBlack);
    subPt->Draw("E0");
    gPad->RedrawAxis();
    c6->Print("subPt.pdf");

    // Draw muonSize
    auto c7 = new TCanvas("c7", "", 200, 10, 800, 600);
    MuonSize->GetXaxis()->SetTitle("# of #mu");
    MuonSize->GetYaxis()->SetTitle("N/bin");
    MuonSize->SetMarkerSize(1.);
    MuonSize->SetMarkerStyle(22.);
    MuonSize->SetMarkerColor(kBlack);
    MuonSize->SetLineColor(kBlack);                                                               
    MuonSize->Draw("E0");
    gPad->RedrawAxis();
    c7->Print("MuonSize.pdf");


    // Draw MuonScatter
    auto c8 = new TCanvas("c8", "", 200, 10, 800, 600);
    c8->SetRightMargin(0.18);
    MuonScatter->GetXaxis()->SetTitle("leadingPt [GeV]");
    MuonScatter->GetYaxis()->SetTitle("subPt [GeV]");
    MuonScatter->GetXaxis()->SetRangeUser(0.0, 150);
    MuonScatter->GetYaxis()->SetRangeUser(0.0, 150);   
    MuonScatter->SetMarkerSize(1.);
    MuonScatter->SetMarkerStyle(22.);
    MuonScatter->SetMarkerColor(kBlack);
    MuonScatter->SetLineColor(kBlack);
    MuonScatter->Draw("colz");
    gPad->RedrawAxis();
    c8->Print("MuonScatter.pdf");

    // Draw Z mass
    auto c51 = new TCanvas("c51", "", 200, 10, 800, 600);
    Z_mass->GetXaxis()->SetTitle("Invariant Z mass [GeV]");
    Z_mass->GetYaxis()->SetTitle("N/bin");
    Z_mass->SetMarkerSize(1.);
    Z_mass->SetMarkerStyle(22.);
    Z_mass->SetMarkerColor(kBlack);
    Z_mass->SetLineColor(kBlack);
    Z_mass->Draw("E0");
    //Z_mass->Fit("gaus");
    //Z_mass->Scale(1/MuonSize->GetEntries());
    //Z_mass->Scale(1.9);
    TF1 *fitB6 = new TF1("fitB6",bwconv,70,110,4);
    fitB6->SetLineColor(kBlue);
    fitB6->SetParameters(91,6,7000,3);
    fitB6->SetParNames("Mass", "fwhm", "Area", "Sigma");
    //fitB6->SetParLimits(4, 91, 92);
    //fitB6->SetParLimits(5, 0, 7);
    Z_mass->Fit("fitB6","+","",70,110);
    gPad->RedrawAxis();
    c51->Print("Z_massOrg.pdf");


    // Draw Z mass
    auto c10 = new TCanvas("c10", "", 200, 10, 800, 600);
    Z_mass1->GetXaxis()->SetTitle("Ztautau [GeV]");
    Z_mass1->GetYaxis()->SetTitle("N/bin");
    Z_mass1->SetMarkerSize(1.);
    Z_mass1->SetMarkerStyle(22.);
    Z_mass1->SetMarkerColor(kBlack);
    Z_mass1->SetLineColor(kBlack);
    Z_mass1->Draw();
    Z_mass1->Scale(1/MuonSize1->GetEntries());
    Z_mass1->Scale(1.9);
    Z_mass1->Scale(542);
    cout << "Efficiency of Z-tautau: " << Z_mass1->GetEntries()/MuonSize1->GetEntries() << endl;
    cout << "number of Z-tau events: " << MuonSize1->GetEntries() << endl;
    gPad->RedrawAxis();
    c10->Print("Z_mass1.pdf");

    // Draw Z mass
    auto c11 = new TCanvas("c11", "", 200, 10, 800, 600);
    Z_mass2->GetXaxis()->SetTitle("Diboson WZlvvv [GeV]");
    Z_mass2->GetYaxis()->SetTitle("N/bin");
    Z_mass2->SetMarkerSize(1.);
    Z_mass2->SetMarkerStyle(22.);
    Z_mass2->SetMarkerColor(kBlack);
    Z_mass2->SetLineColor(kBlack);
    Z_mass2->Draw("E0");
    Z_mass2->Scale(1/MuonSize2->GetEntries());
    Z_mass2->Scale(0.00278);
    cout << "Efficiency of Diboson WZlvvv: " << Z_mass2->GetEntries()/MuonSize2->GetEntries() << endl;
    cout << "number of Diboson WZlvvv events: " << MuonSize2->GetEntries() << endl;
    gPad->RedrawAxis();
    c11->Print("Z_mass2.pdf");

    // Draw Z mass
    auto c25 = new TCanvas("c25", "", 200, 10, 800, 600);
    Z_mass12->GetXaxis()->SetTitle("Diboson WWlvlv [GeV]");
    Z_mass12->GetYaxis()->SetTitle("N/bin");
    Z_mass12->SetMarkerSize(1.);
    Z_mass12->SetMarkerStyle(22.);
    Z_mass12->SetMarkerColor(kBlack);
    Z_mass12->SetLineColor(kBlack);
    Z_mass12->Draw("E0");
    Z_mass12->Scale(1/MuonSize12->GetEntries());
    Z_mass12->Scale(0.0106);
    cout << "Efficiency of Diboson WWlvlv: " << Z_mass12->GetEntries()/MuonSize12->GetEntries() << endl;
    cout << "number of WWlvlv events: " << MuonSize12->GetEntries() << endl;
    gPad->RedrawAxis();
    c25->Print("Z_mass12.pdf");

    // Draw Z mass
    auto c26 = new TCanvas("c26", "", 200, 10, 800, 600);
    Z_mass13->GetXaxis()->SetTitle("Diboson WZlvll [GeV]");
    Z_mass13->GetYaxis()->SetTitle("N/bin");
    Z_mass13->SetMarkerSize(1.);
    Z_mass13->SetMarkerStyle(22.);
    Z_mass13->SetMarkerColor(kBlack);
    Z_mass13->SetLineColor(kBlack);
    Z_mass13->Draw("E0");
    Z_mass13->Scale(1/MuonSize13->GetEntries());
    Z_mass13->Scale(0.00451);
    cout << "Efficiency of Diboson WZlvll: " << Z_mass13->GetEntries()/MuonSize13->GetEntries() << endl;
    cout << "number of WZlvll events: " << MuonSize13->GetEntries() << endl;
    gPad->RedrawAxis();
    c26->Print("Z_mass13.pdf");

    // Draw Z mass
    auto c27 = new TCanvas("c27", "", 200, 10, 800, 600);
    Z_mass14->GetXaxis()->SetTitle("Diboson ZZllll [GeV]");
    Z_mass14->GetYaxis()->SetTitle("N/bin");
    Z_mass14->SetMarkerSize(1.);
    Z_mass14->SetMarkerStyle(22.);
    Z_mass14->SetMarkerColor(kBlack);
    Z_mass14->SetLineColor(kBlack);
    Z_mass14->Draw("E0");
    Z_mass14->Scale(1/MuonSize14->GetEntries());
    Z_mass14->Scale(0.00127);
    cout << "Efficiency of Diboson ZZllll: " << Z_mass14->GetEntries()/MuonSize14->GetEntries() << endl;
    cout << "number of ZZllll events: " << MuonSize14->GetEntries() << endl;
    gPad->RedrawAxis();
    c27->Print("Z_mass14.pdf");

    // Draw Z mass
    auto c28 = new TCanvas("c28", "", 200, 10, 800, 600);
    TH1D* Diboson = new TH1D ("Diboson","",150,0,150);
    Z_mass15->Scale(1/MuonSize15->GetEntries());
    Z_mass15->Scale(0.000549);
    Diboson = (TH1D*) Z_mass15->Clone();
    Diboson->GetXaxis()->SetTitle("Diboson [GeV]");
    Diboson->GetYaxis()->SetTitle("N/bin");
    Diboson->SetMarkerSize(1.);
    Diboson->SetMarkerStyle(22.);
    Diboson->SetMarkerColor(kBlack);
    Diboson->SetLineColor(kBlack);
    Diboson->Add(Z_mass2.GetPtr());
    Diboson->Add(Z_mass12.GetPtr());
    Diboson->Add(Z_mass13.GetPtr());
    Diboson->Add(Z_mass14.GetPtr());
    Diboson->Scale(542);
    Diboson->Draw();
    cout << "Efficiency of Diboson ZZvvvv: " << Z_mass15->GetEntries()/MuonSize15->GetEntries() << endl;
    cout << "number of ZZvvvv events: " << MuonSize15->GetEntries() << endl;
    gPad->RedrawAxis();
    c28->Print("Diboson.pdf");

    // Draw Z mass
    auto c12 = new TCanvas("c12", "", 200, 10, 800, 600);
    Z_mass3->GetXaxis()->SetTitle("Single top [GeV]");
    Z_mass3->GetYaxis()->SetTitle("N/bin");
    Z_mass3->SetMarkerSize(1.);
    Z_mass3->SetMarkerStyle(22.);
    Z_mass3->SetMarkerColor(kBlack);
    Z_mass3->SetLineColor(kBlack);
    Z_mass3->Draw();
    Z_mass3->Scale(1/MuonSize3->GetEntries());
    Z_mass3->Scale(0.0437);
    Z_mass3->Scale(542);
    cout << "Efficiency of Single top: " << Z_mass3->GetEntries()/MuonSize3->GetEntries() << endl;
    cout << "number of Single top events: " << MuonSize3->GetEntries() << endl;
    gPad->RedrawAxis();
    c12->Print("Z_mass3.pdf");

    // Draw Z mass
    auto c13 = new TCanvas("c13", "", 200, 10, 800, 600);
    Z_mass4->GetXaxis()->SetTitle("ttbar [GeV]");
    Z_mass4->GetYaxis()->SetTitle("N/bin");
    Z_mass4->SetMarkerSize(1.);
    Z_mass4->SetMarkerStyle(22.);
    Z_mass4->SetMarkerColor(kBlack);
    Z_mass4->SetLineColor(kBlack);
    Z_mass4->Draw();
    Z_mass4->Scale(1/MuonSize4->GetEntries());
    Z_mass4->Scale(0.696);
    Z_mass4->Scale(542);
    cout << "Efficiency of ttbar: " << Z_mass4->GetEntries()/MuonSize4->GetEntries() << endl;
    cout << "number of ttbar events: " << MuonSize4->GetEntries() << endl;
    gPad->RedrawAxis();
    c13->Print("Z_mass4.pdf");

    // Draw Z mass
    auto c14 = new TCanvas("c14", "", 200, 10, 800, 600);
    Z_mass5->GetXaxis()->SetTitle("Drell-Yan 120-180 [GeV]");
    Z_mass5->GetYaxis()->SetTitle("N/bin");
    Z_mass5->SetMarkerSize(1.);
    Z_mass5->SetMarkerStyle(22.);
    Z_mass5->SetMarkerColor(kBlack);
    Z_mass5->SetLineColor(kBlack);
    Z_mass5->Draw("E0");
    Z_mass5->Scale(1/MuonSize5->GetEntries());
    Z_mass5->Scale(0.0175);
    cout << "Efficiency of DY 120-180: " << Z_mass5->GetEntries()/MuonSize5->GetEntries() << endl;
    cout << "number of DY120-180 events: " << MuonSize5->GetEntries() << endl;
    gPad->RedrawAxis();
    c14->Print("Z_mass5.pdf");

    // Draw Z mass
    auto c15 = new TCanvas("c15", "", 200, 10, 800, 600);
    Z_mass6->GetXaxis()->SetTitle("Drell-Yan 180-250 [GeV]");
    Z_mass6->GetYaxis()->SetTitle("N/bin");
    Z_mass6->SetMarkerSize(1.);
    Z_mass6->SetMarkerStyle(22.);
    Z_mass6->SetMarkerColor(kBlack);
    Z_mass6->SetLineColor(kBlack);
    Z_mass6->Draw("E0");
    Z_mass6->Scale(1/MuonSize6->GetEntries());
    Z_mass6->Scale(0.00292);
    cout << "Efficiency of DY 180-250: " << Z_mass6->GetEntries()/MuonSize6->GetEntries() << endl;
    cout << "number of DY180-250 events: " << MuonSize6->GetEntries() << endl;
    gPad->RedrawAxis();
    c15->Print("Z_mass6.pdf");

    // Draw Z mass
    auto c16 = new TCanvas("c16", "", 200, 10, 800, 600);
    Z_mass7->GetXaxis()->SetTitle("Drell-Yan 250-400 [GeV]");
    Z_mass7->GetYaxis()->SetTitle("N/bin");
    Z_mass7->SetMarkerSize(1.);
    Z_mass7->SetMarkerStyle(22.);
    Z_mass7->SetMarkerColor(kBlack);
    Z_mass7->SetLineColor(kBlack);
    Z_mass7->Draw("E0");
    Z_mass7->Scale(1/MuonSize7->GetEntries());
    Z_mass7->Scale(0.00108);
    cout << "Efficiency of DY 250-400: " << Z_mass7->GetEntries()/MuonSize7->GetEntries() << endl;
    cout << "number of DY250-400 events: " << MuonSize7->GetEntries() << endl;
    gPad->RedrawAxis();
    c16->Print("Z_mass7.pdf");

    // Draw Z mass
    auto c17 = new TCanvas("c17", "", 200, 10, 800, 600);
    Z_mass8->GetXaxis()->SetTitle("Drell-Yan 400-600 [GeV]");
    Z_mass8->GetYaxis()->SetTitle("N/bin");
    Z_mass8->SetMarkerSize(1.);
    Z_mass8->SetMarkerStyle(22.);
    Z_mass8->SetMarkerColor(kBlack);
    Z_mass8->SetLineColor(kBlack);
    Z_mass8->Draw("E0");
    Z_mass8->Scale(1/MuonSize8->GetEntries());
    Z_mass8->Scale(0.000196);
    cout << "Efficiency of DY 400-600: " << Z_mass8->GetEntries()/MuonSize8->GetEntries() << endl;
    cout << "number of DY400-600 events: " << MuonSize8->GetEntries() << endl;
    gPad->RedrawAxis();
    c17->Print("Z_mass8.pdf");

    // Draw Z mass
    auto c18 = new TCanvas("c18", "", 200, 10, 800, 600);
    Z_mass9->GetXaxis()->SetTitle("Drell-Yan 800-1000 [GeV]");
    Z_mass9->GetYaxis()->SetTitle("N/bin");
    Z_mass9->SetMarkerSize(1.);
    Z_mass9->SetMarkerStyle(22.);
    Z_mass9->SetMarkerColor(kBlack);
    Z_mass9->SetLineColor(kBlack);
    Z_mass9->Draw("E0");
    Z_mass9->Scale(1/MuonSize9->GetEntries());
    Z_mass9->Scale(0.0000106);
    cout << "Efficiency of DY 800-1000: " << Z_mass9->GetEntries()/MuonSize9->GetEntries() << endl;
    cout << "number of DY800-1000 events: " << MuonSize9->GetEntries() << endl;
    gPad->RedrawAxis();
    c18->Print("Z_mass9.pdf");

    // Draw Z mass
    auto c24 = new TCanvas("c24", "", 200, 10, 800, 600);
    Z_mass11->GetXaxis()->SetTitle("Drell-Yan 600-800 [GeV]");
    Z_mass11->GetYaxis()->SetTitle("N/bin");
    Z_mass11->SetMarkerSize(1.);
    Z_mass11->SetMarkerStyle(22.);
    Z_mass11->SetMarkerColor(kBlack);
    Z_mass11->SetLineColor(kBlack);
    Z_mass11->Draw("E0");
    Z_mass11->Scale(1/MuonSize11->GetEntries());
    Z_mass11->Scale(0.0000374);
    cout << "Efficiency of DY 600-800: " << Z_mass11->GetEntries()/MuonSize11->GetEntries() << endl;
    cout << "number of DY600-800 events: " << MuonSize11->GetEntries() << endl;
    gPad->RedrawAxis();
    c24->Print("Z_mass11.pdf");

    // Draw Z mass
    auto c19 = new TCanvas("c19", "", 200, 10, 800, 600);
    TH1D* DrellYan = new TH1D ("DY","",150,0,150);
    Z_mass10->Scale(1/MuonSize10->GetEntries());
    Z_mass10->Scale(0.00000426);
    DrellYan = (TH1D*) Z_mass10->Clone();
    DrellYan->GetXaxis()->SetTitle("Drell-Yan [GeV]");
    DrellYan->GetYaxis()->SetTitle("N/bin");
    DrellYan->SetMarkerSize(1.);
    DrellYan->SetMarkerStyle(22.);
    DrellYan->SetMarkerColor(kBlack);
    DrellYan->SetLineColor(kBlack);
    DrellYan->Add(Z_mass5.GetPtr());
    DrellYan->Add(Z_mass6.GetPtr());
    DrellYan->Add(Z_mass7.GetPtr());
    DrellYan->Add(Z_mass8.GetPtr());
    DrellYan->Add(Z_mass9.GetPtr());
    DrellYan->Add(Z_mass11.GetPtr());
    DrellYan->Scale(542);
    DrellYan->Draw();
    cout << "Efficiency of DY 1000-1250: " << Z_mass10->GetEntries()/MuonSize10->GetEntries() << endl;
    cout << "number of DY1000-1250 events: " << MuonSize10->GetEntries() << endl;
    gPad->RedrawAxis();
    c19->Print("DrellYan.pdf");

    // Draw Z mass
    auto c29 = new TCanvas("c29", "", 200, 10, 800, 600);
    Z_mass16->GetXaxis()->SetTitle("Dijet JF17 [GeV]");
    Z_mass16->GetYaxis()->SetTitle("N/bin");
    Z_mass16->SetMarkerSize(1.);
    Z_mass16->SetMarkerStyle(22.);
    Z_mass16->SetMarkerColor(kBlack);
    Z_mass16->SetLineColor(kBlack);
    Z_mass16->Draw("E0");
    Z_mass16->Scale(1/MuonSize16->GetEntries());
    Z_mass16->Scale(2430000);
    cout << "Efficiency of Dijet JF17: " << Z_mass16->GetEntries()/MuonSize16->GetEntries() << endl;
    cout << "number of JF17 events: " << MuonSize16->GetEntries() << endl;
    gPad->RedrawAxis();
    c29->Print("Z_mass16.pdf");

    // Draw Z mass
    auto c30 = new TCanvas("c30", "", 200, 10, 800, 600);
    Z_mass17->GetXaxis()->SetTitle("Dijet JF23 [GeV]");
    Z_mass17->GetYaxis()->SetTitle("N/bin");
    Z_mass17->SetMarkerSize(1.);
    Z_mass17->SetMarkerStyle(22.);
    Z_mass17->SetMarkerColor(kBlack);
    Z_mass17->SetLineColor(kBlack);
    Z_mass17->Draw("E0");
    Z_mass17->Scale(1/MuonSize17->GetEntries());
    Z_mass17->Scale(728000);
    cout << "Efficiency of Dijet JF23: " << Z_mass17->GetEntries()/MuonSize17->GetEntries() << endl;
    cout << "number of JF23 events: " << MuonSize17->GetEntries() << endl;
    gPad->RedrawAxis();
    c30->Print("Z_mass17.pdf");

    // Draw Z mass
    auto c31 = new TCanvas("c31", "", 200, 10, 800, 600);
    Z_mass18->GetXaxis()->SetTitle("Dijet JF35 [GeV]");
    Z_mass18->GetYaxis()->SetTitle("N/bin");
    Z_mass18->SetMarkerSize(1.);
    Z_mass18->SetMarkerStyle(22.);
    Z_mass18->SetMarkerColor(kBlack);
    Z_mass18->SetLineColor(kBlack);
    Z_mass18->Draw("E0");
    Z_mass18->Scale(1/MuonSize18->GetEntries());
    Z_mass18->Scale(134000);
    cout << "Efficiency of Dijet JF35: " << Z_mass18->GetEntries()/MuonSize18->GetEntries() << endl;
    cout << "number of JF35 events: " << MuonSize18->GetEntries() << endl;
    gPad->RedrawAxis();
    c31->Print("Z_mass18.pdf");

    // Draw Z mass
    auto c32 = new TCanvas("c32", "", 200, 10, 800, 600);
    TH1D* Dijet = new TH1D ("DJ","",150,0,150);
    Dijet = (TH1D*) Z_mass18->Clone();
    Dijet->GetXaxis()->SetTitle("Dijet [GeV]");
    Dijet->GetYaxis()->SetTitle("N/bin");
    Dijet->SetMarkerSize(1.);
    Dijet->SetMarkerStyle(22.);
    Dijet->SetMarkerColor(kBlack);
    Dijet->SetLineColor(kBlack);
    Dijet->Add(Z_mass16.GetPtr());
    Dijet->Add(Z_mass17.GetPtr());
    Dijet->Draw("E0");
    gPad->RedrawAxis();
    c32->Print("Dijet.pdf");

    // Draw Z_match
    auto c21 = new TCanvas("c21", "", 200, 10, 800, 600);
    Z_match->GetXaxis()->SetTitle("matched Z mass [GeV]");
    Z_match->GetYaxis()->SetTitle("N/bin");
    Z_match->SetMarkerSize(1.);
    Z_match->SetMarkerStyle(22.);
    Z_match->SetMarkerColor(kBlack);
    Z_match->SetLineColor(kBlack);
    Z_match->Draw("E0");
    //Z_match->Fit("gaus");
    TF1 *fitB5 = new TF1("fitB5",bwconv,70,110,4);
    fitB5->SetLineColor(kBlue);
    fitB5->SetParameters(91,6,5000,3);
    fitB5->SetParNames("Mass", "fwhm", "Area", "Sigma");
    fitB5->FixParameter(0, 91.1351);
    fitB5->FixParameter(1,2.55158);
    Z_match->Fit("fitB5","+","",70,110);
    gPad->RedrawAxis();
    c21->Print("Z_match.pdf");

    // Draw min-bias
    auto c22 = new TCanvas("c22", "", 200, 10, 800, 600);
    MinBias->GetXaxis()->SetTitle("Min Bias [GeV]");
    MinBias->GetYaxis()->SetTitle("N/bin");
    MinBias->SetMarkerSize(1.);
    MinBias->SetMarkerStyle(22.);
    MinBias->SetMarkerColor(kBlack);
    MinBias->SetLineColor(kBlack);
    MinBias->Scale(1/MuonSize->GetEntries());
    //MinBias->Scale(78400000);
    MinBias->Draw();
    gPad->RedrawAxis();
    c22->Print("MinBias.pdf");

    // Draw Z_mass vs Z_match
    //auto c23 = new TCanvas("c23", "", 200, 10, 800, 600);
    //Z_mass->GetXaxis()->SetTitle("Z reco Vs truth match [GeV]");
    //Z_mass->GetYaxis()->SetTitle("N/bin");
    //Z_mass->SetMarkerSize(1.);
    //Z_mass->SetMarkerStyle(22.);
    //Z_mass->SetMarkerColor(kBlack);
    //Z_mass->SetLineColor(kBlack);
    //Z_mass->Draw("E0");
    //Z_match->Draw("same p plc pmc");
    //gPad->RedrawAxis();
    //c23->Print("recoTM.pdf");

    // Draw Z_match
    auto c42 = new TCanvas("c42", "", 200, 10, 800, 600);
    Z_match21->GetXaxis()->SetTitle("matched Z mass [GeV]");
    Z_match21->GetYaxis()->SetTitle("N/bin");
    Z_match21->SetMarkerSize(1.);
    Z_match21->SetMarkerStyle(22.);
    Z_match21->SetMarkerColor(kBlack);
    Z_match21->SetLineColor(kBlack);
    Z_match21->Draw();
    //Z_match21->Fit("gaus");
    gPad->RedrawAxis();
    c42->Print("Z_matchSelect.pdf");

    // Draw min-bias
    auto c43 = new TCanvas("c43", "", 200, 10, 800, 600);
    MinBias21->GetXaxis()->SetTitle("Min Bias [GeV]");
    MinBias21->GetYaxis()->SetTitle("N/bin");
    MinBias21->SetMarkerSize(1.);
    MinBias21->SetMarkerStyle(22.);
    MinBias21->SetMarkerColor(kBlack);
    MinBias21->SetLineColor(kBlack);
    //MinBias21->Scale(1/MuonSize->GetEntries());
    //MinBias21->Scale(78400000);
    MinBias21->Draw();
    gPad->RedrawAxis();
    c43->Print("MinBiasSelect.pdf");

    // Draw Z mass
    auto c50 = new TCanvas("c50", "", 200, 10, 800, 600);
    Z_mass21->GetXaxis()->SetTitle("Invariant Z mass [GeV]");
    Z_mass21->GetYaxis()->SetTitle("N/bin");
    Z_mass21->SetMarkerSize(1.);
    Z_mass21->SetMarkerStyle(22.);
    Z_mass21->SetMarkerColor(kBlack);
    Z_mass21->SetLineColor(kBlack);
    Z_mass21->Draw("E0");
    //Z_mass21->Fit("gaus");
    cout << "Efficiency of Z-mumu: " << Z_mass21->GetEntries()/MuonSize21->GetEntries() << endl;
    gPad->RedrawAxis();
    c50->Print("Z_Select.pdf");

    // Draw Z mass
    auto c9 = new TCanvas("c9", "", 200, 10, 800, 600);
    TH1D* Zorg = new TH1D ("orgZ","",150,0,150);
    Z_mass21->Scale(1/MuonSize21->GetEntries());
    Z_mass21->Scale(1.9);
    Z_mass21->Scale(542);
    Zorg = (TH1D*) Z_mass21->Clone();
    Zorg->GetXaxis()->SetTitle("Invariant Z mass [GeV]");
    Zorg->GetYaxis()->SetTitle("N/bin");
    Zorg->GetYaxis()->SetRangeUser(0, 0.12);
    Zorg->Add(Z_mass1.GetPtr());
    Zorg->Add(Diboson);
    Zorg->Add(Z_mass3.GetPtr());
    Zorg->Add(Z_mass4.GetPtr());
    Zorg->Add(DrellYan);
    //Zorg->Add(MinBias.GetPtr());
    Zorg->SetMarkerSize(1.);
    Zorg->SetMarkerStyle(22.);
    Zorg->SetMarkerColor(kBlack);
    gPad->Update();
    //Zorg->Fit("pol3");
    Zorg->Draw();
    Z_mass1->Draw("same plc pmc");
    Diboson->Draw("same plc pmc");
    Z_mass3->Draw("same plc pmc");
    Z_mass4->Draw("same plc pmc");
    DrellYan->Draw("same plc pmc");
    //MinBias->Draw("same plc pmc");
    Double_t signal = 0;
    Double_t background = 0;
    for (Int_t i=82;i<=100;i++) {
      signal = signal + Zorg->GetBinContent(i);
      background = background + Zorg->GetBinContent(i) - Z_mass21->GetBinContent(i); }
    Double_t significance = signal / sqrt(signal + background);
    cout << "signal: " << signal << endl;
    cout << "background: " << background << endl;
    cout << "significance: " << significance << endl;
    //gPad->SetLogy();
    gPad->RedrawAxis();
    Double_t Mass = 91.0;
    Double_t fwhm = 6.0;
    Double_t IntMin = 60.0;
    Double_t IntMax = 110.0;
    Double_t epsrel = 1.e-12;
    TF1 *fitB = new TF1("fitB",bwconv,60,110,4);
    fitB->SetLineColor(kBlue);
    fitB->SetParameters(91,6,7000,3);
    fitB->SetParNames("Mass", "fwhm", "Area", "Sigma");
    Zorg->Fit("fitB","+","",60,110);
    Double_t AreaZ = fitB->Integral(IntMin, IntMax, epsrel);
    cout << "Area under MC Z peak: " << AreaZ << endl;
    TF1 *fex = new TF1("fex",fitExp,0,150,2);
    fex->SetLineColor(kRed);
    fex->SetParameters(-4,-2);
    fex->SetParNames("const","slope");
    Zorg->Fit("fex","+","",0,150);
    Double_t HiggsBack = 0;
    for (Int_t i=116;i<=134;i++) {
      HiggsBack = HiggsBack + Zorg->GetBinContent(i); }
    cout << "Higgs background in MC: " << HiggsBack << endl;
    TLegend *legend = new TLegend(0.2,0.7,0.48,0.9);
    legend->AddEntry(Zorg,"Signal with background","p");
    legend->AddEntry(Z_mass1.GetPtr(),"Ztautau","p");
    legend->AddEntry(Z_mass4.GetPtr(),"ttbar","p");
    //legend->AddEntry(MinBias.GetPtr(),"Min Bias","p");
    legend->AddEntry(Z_mass3.GetPtr(),"SingleTop","p");
    legend->AddEntry(DrellYan,"Drell-Yan","p");
    legend->AddEntry(Diboson,"Diboson","p");
    legend->Draw();
    c9->Print("Z_mass.pdf");

    // Draw Z mass
    auto c33 = new TCanvas("c33", "", 200, 10, 800, 600);
    TH1D* FillPlot = new TH1D ("Filled Plot","",150,0,150);
    TH1D* Z_mass1_draw = new TH1D ("Z_mass1_draw","",150,0,150);
    for(int i = 1; i < Z_mass1->GetXaxis()->GetNbins() + 1; i++){
      Z_mass1_draw->SetBinContent(i, Z_mass1->GetBinContent(i));
    }
    Z_mass1_draw->SetFillColor(kBlue);
    FillPlot = (TH1D*) Z_mass1_draw->Clone();
    THStack hs("hs","test stacked histograms");
    hs.Add(Z_mass1_draw);
    TH1D* Diboson_draw = new TH1D ("Diboson_draw","",150,0,150);
    for(int i = 1; i < Diboson->GetXaxis()->GetNbins() + 1; i++){
      Diboson_draw->SetBinContent(i, Diboson->GetBinContent(i));
    }
    Diboson_draw->SetFillColor(kRed);
    FillPlot->Add(Diboson_draw);
    hs.Add(Diboson_draw);
    TH1D* Z_mass3_draw = new TH1D ("Z_mass3_draw","",150,0,150);
    for(int i = 1; i < Z_mass3->GetXaxis()->GetNbins() + 1; i++){
      Z_mass3_draw->SetBinContent(i, Z_mass3->GetBinContent(i));
    }
    Z_mass3_draw->SetFillColor(kGreen);
    FillPlot->Add(Z_mass3_draw);
    hs.Add(Z_mass3_draw);
    TH1D* Z_mass4_draw = new TH1D ("Z_mass4_draw","",150,0,150);
    for(int i = 1; i < Z_mass4->GetXaxis()->GetNbins() + 1; i++){
      Z_mass4_draw->SetBinContent(i, Z_mass4->GetBinContent(i));
    }
    Z_mass4_draw->SetFillColor(kYellow);
    FillPlot->Add(Z_mass4_draw);
    hs.Add(Z_mass4_draw);
    TH1D* DrellYan_draw = new TH1D ("DrellYan_draw","",150,0,150);
    for(int i = 1; i < DrellYan->GetXaxis()->GetNbins() + 1; i++){
      DrellYan_draw->SetBinContent(i, DrellYan->GetBinContent(i));
    }
    DrellYan_draw->SetFillColor(kOrange);
    FillPlot->Add(DrellYan_draw);
    hs.Add(DrellYan_draw);
    //TH1D* MinBias_draw = new TH1D ("MinBias_draw","",150,0,150);
    //for(int i = 1; i < MinBias->GetXaxis()->GetNbins() + 1; i++){
    //  MinBias_draw->SetBinContent(i, MinBias->GetBinContent(i));
    //}
    //MinBias_draw->SetFillColor(kMagenta);
    //FillPlot->Add(MinBias_draw);
    //hs.Add(MinBias_draw);
    TH1D* Z_mass21_draw = new TH1D ("Z_mass21_draw","",150,0,150);
    for(int i = 1; i < Z_mass21->GetXaxis()->GetNbins() + 1; i++){
      Z_mass21_draw->SetBinContent(i, Z_mass21->GetBinContent(i));
    }
    Z_mass21_draw->SetFillColor(kBlack);
    FillPlot->Add(Z_mass21_draw);
    hs.Add(Z_mass21_draw);
    Z_mass21_draw->GetXaxis()->SetTitle("Invariant Z mass [GeV]");
    Z_mass21_draw->GetYaxis()->SetRangeUser(0.0, 0.12);
    Z_mass21_draw->GetYaxis()->SetTitle("N/bin");
    gPad->Update();
    Z_mass21_draw->Draw();
    hs.Draw("same");
    //gPad->SetLogy();
    gPad->RedrawAxis();
    TLegend *legend2 = new TLegend(0.2,0.7,0.48,0.9);
    legend2->AddEntry(Z_mass21_draw,"Signal","f");
    legend2->AddEntry(Z_mass1_draw,"Ztautau","f");
    legend2->AddEntry(Z_mass4_draw,"ttbar","f");
    //legend2->AddEntry(MinBias_draw,"Min Bias","f");
    legend2->AddEntry(Z_mass3_draw,"SingleTop","f");
    legend2->AddEntry(DrellYan_draw,"Drell-Yan","f");
    legend2->AddEntry(Diboson_draw,"Diboson","f");
    legend2->Draw();
    c33->Print("FillPlot.pdf");

    // Draw MuonEta
    auto c34 = new TCanvas("c34", "", 200, 10, 800, 600);
    TH1D* etaComp = new TH1D ("etaComp","",40,-2.7,2.7);
    etaComp = (TH1D*) MuonEta20->Clone();
    etaComp->GetXaxis()->SetTitle("#eta");
    etaComp->GetYaxis()->SetTitle("N/bin");
    MuonEta20->GetYaxis()->SetRangeUser(0, 6500);
    etaComp->SetMarkerSize(1.);
    etaComp->SetMarkerStyle(22.);
    etaComp->SetMarkerColor(kBlack);
    etaComp->SetLineColor(kBlack);
    etaComp->Add(MuonEta22.GetPtr());
    etaComp->Add(MuonEta23.GetPtr());
    etaComp->Add(MuonEta25.GetPtr());
    etaComp->Add(MuonEta26.GetPtr());
    etaComp->Add(MuonEta27.GetPtr());
    etaComp->Add(MuonEta28.GetPtr());
    etaComp->Add(MuonEta29.GetPtr());
    etaComp->Add(MuonEta30.GetPtr());
    etaComp->Add(MuonEta31.GetPtr());
    etaComp->Add(MuonEta32.GetPtr());
    etaComp->Add(MuonEta33.GetPtr());
    etaComp->Add(MuonEta34.GetPtr());
    etaComp->Add(MuonEta35.GetPtr());
    etaComp->Add(MuonEta36.GetPtr());
    etaComp->Draw("E0");
    gPad->RedrawAxis();
    c34->Print("DataEtaComplete.pdf");

    // Draw MuonPhi
    auto c35 = new TCanvas("c35", "", 200, 10, 800, 600);
    TH1D* phiComp = new TH1D ("phiComp","",40,-3.14,3.14);
    phiComp = (TH1D*) MuonPhi20->Clone();
    phiComp->GetXaxis()->SetTitle("#phi");
    phiComp->GetYaxis()->SetTitle("N/bin");
    MuonPhi20->GetYaxis()->SetRangeUser(0, 5000);
    phiComp->SetMarkerSize(1.);
    phiComp->SetMarkerStyle(22.);
    phiComp->SetMarkerColor(kBlack);
    phiComp->SetLineColor(kBlack);
    phiComp->Add(MuonPhi22.GetPtr());
    phiComp->Add(MuonPhi23.GetPtr());
    phiComp->Add(MuonPhi25.GetPtr());
    phiComp->Add(MuonPhi26.GetPtr());
    phiComp->Add(MuonPhi27.GetPtr());
    phiComp->Add(MuonPhi28.GetPtr());
    phiComp->Add(MuonPhi29.GetPtr());
    phiComp->Add(MuonPhi30.GetPtr());
    phiComp->Add(MuonPhi31.GetPtr());
    phiComp->Add(MuonPhi32.GetPtr());
    phiComp->Add(MuonPhi33.GetPtr());
    phiComp->Add(MuonPhi34.GetPtr());
    phiComp->Add(MuonPhi35.GetPtr());
    phiComp->Add(MuonPhi36.GetPtr());
    phiComp->Draw("E0");
    gPad->RedrawAxis();
    c35->Print("DataPhiComplete.pdf");

    // Draw MuonPt
    auto c36 = new TCanvas("c36", "", 200, 10, 800, 600);
    MuonPt20->GetXaxis()->SetTitle("Pt [GeV]");
    MuonPt20->GetYaxis()->SetTitle("N/bin");
    MuonPt20->SetMarkerSize(1.);
    MuonPt20->SetMarkerStyle(22.);
    MuonPt20->SetMarkerColor(kBlack);
    MuonPt20->SetLineColor(kBlack);
    MuonPt20->Draw("E0");
    gPad->RedrawAxis();
    c36->Print("DataPt.pdf");

    // Draw leading muon Pt
    auto c37 = new TCanvas("c37", "", 200, 10, 800, 600);
    leadingPt20->GetXaxis()->SetTitle("leading Pt [GeV]");
    leadingPt20->GetYaxis()->SetTitle("N/bin");
    leadingPt20->SetMarkerSize(1.);
    leadingPt20->SetMarkerStyle(22.);
    leadingPt20->SetMarkerColor(kBlack);
    leadingPt20->SetLineColor(kBlack);
    leadingPt20->Draw("E0");
    gPad->RedrawAxis();
    c37->Print("leadingData.pdf");

    // Draw sub-leading Pt
    auto c38 = new TCanvas("c38", "", 200, 10, 800, 600);
    subPt20->GetXaxis()->SetTitle("sub Pt [GeV]");
    subPt20->GetYaxis()->SetTitle("N/bin");
    subPt20->SetMarkerSize(1.);
    subPt20->SetMarkerStyle(22.);
    subPt20->SetMarkerColor(kBlack);
    subPt20->SetLineColor(kBlack);
    subPt20->Draw("E0");
    gPad->RedrawAxis();
    c38->Print("subData.pdf");

    // Draw muonSize
    auto c39 = new TCanvas("c39", "", 200, 10, 800, 600);
    MuonSize20->GetXaxis()->SetTitle("# of #mu");
    MuonSize20->GetYaxis()->SetTitle("N/bin");
    MuonSize20->SetMarkerSize(1.);
    MuonSize20->SetMarkerStyle(22.);
    MuonSize20->SetMarkerColor(kBlack);
    MuonSize20->SetLineColor(kBlack);
    MuonSize20->Draw("E0");
    gPad->RedrawAxis();
    c39->Print("DataSize.pdf");

    // Draw Z mass
    auto c40 = new TCanvas("c40", "", 200, 10, 800, 600);
    Z_mass20->GetXaxis()->SetTitle("Z-mass [GeV]");
    Z_mass20->GetYaxis()->SetTitle("N/bin");
    Z_mass20->SetMarkerSize(1.);
    Z_mass20->SetMarkerStyle(22.);
    Z_mass20->SetMarkerColor(kBlack);
    Z_mass20->SetLineColor(kBlack);
    Z_mass20->Draw("E0");
    gPad->RedrawAxis();
    Double_t Zyield = 0;
    for (Int_t i=82;i<=100;i++) {
      Zyield = Zyield + Z_mass20->GetBinContent(i); }
    cout << "Z yield from data: " << Zyield << endl;
    // fit on data
    //Double_t fwhm1 = 6.0;
    TF1 *fitB1 = new TF1("fitB1",bwconv,60,110,4);
    fitB1->SetLineColor(kBlue);
    fitB1->SetParameters(91,6,500,3);
    fitB1->SetParNames("Mass", "fwhm", "Area", "Sigma");
    Z_mass20->Fit("fitB1","+","",60,110);
    Double_t AreaZd = fitB1->Integral(IntMin, IntMax, epsrel);
    cout << "Area under data Z peak: " << AreaZd << endl;
    c40->Print("Z_Data.pdf");

    // Draw Z mass
    auto c52 = new TCanvas("c52", "", 200, 10, 800, 600);
    TH1D* Data = new TH1D ("Data","",150,0,150);
    cout << "Number of data1 events: " << MuonSize20->GetEntries() << endl;
    cout << "Number of data1 reco's: " << Z_mass20->GetEntries() << endl;
    cout << "Number of data2 events: " << MuonSize22->GetEntries() << endl;
    cout << "Number of data2 reco's: " << Z_mass22->GetEntries() << endl;
    cout << "Number of data3 events: " << MuonSize23->GetEntries() << endl;
    cout << "Number of data3 reco's: " << Z_mass23->GetEntries() << endl;
    cout << "Number of data4 events: " << MuonSize25->GetEntries() << endl;
    cout << "Number of data4 reco's: " << Z_mass25->GetEntries() << endl;
    cout << "Number of data5 events: " << MuonSize26->GetEntries() << endl;
    cout << "Number of data5 reco's: " << Z_mass26->GetEntries() << endl;
    cout << "Number of data6 events: " << MuonSize27->GetEntries() << endl;
    cout << "Number of data6 reco's: " << Z_mass27->GetEntries() << endl;
    cout << "Number of data7 events: " << MuonSize28->GetEntries() << endl;
    cout << "Number of data7 reco's: " << Z_mass28->GetEntries() << endl;
    cout << "Number of data8 events: " << MuonSize29->GetEntries() << endl;
    cout << "Number of data8 reco's: " << Z_mass29->GetEntries() << endl;
    cout << "Number of data9 events: " << MuonSize30->GetEntries() << endl;
    cout << "Number of data9 reco's: " << Z_mass30->GetEntries() << endl;
    cout << "Number of data10 events: " << MuonSize31->GetEntries() << endl;
    cout << "Number of data10 reco's: " << Z_mass31->GetEntries() << endl;
    cout << "Number of data11 events: " << MuonSize32->GetEntries() << endl;
    cout << "Number of data11 reco's: " << Z_mass32->GetEntries() << endl;
    cout << "Number of data12 events: " << MuonSize33->GetEntries() << endl;
    cout << "Number of data12 reco's: " << Z_mass33->GetEntries() << endl;
    cout << "Number of data13 events: " << MuonSize34->GetEntries() << endl;
    cout << "Number of data13 reco's: " << Z_mass34->GetEntries() << endl;
    cout << "Number of data14 events: " << MuonSize35->GetEntries() << endl;
    cout << "Number of data14 reco's: " << Z_mass35->GetEntries() << endl;
    cout << "Number of data15 events: " << MuonSize36->GetEntries() << endl;
    cout << "Number of data15 reco's: " << Z_mass36->GetEntries() << endl;
    cout << "Total number of events ran over: " << MuonSize20->GetEntries() + MuonSize22->GetEntries() + MuonSize23->GetEntries() + MuonSize25->GetEntries() + MuonSize26->GetEntries() + MuonSize27->GetEntries() + MuonSize28->GetEntries() + MuonSize29->GetEntries() + MuonSize30->GetEntries() + MuonSize31->GetEntries() + MuonSize32->GetEntries() + MuonSize33->GetEntries() + MuonSize34->GetEntries() + MuonSize35->GetEntries() + MuonSize36->GetEntries() << endl;
    Data = (TH1D*) Z_mass22->Clone();
    Data->GetXaxis()->SetTitle("Invariant Mass [GeV]");
    Data->GetYaxis()->SetTitle("N/bin");
    Data->GetXaxis()->SetRangeUser(110, 150);
    Data->Add(Z_mass20.GetPtr());
    Data->Add(Z_mass23.GetPtr());
    Data->Add(Z_mass25.GetPtr());
    Data->Add(Z_mass26.GetPtr());
    Data->Add(Z_mass27.GetPtr());
    Data->Add(Z_mass28.GetPtr());
    Data->Add(Z_mass29.GetPtr());
    Data->Add(Z_mass30.GetPtr());
    Data->Add(Z_mass31.GetPtr());
    Data->Add(Z_mass32.GetPtr());
    Data->Add(Z_mass33.GetPtr());
    Data->Add(Z_mass34.GetPtr());
    Data->Add(Z_mass35.GetPtr());
    Data->Add(Z_mass36.GetPtr());
    Data->SetMarkerSize(1.);
    Data->SetMarkerStyle(22.);
    Data->SetMarkerColor(kBlack);
    Data->SetLineColor(kBlack);
    //Data->Fit("expo");
    Data->Draw("E0");
    gPad->RedrawAxis();
    Double_t Mass1 = 125;
    Double_t fwhm1 = 6.0;
    TF1 *fex1 = new TF1("fex1",fitExp,110,150,2);
    fex1->SetLineColor(kRed);
    fex1->SetParameters(10,-0.1);
    fex1->SetParNames("const","slope");
    fex1->FixParameter(0,7.75296);
    fex1->FixParameter(1,-0.0607259);
    Data->Fit("fex1","+","",110,150);
    Double_t IntMin1 = 116;
    Double_t IntMax1 = 134;
    Double_t AreaHBack = fex1->Integral(IntMin1, IntMax1, epsrel);
    cout << "background under Higgs peak: " << AreaHBack << endl;
    TF1 *fitComb = new TF1("fitComb",CombFit,110,150,8);
    fitComb->SetLineColor(kBlue);
    fitComb->SetParameters(Data->GetMaximum(),125,3,Data->GetRMS(),125,6,6.6,-0.03);
    fitComb->SetParNames("const","mass","sigma","const","mass","fwhm","const","slope");
    fitComb->SetParLimits(1,120,130);
    fitComb->SetParLimits(2,0,10);
    //fitComb->SetParLimits(3,0,10);
    //fitComb->FixParameter(1, 125);
    //fitComb->FixParameter(2, 2.99092);
    fitComb->FixParameter(3, 0);
    fitComb->FixParameter(4, 125);
    fitComb->FixParameter(5, 6);
    //fitComb->FixParameter(6, 7.69972);
    //fitComb->FixParameter(7, -0.048357);
    Data->Fit("fitComb","+","",110,150);
    Double_t AreaHD = fitComb->Integral(IntMin1, IntMax1, epsrel);
    cout << "Total area under Higgs peak: " << AreaHD << endl;
    cout << "area inside Higgs peak: " << AreaHBack - AreaHD << endl;
    Double_t HiggsYield = 0;
    for (Int_t i=82;i<=100;i++) {
      HiggsYield = HiggsYield + Data->GetBinContent(i); }
    cout << "Number of entries in Z combine region: " << HiggsYield << endl;
    TLegend *legend3 = new TLegend(0.2,0.7,0.48,0.9);
    legend3->AddEntry(fex1,"Background","l");
    legend3->AddEntry(fitComb,"Background + Signal","l");
    legend3->Draw();
    c52->Print("HiggsData.pdf");


    // Draw Z mass
    auto c41 = new TCanvas("c41", "", 200, 10, 800, 600);
    Z_mass20c->GetXaxis()->SetTitle("Z-mass [GeV]");
    Z_mass20c->GetYaxis()->SetTitle("N/bin");
    Z_mass20c->SetMarkerSize(1.);
    Z_mass20c->SetMarkerStyle(22.);
    Z_mass20c->SetMarkerColor(kBlack);
    Z_mass20c->SetLineColor(kBlack);
    Z_mass20c->Draw("E0");
    gPad->RedrawAxis();
    //TF1 *fitB2 = new TF1("fitB2",fitBoth,55,115,6);
    //fitB2->SetLineColor(kBlue);
    //fitB2->SetParameters(Z_mass20c->GetMaximum(),Mass,Z_mass20c->GetRMS(),500,Mass,fwhm);
    //fitB2->SetParNames("con","mass","sigma","con","mass","fwhm");
    //Z_mass20c->Fit("fitB2","+","",60,110);
    //Double_t AreaZdc = fitB2->Integral(IntMin, IntMax, epsrel);
    //cout << "Area under data Z peak: " << AreaZdc << endl;
    c41->Print("Z_DataCut.pdf");

    // Draw MuonEta
    auto c44 = new TCanvas("c44", "", 200, 10, 800, 600);
    MuonEta21->GetXaxis()->SetTitle("#eta");
    MuonEta21->GetYaxis()->SetTitle("N/bin");
    MuonEta21->GetYaxis()->SetRangeUser(0.0, 2500);
    MuonEta21->SetMarkerSize(1.);
    MuonEta21->SetMarkerStyle(22.);
    MuonEta21->SetMarkerColor(kBlack);
    MuonEta21->SetLineColor(kBlack);
    MuonEta21->Draw("E0");
    gPad->RedrawAxis();
    c44->Print("MuonEtaSelect.pdf");

    // Draw MuonPhi
    auto c45 = new TCanvas("c45", "", 200, 10, 800, 600);
    MuonPhi21->GetXaxis()->SetTitle("#phi");
    MuonPhi21->GetYaxis()->SetTitle("N/bin");
    MuonPhi21->GetYaxis()->SetRangeUser(0.0, 2000);
    MuonPhi21->SetMarkerSize(1.);
    MuonPhi21->SetMarkerStyle(22.);
    MuonPhi21->SetMarkerColor(kBlack);
    MuonPhi21->SetLineColor(kBlack);
    MuonPhi21->Draw("E0");
    gPad->RedrawAxis();
    c45->Print("MuonPhiSelect.pdf");

    // Draw MuonPt
    auto c46 = new TCanvas("c46", "", 200, 10, 800, 600);
    MuonPt21->GetXaxis()->SetTitle("Pt [GeV]");
    MuonPt21->GetYaxis()->SetTitle("N/bin");
    MuonPt21->SetMarkerSize(1.);
    MuonPt21->SetMarkerStyle(22.);
    MuonPt21->SetMarkerColor(kBlack);
    MuonPt21->SetLineColor(kBlack);
    MuonPt21->Draw("E0");
    gPad->RedrawAxis();
    c46->Print("MuonPtSelect.pdf");

    // Draw leading muon Pt
    auto c47 = new TCanvas("c47", "", 200, 10, 800, 600);
    leadingPt21->GetXaxis()->SetTitle("leading Pt [GeV]");
    leadingPt21->GetYaxis()->SetTitle("N/bin");
    leadingPt21->SetMarkerSize(1.);
    leadingPt21->SetMarkerStyle(22.);
    leadingPt21->SetMarkerColor(kBlack);
    leadingPt21->SetLineColor(kBlack);
    leadingPt21->Draw("E0");
    gPad->RedrawAxis();
    c47->Print("leadingSelect.pdf");

    // Draw sub-leading Pt
    auto c48 = new TCanvas("c48", "", 200, 10, 800, 600);
    subPt21->GetXaxis()->SetTitle("sub Pt [GeV]");
    subPt21->GetYaxis()->SetTitle("N/bin");
    subPt21->SetMarkerSize(1.);
    subPt21->SetMarkerStyle(22.);
    subPt21->SetMarkerColor(kBlack);
    subPt21->SetLineColor(kBlack);
    subPt21->Draw("E0");
    gPad->RedrawAxis();
    c48->Print("subSelect.pdf");

    // Draw muonSize
    auto c49 = new TCanvas("c49", "", 200, 10, 800, 600);
    MuonSize21->GetXaxis()->SetTitle("# of #mu");
    MuonSize21->GetYaxis()->SetTitle("N/bin");
    MuonSize21->SetMarkerSize(1.);
    MuonSize21->SetMarkerStyle(22.);
    MuonSize21->SetMarkerColor(kBlack);
    MuonSize21->SetLineColor(kBlack);
    MuonSize21->Draw("E0");
    //MuonSize->Draw("same p pmc");
    gPad->RedrawAxis();
    c49->Print("MuonSizeSelect.pdf");

    // Draw muonSize Higgs
    auto c53 = new TCanvas("c53", "", 200, 10, 800, 600);
    MuonSize24->GetXaxis()->SetTitle("# of #mu");
    MuonSize24->GetYaxis()->SetTitle("N/bin");
    MuonSize24->SetMarkerSize(1.);
    MuonSize24->SetMarkerStyle(22.);
    MuonSize24->SetMarkerColor(kBlack);
    MuonSize24->SetLineColor(kBlack);
    MuonSize24->Draw("E0");
    gPad->RedrawAxis();
    c53->Print("MuonSizeHiggs.pdf");

    // Draw MuonEta
    auto c54 = new TCanvas("c54", "", 200, 10, 800, 600);
    MuonEta24->GetXaxis()->SetTitle("#eta");
    MuonEta24->GetYaxis()->SetTitle("N/bin");
    MuonEta24->GetYaxis()->SetRangeUser(0.0, 2500);
    MuonEta24->SetMarkerSize(1.);
    MuonEta24->SetMarkerStyle(22.);
    MuonEta24->SetMarkerColor(kBlack);
    MuonEta24->SetLineColor(kBlack);
    MuonEta24->Draw("E0");
    gPad->RedrawAxis();
    c54->Print("MuonEtaHiggs.pdf");

    // Draw MuonPhi
    auto c55 = new TCanvas("c55", "", 200, 10, 800, 600);
    MuonPhi24->GetXaxis()->SetTitle("#phi");
    MuonPhi24->GetYaxis()->SetTitle("N/bin");
    MuonPhi24->GetYaxis()->SetRangeUser(0.0, 2000);
    MuonPhi24->SetMarkerSize(1.);
    MuonPhi24->SetMarkerStyle(22.);
    MuonPhi24->SetMarkerColor(kBlack);
    MuonPhi24->SetLineColor(kBlack);
    MuonPhi24->Draw("E0");
    gPad->RedrawAxis();
    c55->Print("MuonPhiHiggs.pdf");

    // Draw MuonPt
    auto c56 = new TCanvas("c56", "", 200, 10, 800, 600);
    MuonPt24->GetXaxis()->SetTitle("Pt [GeV]");
    MuonPt24->GetYaxis()->SetTitle("N/bin");
    MuonPt24->SetMarkerSize(1.);
    MuonPt24->SetMarkerStyle(22.);
    MuonPt24->SetMarkerColor(kBlack);
    MuonPt24->SetLineColor(kBlack);
    MuonPt24->Draw("E0");
    gPad->RedrawAxis();
    c56->Print("MuonPtHiggs.pdf");

    // Draw leading muon Pt
    auto c57 = new TCanvas("c57", "", 200, 10, 800, 600);
    leadingPt24->GetXaxis()->SetTitle("leading Pt [GeV]");
    leadingPt24->GetYaxis()->SetTitle("N/bin");
    leadingPt24->SetMarkerSize(1.);
    leadingPt24->SetMarkerStyle(22.);
    leadingPt24->SetMarkerColor(kBlack);
    leadingPt24->SetLineColor(kBlack);
    leadingPt24->Draw("E0");
    gPad->RedrawAxis();
    c57->Print("leadingHiggs.pdf");

    // Draw sub-leading Pt
    auto c58 = new TCanvas("c58", "", 200, 10, 800, 600);
    subPt24->GetXaxis()->SetTitle("sub Pt [GeV]");
    subPt24->GetYaxis()->SetTitle("N/bin");
    subPt24->SetMarkerSize(1.);
    subPt24->SetMarkerStyle(22.);
    subPt24->SetMarkerColor(kBlack);
    subPt24->SetLineColor(kBlack);
    subPt24->Draw("E0");
    gPad->RedrawAxis();
    c58->Print("subHiggs.pdf");

    // Draw Z mass
    auto c59 = new TCanvas("c59", "", 200, 10, 800, 600);
    HiggsMass->GetXaxis()->SetTitle("Invariant Higgs Mass [GeV]");
    HiggsMass->GetYaxis()->SetTitle("N/bin");
    HiggsMass->GetYaxis()->SetRangeUser(0.0, 1200);
    HiggsMass->SetMarkerSize(1.);
    HiggsMass->SetMarkerStyle(22.);
    HiggsMass->SetMarkerColor(kBlack);
    HiggsMass->SetLineColor(kBlack);
    HiggsMass->Draw("E0");
    //HiggsMass->Scale(1/MuonSize24->GetEntries());
    //HiggsMass->Scale(0.0283);
    HiggsMass->Fit("gaus");
    gPad->RedrawAxis();
    cout << "Efficiency of ggH125mumu: " << HiggsMass->GetEntries()/MuonSize24->GetEntries() << endl;
    Double_t HiggsSignal = 0;
    for (Int_t i=116;i<=134;i++) {
      HiggsSignal = HiggsSignal + HiggsMass->GetBinContent(i); }
    cout << "Higgs signal MC: " << HiggsSignal << endl;
    c59->Print("HiggsMass.pdf");

    // Draw truthZ
    auto c60 = new TCanvas("c60", "", 200, 10, 800, 600);
    truthZ->GetXaxis()->SetTitle("true Z mass [GeV]");
    truthZ->GetYaxis()->SetTitle("N/bin");
    truthZ->SetMarkerSize(1.);
    truthZ->SetMarkerStyle(22.);
    truthZ->SetMarkerColor(kBlack);
    truthZ->SetLineColor(kBlack);
    truthZ->Draw("E0");
    //truthZ->Fit("gaus");
    gPad->RedrawAxis();
    TF1 *fitB2 = new TF1("fitB2",bwconv,70,110,4);
    fitB2->SetLineColor(kBlue);
    fitB2->SetParameters(91,4,10000,2);
    fitB2->SetParNames("Mass", "fwhm", "Area", "Sigma");
    //fitB2->SetParLimits(4, 91, 92);
    //fitB2->SetParLimits(5, 0, 7);
    truthZ->Fit("fitB2","+","",70,110);
    c60->Print("truthZ.pdf");

    // Draw children
    auto c61 = new TCanvas("c61", "", 200, 10, 800, 600);
    children->GetXaxis()->SetTitle("true #mu Pt [GeV]");
    children->GetYaxis()->SetTitle("N/bin");
    children->SetMarkerSize(1.);
    children->SetMarkerStyle(22.);
    children->SetMarkerColor(kBlack);
    children->SetLineColor(kBlack);
    children->Draw("E0");
    //mChildren->Draw("same");
    gPad->RedrawAxis();
    c61->Print("children.pdf");

    // Draw mChildren
    auto c62 = new TCanvas("c62", "", 200, 10, 800, 600);
    mChildren->GetXaxis()->SetTitle("#mu truth match Pt [GeV]");
    mChildren->GetYaxis()->SetTitle("N/bin");
    mChildren->SetMarkerSize(1.);
    mChildren->SetMarkerStyle(22.);
    mChildren->SetMarkerColor(kBlack);
    mChildren->SetLineColor(kBlack);
    mChildren->Draw("E0");
    gPad->RedrawAxis();
    c62->Print("mChildren.pdf");

    // Draw truthEta
    auto c63 = new TCanvas("c63", "", 200, 10, 800, 600);
    truthEta->GetXaxis()->SetTitle("true #eta");
    truthEta->GetYaxis()->SetTitle("N/bin");
    truthEta->GetYaxis()->SetRangeUser(0.0, 800);
    truthEta->SetMarkerSize(1.);
    truthEta->SetMarkerStyle(22.);
    truthEta->SetMarkerColor(kBlack);
    truthEta->SetLineColor(kBlack);
    truthEta->Draw("E0");
    gPad->RedrawAxis();
    c63->Print("truthEta.pdf");

    // Draw truthPhi
    auto c64 = new TCanvas("c64", "", 200, 10, 800, 600);
    truthPhi->GetXaxis()->SetTitle("true #phi");
    truthPhi->GetYaxis()->SetTitle("N/bin");
    truthPhi->GetYaxis()->SetRangeUser(0.0, 1300);
    truthPhi->SetMarkerSize(1.);
    truthPhi->SetMarkerStyle(22.);
    truthPhi->SetMarkerColor(kBlack);
    truthPhi->SetLineColor(kBlack);
    truthPhi->Draw("E0");
    gPad->RedrawAxis();
    c64->Print("truthPhi.pdf");  

    // Draw matchEta
    auto c65 = new TCanvas("c65", "", 200, 10, 800, 600);
    matchEta->GetXaxis()->SetTitle("truth matched #eta");
    matchEta->GetYaxis()->SetTitle("N/bin");
    matchEta->GetYaxis()->SetRangeUser(0.0, 1400);
    matchEta->SetMarkerSize(1.);
    matchEta->SetMarkerStyle(22.);
    matchEta->SetMarkerColor(kBlack);
    matchEta->SetLineColor(kBlack);
    matchEta->Draw("E0");
    truthEta->Draw("same");
    gPad->RedrawAxis();
    c65->Print("matchEta.pdf");

    // Draw matchPhi
    auto c66 = new TCanvas("c66", "", 200, 10, 800, 600);
    matchPhi->GetXaxis()->SetTitle("truth matched #phi");
    matchPhi->GetYaxis()->SetTitle("N/bin");
    matchPhi->GetYaxis()->SetRangeUser(0.0, 1300);
    matchPhi->SetMarkerSize(1.);
    matchPhi->SetMarkerStyle(22.);
    matchPhi->SetMarkerColor(kBlack);
    matchPhi->SetLineColor(kBlack);
    matchPhi->Draw("E0");
    //truthPhi->Draw("same");
    gPad->RedrawAxis();
    c66->Print("matchPhi.pdf");

    //draw eta ratio
    auto c80 = new TCanvas("c80", "", 200, 10, 800, 600);
    TH1D *Ratio = new TH1D ("ratioEta","",40,-2.7,2.7);
    Ratio = (TH1D*) matchEta->Clone();
    Ratio->GetXaxis()->SetTitle("match vs truth #eta");
    Ratio->GetXaxis()->SetRangeUser(-2.7, 2.7);
    Ratio->GetYaxis()->SetTitle("efficiency");
    Ratio->GetYaxis()->SetRangeUser(0.0, 1.0);
    Ratio->Sumw2();
    Ratio->Divide(truthEta.GetPtr());
    Ratio->Draw("E0");
    gPad->RedrawAxis();
    c80->Print("ratioEta.pdf");

    //draw phi ratio
    auto c81 = new TCanvas("c81", "", 200, 10, 800, 600);
    TH1D *Ratio2 = new TH1D ("ratioPhi","",40,-3.14,3.14);
    Ratio2 = (TH1D*) mphi_cut->Clone();
    Ratio2->GetXaxis()->SetTitle("match vs truth #phi");
    Ratio2->GetXaxis()->SetRangeUser(-3.14, 3.14);
    Ratio2->GetYaxis()->SetTitle("efficiency");
    Ratio2->GetYaxis()->SetRangeUser(0.0, 1.0);
    Ratio2->Sumw2();
    Ratio2->Divide(tphi_cut.GetPtr());
    Ratio2->Draw("E0");
    gPad->RedrawAxis();
    c81->Print("ratioPhi.pdf");

    //draw pt ratio
    //auto c20 = new TCanvas("c20", "", 200, 10, 800, 600);
    //TH1D *Ratio3 = new TH1D ("ratioPt","",40,0.0,150);
    //Ratio3 = (TH1D*) mpt_cut->Clone();
    //Ratio3->GetXaxis()->SetTitle("match vs truth pt");
    //Ratio3->GetXaxis()->SetRangeUser(0.0, 100);
    //Ratio3->GetYaxis()->SetTitle("efficiency");
    //Ratio3->GetYaxis()->SetRangeUser(0.0, 1.0);
    //Ratio3->Sumw2();
    //Ratio3->Divide(tpt_cut.GetPtr());
    //Ratio3->Draw("E0");
    //gPad->RedrawAxis();
    //c20->Print("ratioPt.pdf");

    // Draw MuonPt with cuts
    //auto c21 = new TCanvas("c21", "", 200, 10, 800, 600);
    //pt_cut->GetXaxis()->SetTitle("pt");
    //pt_cut->GetYaxis()->SetTitle("N/bin");
    //pt_cut->GetYaxis()->SetRangeUser(0.0, 2500);
    //pt_cut->SetMarkerSize(1.);
    //pt_cut->SetMarkerStyle(22.);
    //pt_cut->SetMarkerColor(kBlack);
    //pt_cut->SetLineColor(kBlack);
    //pt_cut->Draw("E0");
    //gPad->RedrawAxis();
    //c21->Print("pt_cut.pdf");

    // Draw MuonEta with cuts
    //auto c22 = new TCanvas("c22", "", 200, 10, 800, 600);
    //eta_cut->GetXaxis()->SetTitle("#eta");
    //eta_cut->GetYaxis()->SetTitle("N/bin");
    //eta_cut->GetYaxis()->SetRangeUser(0.0, 1700);
    //eta_cut->SetMarkerSize(1.);
    //eta_cut->SetMarkerStyle(22.);
    //eta_cut->SetMarkerColor(kBlack);
    //eta_cut->SetLineColor(kBlack);
    //eta_cut->Draw("E0");
    //gPad->RedrawAxis();
    //c22->Print("eta_cut.pdf");

    // Draw matchPhi with cuts
    //auto c23 = new TCanvas("c23", "", 200, 10, 800, 600);
    //mphi_cut->GetXaxis()->SetTitle("truth matched #phi");
    //mphi_cut->GetYaxis()->SetTitle("N/bin");
    //mphi_cut->GetYaxis()->SetRangeUser(0.0, 1300);
    //mphi_cut->SetMarkerSize(1.);
    //mphi_cut->SetMarkerStyle(22.);
    //mphi_cut->SetMarkerColor(kBlack);
    //mphi_cut->SetLineColor(kBlack);
    //mphi_cut->Draw("E0");
    //gPad->RedrawAxis();
    //c23->Print("mphi_cut.pdf");

    // Draw truthPhi with cuts
    //auto c24 = new TCanvas("c24", "", 200, 10, 800, 600);
    //tphi_cut->GetXaxis()->SetTitle("true #phi");
    //tphi_cut->GetYaxis()->SetTitle("N/bin");
    //tphi_cut->GetYaxis()->SetRangeUser(0.0, 1300);
    //tphi_cut->SetMarkerSize(1.);
    //tphi_cut->SetMarkerStyle(22.);
    //tphi_cut->SetMarkerColor(kBlack);
    //tphi_cut->SetLineColor(kBlack);
    //tphi_cut->Draw("E0");
    //gPad->RedrawAxis();
    //c24->Print("tphi_cut.pdf");

    // Draw mChildren with cuts
    //auto c25 = new TCanvas("c25", "", 200, 10, 800, 600);
    //mpt_cut->GetXaxis()->SetTitle("#mu truth match Pt [GeV]");
    //mpt_cut->GetYaxis()->SetTitle("N/bin");
    //mpt_cut->SetMarkerSize(1.);
    //mpt_cut->SetMarkerStyle(22.);
    //mpt_cut->SetMarkerColor(kBlack);
    //mpt_cut->SetLineColor(kBlack);
    //mpt_cut->Draw("E0");
    //gPad->RedrawAxis();
    //c25->Print("mpt_cut.pdf");

    // Draw children with cuts
    //auto c26 = new TCanvas("c26", "", 200, 10, 800, 600);
    //tpt_cut->GetXaxis()->SetTitle("true #mu Pt [GeV]");
    //tpt_cut->GetYaxis()->SetTitle("N/bin");
    //tpt_cut->SetMarkerSize(1.);
    //tpt_cut->SetMarkerStyle(22.);
    //tpt_cut->SetMarkerColor(kBlack);
    //tpt_cut->SetLineColor(kBlack);
    //tpt_cut->Draw("E0");
    //gPad->RedrawAxis();
    //c26->Print("tpt_cut.pdf");

    // Draw pt resolution
    auto c67 = new TCanvas("c67", "", 200, 10, 800, 600);
    ptRes->GetXaxis()->SetTitle("1/truth - 1/match pt [GeV]");
    ptRes->GetYaxis()->SetTitle("");
    ptRes->GetYaxis()->SetRangeUser(0.0, 6500);
    ptRes->SetMarkerSize(1.);
    ptRes->SetMarkerStyle(22.);
    ptRes->SetMarkerColor(kBlack);
    ptRes->SetLineColor(kBlack);
    ptRes->Draw("E0");
    ptRes->Fit("gaus");
    gPad->RedrawAxis();
    c67->Print("ptRes.pdf");

    // Draw eta resolution
    auto c68 = new TCanvas("c68", "", 200, 10, 800, 600);
    etaRes->GetXaxis()->SetTitle("truth - match #eta");
    etaRes->GetYaxis()->SetTitle("");
    //etaRes->GetYaxis()->SetRangeUser(0.0, 1700);
    etaRes->SetMarkerSize(1.);
    etaRes->SetMarkerStyle(22.);
    etaRes->SetMarkerColor(kBlack);
    etaRes->SetLineColor(kBlack);
    etaRes->Draw("E0");
    gPad->RedrawAxis();
    c68->Print("etaRes.pdf");

    // Draw MuonScatter
    //auto c29 = new TCanvas("c29", "", 200, 10, 800, 600);
    //c29->SetRightMargin(0.18);
    //scatterTruth->GetXaxis()->SetTitle("truth Pt [GeV]");
    //scatterTruth->GetYaxis()->SetTitle("truth #eta");
    //scatterTruth->GetXaxis()->SetRangeUser(0.0, 150);
    //scatterTruth->GetYaxis()->SetRangeUser(-2.7, 2.7);
    //scatterTruth->SetMarkerSize(1.);
    //scatterTruth->SetMarkerStyle(22.);
    //scatterTruth->SetMarkerColor(kBlack);
    //scatterTruth->SetLineColor(kBlack);
    //scatterTruth->Draw("colz");
    //gPad->RedrawAxis();
    //c29->Print("scatterTruth.pdf");

    // Draw MuonScatter
    //auto c30 = new TCanvas("c30", "", 200, 10, 800, 600);
    //c30->SetRightMargin(0.18);
    //scatterMatch->GetXaxis()->SetTitle("match Pt [GeV]");
    //scatterMatch->GetYaxis()->SetTitle("match #eta");
    //scatterMatch->GetXaxis()->SetRangeUser(0.0, 150);
    //scatterMatch->GetYaxis()->SetRangeUser(-2.7, 2.7);
    //scatterMatch->SetMarkerSize(1.);
    //scatterMatch->SetMarkerStyle(22.);
    //scatterMatch->SetMarkerColor(kBlack);
    //scatterMatch->SetLineColor(kBlack);
    //scatterMatch->Draw("colz");
    //gPad->RedrawAxis();
    //c30->Print("scatterMatch.pdf");

    // Draw pt resolution 0-10 GeV
    auto c69 = new TCanvas("c69", "", 200, 10, 800, 600);
    auto ptBin1 = d.Define("ptBin1","ptRes[mChildren > 0 && mChildren < 10]");
    auto ptRes1 = ptBin1.Histo1D({"ptRes1","",40,-0.02,0.02},{"ptBin1"});
    ptRes1->GetXaxis()->SetTitle("1/truth - 1/match pt [0-10 GeV]");
    ptRes1->GetYaxis()->SetTitle("");
    ptRes1->SetMarkerSize(1.);
    ptRes1->SetMarkerStyle(22.);
    ptRes1->SetMarkerColor(kBlack);
    ptRes1->SetLineColor(kBlack);
    ptRes1->Draw("E0");
    ptRes1->Fit("gaus");
    Double_t width1 = ptRes1->GetStdDev();
    Double_t error1 = ptRes1->GetStdDevError();
    gPad->RedrawAxis();
    c69->Print("ptRes1.pdf");

    // Draw pt resolution 10-20 GeV
    auto c70 = new TCanvas("c70", "", 200, 10, 800, 600);
    auto ptBin2 = d.Define("ptBin2","ptRes[mChildren > 10 && mChildren < 20]");
    auto ptRes2 = ptBin2.Histo1D({"ptRes2","",40,-0.02,0.02},{"ptBin2"});
    ptRes2->GetXaxis()->SetTitle("1/truth - 1/match pt [10-20 GeV]");
    ptRes2->GetYaxis()->SetTitle("");
    ptRes2->SetMarkerSize(1.);
    ptRes2->SetMarkerStyle(22.);
    ptRes2->SetMarkerColor(kBlack);
    ptRes2->SetLineColor(kBlack);
    ptRes2->Draw("E0");
    ptRes2->Fit("gaus");
    Double_t width2 = ptRes2->GetStdDev();
    Double_t error2 = ptRes2->GetStdDevError();
    gPad->RedrawAxis();
    c70->Print("ptRes2.pdf");

    // Draw pt resolution 20-30 GeV
    auto c71 = new TCanvas("c71", "", 200, 10, 800, 600);
    auto ptBin3 = d.Define("ptBin3","ptRes[mChildren > 20 && mChildren < 30]");
    auto ptRes3 = ptBin3.Histo1D({"ptRes3","",40,-0.02,0.02},{"ptBin3"});
    ptRes3->GetXaxis()->SetTitle("1/truth - 1/match pt [20-30 GeV]");
    ptRes3->GetYaxis()->SetTitle("");
    ptRes3->SetMarkerSize(1.);
    ptRes3->SetMarkerStyle(22.);
    ptRes3->SetMarkerColor(kBlack);
    ptRes3->SetLineColor(kBlack);
    ptRes3->Draw("E0");
    ptRes3->Fit("gaus");
    Double_t width3 = ptRes3->GetStdDev();
    Double_t error3 = ptRes3->GetStdDevError();
    gPad->RedrawAxis();
    c71->Print("ptRes3.pdf");

    // Draw pt resolution 30-40 GeV
    auto c72 = new TCanvas("c72", "", 200, 10, 800, 600);
    auto ptBin4 = d.Define("ptBin4","ptRes[mChildren > 30 && mChildren < 40]");
    auto ptRes4 = ptBin4.Histo1D({"ptRes4","",40,-0.02,0.02},{"ptBin4"});
    ptRes4->GetXaxis()->SetTitle("1/truth - 1/match pt [30-40 GeV]");
    ptRes4->GetYaxis()->SetTitle("");
    ptRes4->SetMarkerSize(1.);
    ptRes4->SetMarkerStyle(22.);
    ptRes4->SetMarkerColor(kBlack);
    ptRes4->SetLineColor(kBlack);
    ptRes4->Draw("E0");
    ptRes4->Fit("gaus");
    Double_t width4 = ptRes4->GetStdDev();
    Double_t error4 = ptRes4->GetStdDevError();
    gPad->RedrawAxis();
    c72->Print("ptRes4.pdf");

    // Draw pt resolution 40-50 GeV
    auto c73 = new TCanvas("c73", "", 200, 10, 800, 600);
    auto ptBin5 = d.Define("ptBin5","ptRes[mChildren > 40 && mChildren < 50]");
    auto ptRes5 = ptBin5.Histo1D({"ptRes5","",40,-0.02,0.02},{"ptBin5"});
    ptRes5->GetXaxis()->SetTitle("1/truth - 1/match pt [40-50 GeV]");
    ptRes5->GetYaxis()->SetTitle("");
    ptRes5->SetMarkerSize(1.);
    ptRes5->SetMarkerStyle(22.);
    ptRes5->SetMarkerColor(kBlack);
    ptRes5->SetLineColor(kBlack);
    ptRes5->Draw("E0");
    ptRes5->Fit("gaus");
    Double_t width5 = ptRes5->GetStdDev();
    Double_t error5 = ptRes5->GetStdDevError();
    gPad->RedrawAxis();
    c73->Print("ptRes5.pdf");

    // Draw pt resolution 50-60 GeV
    auto c74 = new TCanvas("c74", "", 200, 10, 800, 600);
    auto ptBin6 = d.Define("ptBin6","ptRes[mChildren > 50 && mChildren < 60]");
    auto ptRes6 = ptBin6.Histo1D({"ptRes6","",40,-0.02,0.02},{"ptBin6"});
    ptRes6->GetXaxis()->SetTitle("1/truth - 1/match pt [50-60 GeV]");
    ptRes6->GetYaxis()->SetTitle("");
    ptRes6->SetMarkerSize(1.);
    ptRes6->SetMarkerStyle(22.);
    ptRes6->SetMarkerColor(kBlack);
    ptRes6->SetLineColor(kBlack);
    ptRes6->Draw("E0");
    ptRes6->Fit("gaus");
    Double_t width6 = ptRes6->GetStdDev();
    Double_t error6 = ptRes6->GetStdDevError();
    gPad->RedrawAxis();
    c74->Print("ptRes6.pdf");

    // Draw pt resolution 60-70 GeV
    auto c75 = new TCanvas("c75", "", 200, 10, 800, 600);
    auto ptBin7 = d.Define("ptBin7","ptRes[mChildren > 60 && mChildren < 70]");
    auto ptRes7 = ptBin7.Histo1D({"ptRes7","",40,-0.02,0.02},{"ptBin7"});
    ptRes7->GetXaxis()->SetTitle("1/truth - 1/match pt [60-70 GeV]");
    ptRes7->GetYaxis()->SetTitle("");
    ptRes7->SetMarkerSize(1.);
    ptRes7->SetMarkerStyle(22.);
    ptRes7->SetMarkerColor(kBlack);
    ptRes7->SetLineColor(kBlack);
    ptRes7->Draw("E0");
    ptRes7->Fit("gaus");
    Double_t width7 = ptRes7->GetStdDev();
    Double_t error7 = ptRes7->GetStdDevError();
    gPad->RedrawAxis();
    c75->Print("ptRes7.pdf");

    // Draw pt resolution 70-80 GeV
    auto c76 = new TCanvas("c76", "", 200, 10, 800, 600);
    auto ptBin8 = d.Define("ptBin8","ptRes[mChildren > 70 && mChildren < 80]");
    auto ptRes8 = ptBin8.Histo1D({"ptRes8","",40,-0.02,0.02},{"ptBin8"});
    ptRes8->GetXaxis()->SetTitle("1/truth - 1/match pt [70-80 GeV]");
    ptRes8->GetYaxis()->SetTitle("");
    ptRes8->SetMarkerSize(1.);
    ptRes8->SetMarkerStyle(22.);
    ptRes8->SetMarkerColor(kBlack);
    ptRes8->SetLineColor(kBlack);
    ptRes8->Draw("E0");
    ptRes8->Fit("gaus");
    Double_t width8 = ptRes8->GetStdDev();
    Double_t error8 = ptRes8->GetStdDevError();
    gPad->RedrawAxis();
    c76->Print("ptRes8.pdf");

    // Draw pt resolution 80-90 GeV
    auto c77 = new TCanvas("c77", "", 200, 10, 800, 600);
    auto ptBin9 = d.Define("ptBin9","ptRes[mChildren > 80 && mChildren < 90]");
    auto ptRes9 = ptBin9.Histo1D({"ptRes9","",40,-0.02,0.02},{"ptBin9"});
    ptRes9->GetXaxis()->SetTitle("1/truth - 1/match pt [80-90 GeV]");
    ptRes9->GetYaxis()->SetTitle("");
    ptRes9->SetMarkerSize(1.);
    ptRes9->SetMarkerStyle(22.);
    ptRes9->SetMarkerColor(kBlack);
    ptRes9->SetLineColor(kBlack);
    ptRes9->Draw("E0");
    ptRes9->Fit("gaus");
    Double_t width9 = ptRes9->GetStdDev();
    Double_t error9 = ptRes9->GetStdDevError();
    gPad->RedrawAxis();
    c77->Print("ptRes9.pdf");

    // Draw pt resolution 90-100 GeV
    auto c78 = new TCanvas("c78", "", 200, 10, 800, 600);
    auto ptBin10 = d.Define("ptBin10","ptRes[mChildren > 90 && mChildren < 100]");
    auto ptRes10 = ptBin10.Histo1D({"ptRes10","",40,-0.02,0.02},{"ptBin10"});
    ptRes10->GetXaxis()->SetTitle("1/truth - 1/match pt [90-100 GeV]");
    ptRes10->GetYaxis()->SetTitle("");
    ptRes10->SetMarkerSize(1.);
    ptRes10->SetMarkerStyle(22.);
    ptRes10->SetMarkerColor(kBlack);
    ptRes10->SetLineColor(kBlack);
    ptRes10->Draw("E0");
    ptRes10->Fit("gaus");
    Double_t width10 = ptRes10->GetStdDev();
    Double_t error10 = ptRes10->GetStdDevError();
    gPad->RedrawAxis();
    c78->Print("ptRes10.pdf");

    // plot std devs
    auto c79 = new TCanvas("c79", "", 200, 10, 800, 600);
    TH1D *TRes = new TH1D ("ResPt","",10,0.0,100);
    TRes->SetBinContent(1, 5, width1*5);
    TRes->SetBinError(1, error1*5);
    TRes->SetBinContent(2, 15, width2*15);
    TRes->SetBinError(2, error2*15);
    TRes->SetBinContent(3, 25, width3*25);
    TRes->SetBinError(3, error3*25);
    TRes->SetBinContent(4, 35, width4*35);
    TRes->SetBinError(4, error4*35);
    TRes->SetBinContent(5, 45, width5*45);
    TRes->SetBinError(5, error5*45);
    TRes->SetBinContent(6, 55, width6*55);
    TRes->SetBinError(6, error6*55);
    TRes->SetBinContent(7, 65, width7*65);
    TRes->SetBinError(7, error7*65);
    TRes->SetBinContent(8, 75, width8*75);
    TRes->SetBinError(8, error8*75);
    TRes->SetBinContent(9, 85, width9*85);
    TRes->SetBinError(9, error9*85);
    TRes->SetBinContent(10, 95, width10*95);
    TRes->SetBinError(10, error10*95);
    TRes->GetXaxis()->SetTitle("Pt [GeV]");
    TRes->GetYaxis()->SetTitle("#sigma(pt)/pt");
    TRes->GetYaxis()->SetRangeUser(0,0.19);
    TRes->SetMarkerSize(1.);
    TRes->SetMarkerStyle(22.);
    TRes->SetMarkerColor(kBlack);
    TRes->SetLineColor(kBlack);
    TRes->Draw("E0");
    gPad->RedrawAxis();
    c79->Print("StdDev.pdf");
}
