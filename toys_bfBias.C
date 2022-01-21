#include "/scratchfs/atlas/huirun/atlaswork/ATLAS_style/atlasrootstyle/AtlasUtils.h"
#include "/scratchfs/atlas/huirun/atlaswork/ATLAS_style/atlasrootstyle/AtlasStyle.h"
#include "/scratchfs/atlas/huirun/atlaswork/ATLAS_style/atlasrootstyle/AtlasLabels.h"

#ifdef __CLING__
#include "/scratchfs/atlas/huirun/atlaswork/ATLAS_style/atlasrootstyle/AtlasLabels.C"
#include "/scratchfs/atlas/huirun/atlaswork/ATLAS_style/atlasrootstyle/AtlasUtils.C"
#endif

#include "utils.h"
#include "/afs/ihep.ac.cn/users/g/guofy/HggTwoSidedCBPdf.cxx"
#include "/afs/ihep.ac.cn/users/g/guofy/HggTwoSidedCBPdf.h"

#include <string>

#include "Riostream.h"

#include "RooDSCBShape.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"

#include "RooProdPdf.h"

using namespace RooFit;
using namespace std;

ClassImp(RooDSCBShape)

RooDSCBShape::RooDSCBShape(const char *name, const char *title,
  RooAbsReal& _x,
  RooAbsReal& _mu,
  RooAbsReal& _sig,
  RooAbsReal& _a1,
  RooAbsReal& _n1,
  RooAbsReal& _a2,
  RooAbsReal& _n2) :
  RooAbsPdf(name,title),
  x("x","x",this,_x),
  mu("mu","mu",this,_mu),
  sig("sig","sig",this,_sig),
  a1("a1","a1",this,_a1), //a1 must be > 0
  n1("n1","n1",this,_n1),
  a2("a2","a2",this,_a2), //a2 must be > 0
  n2("n2","n2",this,_n2)
  {
  }


  RooDSCBShape::RooDSCBShape(const RooDSCBShape& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  mu("mu",this,other.mu),
  sig("sig",this,other.sig),
  a1("a1",this,other.a1),
  n1("n1",this,other.n1),
  a2("a2",this,other.a2),
  n2("n2",this,other.n2)
  {
  }



  Double_t RooDSCBShape::evaluate() const
  {
    double u   = (x-mu)/sig;
    double A1  = TMath::Power(n1/TMath::Abs(a1),n1)*TMath::Exp(-a1*a1/2);
    double A2  = TMath::Power(n2/TMath::Abs(a2),n2)*TMath::Exp(-a2*a2/2);
    double B1  = n1/a1 - a1;
    double B2  = n2/a2 - a2;

    double result(1);
    if      (u<-TMath::Abs(a1)) result *= A1*TMath::Power(B1-u,-n1);
    else if (u<TMath::Abs(a2))  result *= TMath::Exp(-u*u/2);
    else                        result *= A2*TMath::Power(B2+u,-n2);
    return result;
  }


  Int_t RooDSCBShape::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char*) const
  {
    if (matchArgs(allVars,analVars,x)) return 1 ;
    return 0 ;
  }

  Double_t RooDSCBShape::analyticalIntegral(Int_t code, const char* r) const
  {
    double umin = (x.min(r) - mu) / sig;
    double umax = (x.max(r) - mu) / sig;
    R__ASSERT(code==1);

    double integral = 0.;

    integral += IntPwLw(TMath::Max(-umax, TMath::Abs(a1)), TMath::Max(-umin, TMath::Abs(a1)), a1, n1);
    integral += IntGaus(TMath::Max(umin, -TMath::Abs(a1)), TMath::Min(umax, TMath::Abs(a2)));
    integral += IntPwLw(TMath::Max(umin, TMath::Abs(a2)), TMath::Max(umax, TMath::Abs(a2)), a2, n2);

    return sig * integral;
  }


  double RooDSCBShape::IntGaus(double x0, double x1) const
  {
    static const double rootPiBy2 = TMath::Sqrt(TMath::PiOver2());

    if (x0 >= x1) return 0; // needed in case umin > a2
    if (x0*x1<0) // they are at different side of zero
    {
      return rootPiBy2 * ( TMath::Erf(TMath::Abs(x1) / TMath::Sqrt2()) + TMath::Erf(TMath::Abs(x0) / TMath::Sqrt2()) );
    }
    else //They are at the same side of zero
    {
      return rootPiBy2 * TMath::Abs( TMath::Erf(TMath::Abs(x1) / TMath::Sqrt2()) - TMath::Erf(TMath::Abs(x0) / TMath::Sqrt2()) );
    }

  }

  double RooDSCBShape::IntPwLw(double x0, double x1, double alpha, double n) const
  {

    if (x0 == x1) return 0; // already implicit below but so it's clear

    bool useLog = false;
    if(fabs(n - 1.0) < 1.0e-05)
    useLog = true;

    double A  = TMath::Power(n/TMath::Abs(alpha),n)*TMath::Exp(-alpha*alpha/2);
    double B  = n/TMath::Abs(alpha) - TMath::Abs(alpha);

    double result = 0.;
    if(useLog)
    {
      result = A * ( TMath::Log(B + x1) - TMath::Log(B + x0));
    }
    else
    {
      result = A / (1. - n) * ( TMath::Power(B + x1, 1. - n) - TMath::Power(B + x0, 1. - n) );
    }
    return result;
  }

TH1F *scaleHistXaxis(TH1F* hist, float scale){
  int nbins = hist->GetNbinsX();
  float min = hist->GetXaxis()->GetXmin();
  float max = hist->GetXaxis()->GetXmax();
  cout<<"nbins, min, max: "<<nbins<<", "<<min<<", "<<max<<endl;
  TH1F *hout = new TH1F("hout", "", nbins, min*scale, max*scale);
  for (int i = 1; i <= nbins; i++){
    float bcon = hist->GetBinContent(i);
    float berr = hist->GetBinError(i);
    hout->SetBinContent(i, bcon);
    hout->SetBinError(i, berr);
  }
  return hout;
}

void toys_bfBias()
{
  SetAtlasStyle();

  char *cf_cats = (char*)"cats.cfg";
  char *cf_bins = (char*)"binnings.cfg";
  map<TString, string> catCuts;
  map<TString, vector<float>> catBins;
  getCatCuts(cf_cats, catCuts); for(auto c : catCuts) cout<<c.first<<c.second<<endl;
  getCatBinning(cf_bins, catBins); for(auto c : catBins) { cout<<"dafad"<<c.first<<endl; for(auto b : c.second)  cout<<b<<", "; cout<<endl<<endl; }

  map<TString, TString> bins_name;
  bins_name["b1"] = "bin0";
  bins_name["b2"] = "bin1";
  bins_name["b3"] = "bin2";
  bins_name["b4"] = "bin3";
  bins_name["b5"] = "bin4";
  bins_name["b6"] = "bin5";
  //bins_name["b7"] = "bin6";
  //bins_name["b8"] = "bin7";
  //bins_name["b9"] = "bin8";
  //bins_name["b10"] = "bin9";

  enum Para{B, C, POWA, YIELD};
  map<TString, vector<float>> m_para;

  double myy_out, oo_out, weight_out;

  TString funcname = "ExpPoly2";

  TFile *f_out = new TFile("tree_bfBias_Asimov_"+funcname+".root", "recreate");
  TFile *f_out_toy = new TFile("tree_bfBias_toys.root", "recreate");
  //TFile *f_out = new TFile("tree_asimov_minus_0010.root", "recreate");

for(auto cat : catCuts){
  TString catName = cat.first; cout<<"=== "<<catName<<" ==="<<endl;
  string catCut = cat.second;

TString bdtCat;
Ssiz_t from = 0;
catName.Tokenize(bdtCat, from, "[_]");


    TCanvas *canv = new TCanvas("c", "canvas", 800, 600);

    cout<<"=============="<<catName<<"=============="<<endl;

    TTree *t_out = new TTree(Form("CollectionTree_%s", catName.Data()), "");
 
    //t_out->Branch("m_yy_SR", &myy_out, "m_yy_SR/D");
    t_out->Branch("m_yy", &myy_out, "m_yy/D");
    //t_out->Branch("oo1", &oo_out, "oo1/D");
    t_out->Branch("weight", &weight_out, "weight/D");
  
//    TFile *ft_SB = new TFile(Form("tree_dataSB_OO_%s.root", catName.Data()), "read");
//    TTree *t_SB = (TTree*) ft_SB->Get("CollectionTree");
//    TFile *fh_SB = new TFile(Form("hists_dataSB_%s.root", catName.Data()), "read");
    TFile *fh_SB = new TFile("../bkgTemplate/template_"+bdtCat+".root", "read");
    TH1F *h_SB = (TH1F*) fh_SB->Get("template_uncer_"+catName);
    float n_SB = h_SB->Integral();

    TH1F *h_SB_MeV = scaleHistXaxis(h_SB, 1000);
 
    TFile *f_ggh = new TFile(Form("MC16Xs/fullrun2/hists_ggF_%s.root", catName.Data()), "read");
    TTree *t_ggh = (TTree*) f_ggh->Get("CollectionTree");
    TH1F *h_ggh = (TH1F*) f_ggh->Get("m_yy_ggF");
    float n_ggh = h_ggh->Integral();
  
    TFile *f_vbf = new TFile(Form("MC16Xs/fullrun2/hists_VBF_%s.root", catName.Data()), "read");
    //TFile *f_vbf = new TFile("hists_vbf_minus_0010_"+bin->first+".root", "read");
    TTree *t_vbf = (TTree*) f_vbf->Get("CollectionTree");
    TH1F *h_vbf = (TH1F*) f_vbf->Get("m_yy_VBF");
    float n_vbf = h_vbf->Integral();
  
    TH1F h_sig = *h_ggh+*h_vbf;
    float n_sig = n_ggh+n_vbf;
  
    float SR_up = 120000;
    float SR_lo = 130000;
  
    RooRealVar m_yy("m_yy", "m_yy", 105000, 160000);
    RooRealVar oo1("oo1", "oo1", -20, 20);
    RooRealVar weight("weight", "weight", -10, 100);
    RooPlot* oo_frame = oo1.frame(Title("oo projection of DSCB(oo)*polyexp(m_yy)")) ;
  
    RooDataHist dh_SB("dh_SB", "", m_yy, Import(*h_SB_MeV));// !!! order
//    RooDataHist dh_sig("dh_sig", "", m_yy, Import(h_sig));
 
//    RooDataSet ds_SB("ds_SB","", RooArgSet(m_yy, weight), WeightVar(weight), Import(*t_SB));
    RooDataSet ds_ggh("ds_ggh","", RooArgSet(m_yy, weight), WeightVar(weight), Import(*t_ggh));
    RooDataSet ds_vbf("ds_vbf","", RooArgSet(m_yy, weight), WeightVar(weight), Import(*t_vbf));
    ds_ggh.append(ds_vbf);
    RooDataSet ds_sig = ds_ggh;
 
    RooRealVar b("b", "", -12, 12);
    RooRealVar c("c", "", -15, 15);
    RooRealVar powa("powa", "", -100, 100);
    //RooGenericPdf polyexp("polyexp", "fit sideband m_yy", "pow(m_yy/1000/125, powa)", RooArgSet(m_yy, powa));// Pow
    //RooGenericPdf polyexp("polyexp", "fit sideband m_yy", "exp(b*(m_yy/1000-125)/125)", RooArgSet(b, m_yy));// Exp
    RooGenericPdf polyexp("polyexp", "fit sideband m_yy", "exp(b*(m_yy/1000-125)/125+c*((m_yy/1000-125)/125)*((m_yy/1000-125)/125))", RooArgSet(b, c, m_yy));// ExpPoly2
  
    m_yy.setRange("SR", SR_up, SR_lo);
    m_yy.setRange("SB1", 105000, SR_up);
    m_yy.setRange("SB2", SR_lo, 160000);
  
    //polyexp.fitTo(dh_SB, Range("SB1,SB2"), Save());
    polyexp.fitTo(dh_SB, Save());
//    polyexp.fitTo(ds_SB, Range("SB1,SB2"), Save());

    m_para[catName].push_back(b.getVal());
    //m_para[bin->first].push_back(0.1);// Exp
    m_para[catName].push_back(c.getVal());// ExpPoly2
    m_para[catName].push_back(powa.getVal());// Pow

    RooPlot* myy_frame = m_yy.frame(Title("fit to m_yy distribution"));
    //ds_SB.plotOn(myy_frame);
    //polyexp.plotOn(myy_frame);
//    myy_frame->Draw();
  
    RooAbsReal *int_SB1 = polyexp.createIntegral(m_yy, NormSet(m_yy), Range("SB1"));
    RooAbsReal *int_SB2 = polyexp.createIntegral(m_yy, NormSet(m_yy), Range("SB2"));
    RooAbsReal *int_SR = polyexp.createIntegral(m_yy, NormSet(m_yy), Range("SR"));
    float frac_SR_SB = int_SR->getVal()/(int_SB1->getVal()+int_SB2->getVal());
    float frac_SR = int_SR->getVal()/(int_SR->getVal()+int_SB1->getVal()+int_SB2->getVal());
    cout<<"int_SR: "<<int_SR->getVal()<<endl;
    cout<<"int_SB1: "<<int_SB1->getVal()<<endl;
    cout<<"int_SB2: "<<int_SB2->getVal()<<endl;
  
    float n_SR = n_SB*frac_SR;cout<<"n_SR: "<<n_SR<<endl;cout<<"n_SB: "<<n_SB<<endl;
    float n_bkg = n_SB;
 
    m_para[catName].push_back(n_bkg);
 
    RooRealVar mean("mean","mean",120000,130000);// MeV or GeV
    RooRealVar sigma("sigma","sigma",700,10000);
    RooRealVar n1("n1","",0,25);
    RooRealVar alpha1("alpha1","",0,3);
    RooRealVar n2("n2","",0,45);
    RooRealVar alpha2("alpha2","",0,3);
    HggTwoSidedCBPdf DSCB_myy("sig","signal component",m_yy,mean,sigma,alpha1,n1,alpha2,n2);
  
    //DSCB_myy.fitTo(dh_sig, Save());
    DSCB_myy.fitTo(ds_sig, Save());
    ds_sig.plotOn(myy_frame);
    DSCB_myy.plotOn(myy_frame);
    myy_frame->Draw();

    canv->SaveAs("plotFit/"+catName+".png");
    delete canv;

    RooRealVar N_sig("N_sig","", n_sig);
    RooRealVar N_bkg("N_bkg","", n_bkg);
    RooAddPdf model_asiv("model_asiv","", RooArgList(DSCB_myy, polyexp), RooArgList(N_sig, N_bkg));

//    ofstream ofdata("data_OO_"+bin->first+".txt",ios::out);
//    if(!ofdata){
//      ofdata.close();
//      cout<<"can't open file for writing"<<endl;
//    }
 
    TCanvas *c1 = new TCanvas("c1", "canvas", 800, 600);
    RooPlot* myyFr1 = m_yy.frame(Title("fit to m_yy distribution"));
    model_asiv.plotOn(myyFr1);
    myyFr1->Draw();
    float frMax = myyFr1->GetMaximum();
    cout<<"maximum of myy: "<<frMax<<endl;
    c1->SaveAs("plotFit/model_asiv_"+catName+".png");
 
    int n_points = 220;
    double interval = (160000.-105000.)/n_points;

    TH1F *hout_asi = new TH1F(Form("Asi_%s", catName.Data()), "", n_points, 105000, 160000);
  
    RooArgSet nset(m_yy);
  
    double intPdf = 0;
    for(int i=0;i<n_points;i++){
      //myy_out = 105+interval*i;//cout<<"m_yy:"<<myy_out<<endl;
      //m_yy.setVal(myy_out);
//      myy_out = 1000*(105+interval*i);//cout<<"m_yy:"<<myy_out<<endl;
      myy_out = 105000+interval*i;//cout<<"m_yy:"<<myy_out<<endl;
      //if(myy_out>130000||myy_out<120000) continue;
      m_yy.setVal(myy_out);
      double pdf_val = model_asiv.getVal(&nset);//cout<<"pdf: "<<pdf_val<<endl;
      weight_out = pdf_val*interval*(n_sig+n_bkg);//cout<<"weight: "<<weight_out<<endl;
      t_out->Fill();
      hout_asi->Fill(myy_out, weight_out);
      intPdf += pdf_val*interval;

//      for(int i = 0; i < (int)(weight_out+0.5); i++){
//        ofdata<<myy_out<<endl;
//      }
    }
    f_out->cd();
    t_out->Write();
    hout_asi->Write();
    delete t_out;

    cout<<"intPdf: "<<intPdf<<endl;
    cout<<"n_ggh, n_vbf: "<<n_ggh<<", "<<n_vbf<<endl;
    cout<<"n_SB, n_sig: "<<n_SB<<", "<<n_sig<<endl;
    cout<<"n_SR: "<<n_SR<<endl;

    int nToys = 50;
    int itoytodraw = 1;

    TRandom rdm;
    TRandom rdm1;
    TRandom rdm2;
    rdm.SetSeed(10*n_sig*n_bkg);
    rdm1.SetSeed(n_sig/n_bkg);
    rdm2.SetSeed(n_bkg/n_sig);
    for(int i = 0; i < nToys; i++){
      TH1F *hout_toy = new TH1F(Form("toy_%s_%i", catName.Data(), i), "", n_points, 105000, 160000);

      //int rNtot = rdm.Poisson(n_sig+n_bkg);
      //int rNtot = floor(n_sig+n_bkg) + (rdm.Uniform() > 0.5 ? 1 : 0);
      int rNtot = floor(n_sig+n_bkg) + rdm.Poisson((n_sig+n_bkg)-floor(n_sig+n_bkg));
      cout<<catName<<" "<<i<<" toy, Ntot RNtot: "<<n_sig+n_bkg<<", "<<rNtot<<endl;
      int count = 0;
      while(count < rNtot){
        float rX = rdm1.Uniform();
        float rY = rdm2.Uniform();
  
        float myy = 55000*rX+105000;
        float pdf = frMax*rY;
  
        m_yy.setVal(myy);
        float pdf_val = model_asiv.getVal(&nset);
  
        if(pdf_val < pdf) continue;
        hout_toy->Fill(myy);
        count++;
      }
      if (itoytodraw == i && funcname == "ExpPoly2") {
        RooDataHist dh_toy("dh_toy", "", m_yy, Import(*hout_toy));// !!! order
        TCanvas *c2 = new TCanvas("c1", "canvas", 800, 600);
        RooPlot* myyFr2 = m_yy.frame();
        myyFr2->GetXaxis()->SetTitle("m_{#gamma#gamma} [MeV]");
        myyFr2->SetTitle(Form("toy %i "+catName, i));
        dh_toy.plotOn(myyFr2, Binning(55, 105000, 160000));
        model_asiv.plotOn(myyFr2);
        myyFr2->Draw();
        //model_asiv.plotOn(myyFr2, Normalization(n_sig+n_bkg, RooAbsReal::NumEvent));
        //myyFr2->addTH1(hout_toy);
        //hout_toy->Draw("same e");
        cout<<n_sig+n_bkg<<", "<<hout_toy->Integral()<<endl;

        myText(0.22, 0.88, 1, Form(catName+", toy %i from "+funcname, i));

        c2->SaveAs(Form("plotFit/toy%i_"+catName+"_"+funcname+".png", i));
      }

      f_out_toy->cd();
      hout_toy->Write();
    }

    if (funcname == "ExpPoly2") {
      TCanvas *c3 = new TCanvas("c1", "canvas", 800, 600);
      RooPlot* myyFr3 = m_yy.frame();
      myyFr3->GetXaxis()->SetTitle("m_{#gamma#gamma} [MeV]");
      dh_SB.plotOn(myyFr3, Binning(55, 105000, 160000));
      polyexp.plotOn(myyFr3);
      myyFr3->Draw();

      myText(0.22, 0.88, 1, (catName+", background fit with "+funcname).Data());

      c3->SaveAs("plotFit/bkgfit_"+catName+"_"+funcname+".png");

      TCanvas *c4 = new TCanvas("c1", "canvas", 800, 600);
      RooPlot* myyFr4 = m_yy.frame();
      myyFr4->GetXaxis()->SetTitle("m_{#gamma#gamma} [MeV]");
      ds_sig.plotOn(myyFr4, Binning(55, 105000, 160000));
      DSCB_myy.plotOn(myyFr4);
      myyFr4->Draw();

      myText(0.22, 0.88, 1, (catName+", signal fit with DSCB").Data());

      c4->SaveAs("plotFit/sigfit_"+catName+"_DSCB.png");
    }

    if(catName.Contains("LT")){
      RooDataHist dh_asi("dh_asi", "", m_yy, Import(*hout_asi));// !!! order
      TCanvas *c5 = new TCanvas("c1", "canvas", 800, 600);
      RooPlot* myyFr5 = m_yy.frame();
      myyFr5->GetXaxis()->SetTitle("m_{#gamma#gamma} [MeV]");
      dh_asi.plotOn(myyFr5, Binning(55, 105000, 160000));
      model_asiv.plotOn(myyFr5);
      myyFr5->Draw();

      myText(0.22, 0.88, 1, (catName+", Asimov data from "+funcname).Data());

      c5->SaveAs("plotFit/Asimov_"+catName+"_DSCB_"+funcname+".png");
    }

//    ofdata.close();
}
//  ofstream ofpara("para_bkg.tmp", ios::out);
//  if(!ofpara){
//    ofpara.close();
//    cout<<"error can't open file for record"<<endl;
//  }
//
//  for(auto p = m_para.begin(); p != m_para.end(); p++){
//    ofpara<<p->first<<","<<p->second[B]<<","<<p->second[C]<<","<<p->second[POWA]<<","<<p->second[YIELD]<<endl;
//  }
//
//  for(auto p = m_para.begin(); p != m_para.end(); p++){
//    cout<<p->first<<","<<p->second[B]<<","<<p->second[C]<<","<<p->second[POWA]<<","<<p->second[YIELD]<<endl;
//  }

  //RooPlot* myy_frame = m_yy.frame(Title("sideband data (polynomial exponential)"));// !!! order
  //dh_SB.plotOn(myy_frame);
  //polyexp.plotOn(myy_frame);
  //polyexp.plotOn(myy_frame, Range(105,160), LineStyle(kDashed));
  //myy_frame->Draw();

  delete f_out;
}
