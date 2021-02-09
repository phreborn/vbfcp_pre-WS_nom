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

void getBkgParaYield()
{
  char *cf_cats = (char*)"cats.cfg";
  char *cf_bins = (char*)"binnings.cfg";
  map<TString, string> catCuts;
  map<TString, vector<float>> catBins;
  getCatCuts(cf_cats, catCuts); for(auto c : catCuts) cout<<c.first<<c.second<<endl;
  getCatBinning(cf_bins, catBins); for(auto c : catBins) { cout<<"dafad"<<c.first<<endl; for(auto b : c.second)  cout<<b<<", "; cout<<endl<<endl; }

  map<TString, pair<float, float>> bins;
  //bins["b1"] = make_pair(-999999999, -3);
  //bins["b2"] = make_pair(-3, -1.5);
  //bins["b3"] = make_pair(-1.5, -1);
  //bins["b4"] = make_pair(-1, -0.5);
  //bins["b5"] = make_pair(-0.5, 0);
  //bins["b6"] = make_pair(0, 0.5);
  //bins["b7"] = make_pair(0.5, 1);
  //bins["b8"] = make_pair(1, 1.5);
  //bins["b9"] = make_pair(1.5, 3);
  //bins["b10"] = make_pair(3, 99999999);
  bins["b1"] = make_pair(-999999999, -2);
  bins["b2"] = make_pair(-2, -1);
  bins["b3"] = make_pair(-1, 0);
  bins["b4"] = make_pair(0, 1);
  bins["b5"] = make_pair(1, 2);
  bins["b6"] = make_pair(2, 99999999);

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

  enum Para{B, C, YIELD};
  map<TString, vector<float>> m_para;

  double myy_out, oo_out, weight_out;

  TFile *f_out = new TFile("tree_asimov.root", "recreate");
  //TFile *f_out = new TFile("tree_asimov_minus_0010.root", "recreate");

for(auto cat : catCuts){
  TString catName = cat.first; cout<<"=== "<<catName<<" ==="<<endl;
  string catCut = cat.second;

  vector<float> binEdges = catBins[catName];
  for(int i = 0; i < binEdges.size()-1; i++){
    TCanvas *canv = new TCanvas("c", "canvas", 800, 600);

    TString binName = Form("b%i", i+1); cout<<binName<<endl;

    double b_l = binEdges.at(i);
    double b_r = binEdges.at(i+1);

    cout<<"=============="<<catName<<"_"<<binName<<"=============="<<endl;

    TTree *t_out = new TTree("CollectionTree10"+bins_name[binName], "");
  
    //t_out->Branch("m_yy_SR", &myy_out, "m_yy_SR/D");
    t_out->Branch("m_yy", &myy_out, "m_yy/D");
    //t_out->Branch("oo1", &oo_out, "oo1/D");
    t_out->Branch("weight", &weight_out, "weight/D");
  
    TFile *f_SB = new TFile(Form("hists_dataSB_%s_%s.root", catName.Data(), binName.Data()), "read");
    TTree *t_SB = (TTree*) f_SB->Get("CollectionTree");
    TH1F *h_SB = (TH1F*) f_SB->Get("m_yy_dataSB");
    float n_SB = h_SB->Integral();
  
    TFile *f_ggh = new TFile(Form("hists_ggh_%s_%s.root", catName.Data(), binName.Data()), "read");
    TTree *t_ggh = (TTree*) f_ggh->Get("CollectionTree");
    TH1F *h_ggh = (TH1F*) f_ggh->Get("m_yy_ggh");
    float n_ggh = h_ggh->Integral();
  
    TFile *f_vbf = new TFile(Form("hists_vbf_%s_%s.root", catName.Data(), binName.Data()), "read");
    //TFile *f_vbf = new TFile("hists_vbf_minus_0010_"+bin->first+".root", "read");
    TTree *t_vbf = (TTree*) f_vbf->Get("CollectionTree");
    TH1F *h_vbf = (TH1F*) f_vbf->Get("m_yy_vbf");
    float n_vbf = h_vbf->Integral();
  
    TH1F h_sig = *h_ggh+*h_vbf;
    float n_sig = n_ggh+n_vbf;
  
    float SR_up = 120000;
    float SR_lo = 130000;
  
    RooRealVar m_yy("m_yy", "m_yy", 105000, 160000);
    RooRealVar oo1("oo1", "oo1", -20, 20);
    RooRealVar weight("weight", "weight", -10, 100);
    RooPlot* oo_frame = oo1.frame(Title("oo projection of DSCB(oo)*polyexp(m_yy)")) ;
  
//    RooDataHist dh_SB("dh_SB", "", m_yy, Import(*h_SB));// !!! order
//    RooDataHist dh_sig("dh_sig", "", m_yy, Import(h_sig));
 
    RooDataSet ds_SB("ds_SB","", RooArgSet(m_yy, weight), WeightVar(weight), Import(*t_SB));
    RooDataSet ds_ggh("ds_ggh","", RooArgSet(m_yy, weight), WeightVar(weight), Import(*t_ggh));
    RooDataSet ds_vbf("ds_vbf","", RooArgSet(m_yy, weight), WeightVar(weight), Import(*t_vbf));
    ds_ggh.append(ds_vbf);
    RooDataSet ds_sig = ds_ggh;
 
    RooRealVar b("b", "", -12, 12);
    RooRealVar c("c", "", -15, 15);
    //RooGenericPdf polyexp("polyexp", "fit sideband m_yy", "exp(b*(m_yy/1000-125)/125)", RooArgSet(b, m_yy));// Exp
    RooGenericPdf polyexp("polyexp", "fit sideband m_yy", "exp(b*(m_yy/1000-125)/125+c*((m_yy/1000-125)/125)*((m_yy/1000-125)/125))", RooArgSet(b, c, m_yy));// ExpPoly2
  
    m_yy.setRange("SR", SR_up, SR_lo);
    m_yy.setRange("SB1", 105000, SR_up);
    m_yy.setRange("SB2", SR_lo, 160000);
  
    //polyexp.fitTo(dh_SB, Range("SB1,SB2"), Save());
    polyexp.fitTo(ds_SB, Range("SB1,SB2"), Save());

    m_para[catName+"_"+binName].push_back(b.getVal());
    //m_para[bin->first].push_back(0.1);// Exp
    m_para[catName+"_"+binName].push_back(c.getVal());// ExpPoly2

    RooPlot* myy_frame = m_yy.frame(Title("fit to m_yy distribution"));
    //ds_SB.plotOn(myy_frame);
    //polyexp.plotOn(myy_frame);
//    myy_frame->Draw();
  
    RooAbsReal *int_SB1 = polyexp.createIntegral(m_yy, NormSet(m_yy), Range("SB1"));
    RooAbsReal *int_SB2 = polyexp.createIntegral(m_yy, NormSet(m_yy), Range("SB2"));
    RooAbsReal *int_SR = polyexp.createIntegral(m_yy, NormSet(m_yy), Range("SR"));
    float frac_SR_SB = int_SR->getVal()/(int_SB1->getVal()+int_SB2->getVal());
    cout<<"int_SR: "<<int_SR->getVal()<<endl;
    cout<<"int_SB1: "<<int_SB1->getVal()<<endl;
    cout<<"int_SB2: "<<int_SB2->getVal()<<endl;
  
    float n_SR = n_SB*frac_SR_SB;cout<<"n_SR: "<<n_SR<<endl;cout<<"n_SB: "<<n_SB<<endl;
    float n_bkg = n_SB+n_SR;
 
    m_para[catName+"_"+binName].push_back(n_bkg);
 
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

    canv->SaveAs("plotFit/"+catName+"_"+binName+".png");
    delete canv;

    RooRealVar N_sig("N_sig","", n_sig);
    RooRealVar N_bkg("N_bkg","", n_bkg);
    RooAddPdf model_asiv("model_asiv","", RooArgList(DSCB_myy, polyexp), RooArgList(N_sig, N_bkg));

//    ofstream ofdata("data_OO_"+bin->first+".txt",ios::out);
//    if(!ofdata){
//      ofdata.close();
//      cout<<"can't open file for writing"<<endl;
//    }
  
    int n_points = 550;
    double interval = (160000.-105000.)/n_points;
  
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
      intPdf += pdf_val*interval;

//      for(int i = 0; i < (int)(weight_out+0.5); i++){
//        ofdata<<myy_out<<endl;
//      }
    }
    cout<<"intPdf: "<<intPdf<<endl;
    cout<<"n_ggh, n_vbf: "<<n_ggh<<", "<<n_vbf<<endl;
    cout<<"n_SB, n_sig: "<<n_SB<<", "<<n_sig<<endl;
    cout<<"n_SR: "<<n_SR<<endl;

//    ofdata.close();

    f_out->cd();
    t_out->Write();

    //f_out->Close();

    delete t_out;
  }// end bin
}
  ofstream ofpara("para_bkg.csv", ios::out);
  if(!ofpara){
    ofpara.close();
    cout<<"error can't open file for record"<<endl;
  }

  for(auto p = m_para.begin(); p != m_para.end(); p++){
    ofpara<<p->first<<","<<p->second[B]<<","<<p->second[C]<<","<<p->second[YIELD]<<endl;
  }

  for(auto p = m_para.begin(); p != m_para.end(); p++){
    cout<<p->first<<","<<p->second[B]<<","<<p->second[C]<<","<<p->second[YIELD]<<endl;
  }

  //RooPlot* myy_frame = m_yy.frame(Title("sideband data (polynomial exponential)"));// !!! order
  //dh_SB.plotOn(myy_frame);
  //polyexp.plotOn(myy_frame);
  //polyexp.plotOn(myy_frame, Range(105,160), LineStyle(kDashed));
  //myy_frame->Draw();

  delete f_out;
}
