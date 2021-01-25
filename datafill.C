#include "utils.h"

void datafill()
{

  char *cf_cats = (char*)"cats.cfg";
  char *cf_bins = (char*)"binnings.cfg";
  map<TString, string> catCuts;
  map<TString, vector<float>> catBins;
  getCatCuts(cf_cats, catCuts); for(auto c : catCuts) cout<<c.first<<c.second<<endl;
  getCatBinning(cf_bins, catBins); for(auto c : catBins) { cout<<"dafad"<<c.first<<endl; for(auto b : c.second)  cout<<b<<", "; cout<<endl<<endl; }

  string config = "config";
  string baseCuts = "";
  string blindCut = "";
  readConfigFile(config.data(), "cuts", baseCuts);
  readConfigFile(config.data(), "blindSel", blindCut);

  bool doBlind = true;
  string blind_tmp = "";
  readConfigFile(config.data(), "doBlind", blind_tmp);
  if(blind_tmp.find("TRUE") != std::string::npos) doBlind = true;
  else doBlind = false;

  map<TString, pair<float, float>> bins;
////  bins[""] = make_pair(-999999999, 999999999);
//  bins["b1"] = make_pair(-999999999, -3);
//  bins["b2"] = make_pair(-3, -1.5);
//  bins["b3"] = make_pair(-1.5, -1);
//  bins["b4"] = make_pair(-1, -0.5);
//  bins["b5"] = make_pair(-0.5, 0);
//  bins["b6"] = make_pair(0, 0.5);
//  bins["b7"] = make_pair(0.5, 1);
//  bins["b8"] = make_pair(1, 1.5);
//  bins["b9"] = make_pair(1.5, 3);
//  bins["b10"] = make_pair(3, 99999999);
//  bins["b1"] = make_pair(-999999999, -2);
//  bins["b2"] = make_pair(-2, -1);
//  bins["b3"] = make_pair(-1, 0);
//  bins["b4"] = make_pair(0, 1);
//  bins["b5"] = make_pair(1, 2);
//  bins["b6"] = make_pair(2, 99999999);
  bins["b1"] = make_pair(-999999999, -2);
  bins["b2"] = make_pair(-2, -1);
  bins["b3"] = make_pair(-1, 0);
  bins["b4"] = make_pair(0, 1);
  bins["b5"] = make_pair(1, 2);
  bins["b6"] = make_pair(2, 99999999);

  std::map<TString, std::vector<float>> cats;
  cats["TT"] = {0.14, 1., 0.23, 1.};
  cats["TL"] = {0.14, 1., -1., 0.23};
  cats["LT"] = {-1, 0.14, 0.05, 1.};
  cats["LL"] = {-1, 0.14, -1., 0.05};

  //TString dirpath = "/scratchfs/bes/chenhr/atlaswork/VBF_CP/ntuples/";
  //TString dirpath = "/workfs/bes/chenhr/atlas/vbf_cp/ntuples/";//
  TString dirpath = "/scratchfs/bes/chenhr/atlaswork/VBF_CP/h025/";
  std::string path_str = dirpath.Data();
  std::vector<std::string> sub_dirs = getDirBinsSortedPath(path_str);

  std::vector<std::string> files(0);

  for(auto d : sub_dirs){
    if(d==".") continue;
    if(d.find("data") == std::string::npos) continue;
    if(d.find(".root") != std::string::npos) continue;
    cout<<"d: "<<path_str+d<<endl;
    std::vector<std::string> fs = getDirBinsSortedPath(path_str+d+"/");
    for(auto f : fs){
      if(f==".") continue;
      if(f.find("slim") != std::string::npos) continue;
      cout<<"f: "<<path_str+"/"+d+"/"+f<<endl;
      files.push_back(path_str+d+"/"+f);
      
    }
  }

  cout<<"\n"<<endl;

  TString sample = "dataSB";
//  TString sample = "fulldata";

  cout<<"processing "<<sample<<".."<<endl;

  Int_t N_j,N_photon,cutflow,Category;
  Float_t m_yy,pT_y1,pT_y2,m_jj_30,DeltaEta_jj,Zepp,oo1,oo2,WeightDtilde1,WeightDtilde2,weight,xsec_kF_eff,total_weights;
  Bool_t isDalitz,isPassed,isPassedIsolation,isPassedPID,isPassedTriggerMatch;

for(auto cat : catCuts){
  TString catName = cat.first; cout<<"=== "<<catName<<" ==="<<endl;
  string catCut = cat.second;

  vector<float> binEdges = catBins[catName];
  for(int i = 0; i < binEdges.size()-1; i++){
    TString binName = Form("b%i", i+1); cout<<binName<<endl;

    double b_l = binEdges.at(i);
    double b_r = binEdges.at(i+1);

    TH1F *h_myy = new TH1F("m_yy_"+sample, "", 110, 105, 160);
    TH1F *h_oo = new TH1F("oo_"+sample, "", 40, -20, 20);

    TTree *t_out = new TTree("CollectionTree", "");

//    double myy_out, oo_out, weight_out;
//    t_out->Branch("m_yy", &myy_out, "m_yy/D");
//    t_out->Branch("oo1", &oo_out, "oo1/D");
//    t_out->Branch("weight", &weight_out, "weight/D");

    TChain ch("output", "output");
  
    for(auto f : files){
      TString filepath = f.data();
//      TFile *f_in = new TFile(filepath,"read");
//      TTree *tree = (TTree*) f_in->Get("output");
  
      cout<<filepath<<endl;

      ch.Add(filepath);
 
//      tree->SetBranchAddress("cutflow", &cutflow);
//      //tree->SetBranchAddress("Category", &Category);//
//      tree->SetBranchAddress("catCoup_XGBoost_ttH", &Category);
//      tree->SetBranchAddress("isDalitz", &isDalitz);
//      tree->SetBranchAddress("isPassed", &isPassed);
//      tree->SetBranchAddress("isPassedIsolation", &isPassedIsolation);
//      tree->SetBranchAddress("isPassedPID", &isPassedPID);
//      tree->SetBranchAddress("isPassedTriggerMatch", &isPassedTriggerMatch);
//      tree->SetBranchAddress("m_yy", &m_yy);
//      //tree->SetBranchAddress("N_j", &N_j);//
//      tree->SetBranchAddress("N_j_30", &N_j);
//      tree->SetBranchAddress("N_photon", &N_photon);
//      tree->SetBranchAddress("pT_y1", &pT_y1);
//      tree->SetBranchAddress("pT_y2", &pT_y2);
//      tree->SetBranchAddress("m_jj_30", &m_jj_30);
//      tree->SetBranchAddress("DeltaEta_jj", &DeltaEta_jj);
//      tree->SetBranchAddress("Zepp", &Zepp);
//      tree->SetBranchAddress("oo1", &oo1);
//      tree->SetBranchAddress("oo2", &oo2);
//  
//  
//      Long64_t endentry = tree->GetEntries();
//  
//      for(int i = 0; i < endentry; i++){
//        tree->GetEntry(i);
//        if(i%100000==0&&i!=0) std::cout<<i<<" events processed"<<std::endl;
//        if(i==endentry-1) cout<<endentry<<" events processed"<<endl;
//        if(isDalitz==1) continue;
//        if(isPassed==0) continue;
//        if(N_j<2) continue;
//        if(m_jj_30/1000<400) continue;
//        if(DeltaEta_jj>-2&&DeltaEta_jj<2) continue;
//        if(Zepp>5||Zepp<-5) continue;
//        if(Category<11||Category>14) continue;
//  
//        if(m_yy/1000>120&&m_yy/1000<130) continue;
//
//        if(oo1<b_l||oo1>=b_r) continue; 
//
//        //myy_out = m_yy/1000;
//        myy_out = m_yy;
//        oo_out = oo1;
//        weight_out = 1.;
//        t_out->Fill();
// 
//        h_myy->Fill(m_yy/1000);
//        h_oo->Fill(oo1);
//  
//      }
//      delete f_in;

    }

    string binCut = Form("oo1 >= %f && oo1 < %f", b_l, b_r);
    string allCuts = "";
    if(doBlind) allCuts = Form("%s && %s && %s && %s", baseCuts.data(), catCut.data(), binCut.data(), blindCut.data());
    else allCuts = Form("%s && %s && %s", baseCuts.data(), catCut.data(), binCut.data()); cout<<"all cuts: "<<allCuts<<endl;

    ROOT::RDataFrame df(ch, {"m_yy"});
    auto df_cut = df.Filter(allCuts);
    df_cut.Foreach([&h_myy, &h_oo](float myy, float oo){ h_myy->Fill(myy/1000); h_oo->Fill(oo); }, {"m_yy", "oo1"}); cout<<df_cut.Sum("weight").GetValue()<<endl;
    //df_cut.Snapshot("CollectionTree", Form("tree_%s_OO_%s_%s.root", sample.Data(), catName.Data(), binName.Data()), {"weight", "m_yy", "oo1"});
    df_cut.Snapshot("CollectionTree", Form("hists_%s_%s_%s.root", sample.Data(), catName.Data(), binName.Data()), {"weight", "m_yy", "oo1"});
  
    TFile *f_out = new TFile(Form("hists_%s_%s_%s.root", sample.Data(), catName.Data(), binName.Data()),"update");

    f_out->cd();
    //t_out->Write();
    h_myy->Write();
    h_oo->Write();
    //f_out->Close();
  
    delete h_myy;
    delete h_oo;
  }
}

}
