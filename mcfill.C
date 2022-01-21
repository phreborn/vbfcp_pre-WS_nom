#include "utils.h"

using namespace TMVA::Experimental;

void mcfill()
{
  map<TString,TString> samples_id;
  samples_id["ggF"] = "343981";
  samples_id["VBF"] = "346214";

  map<TString,float> lumi;
  lumi["mc16a"] = 36207.66;
  lumi["mc16d"] = 44307.4;
  lumi["mc16e"] = 58450.1;

  char *cf_cats = (char*)"cats.cfg";
  //char *cf_bins = (char*)"binnings.cfg";
  map<TString, string> catCuts;
  //map<TString, vector<float>> catBins;
  getCatCuts(cf_cats, catCuts); for(auto c : catCuts) cout<<c.first<<c.second<<endl;
  //getCatBinning(cf_bins, catBins); for(auto c : catBins) { cout<<"dafad"<<c.first<<endl; for(auto b : c.second)  cout<<b<<", "; cout<<endl<<endl; }

  string config = "config_sysNtuple";
  string baseCuts = "";
  string blindCut = "";
  readConfigFile(config.data(), "cuts", baseCuts);
  readConfigFile(config.data(), "blindSel", blindCut);

  baseCuts = TString(baseCuts.data()).ReplaceAll("??", "Nominal.").Data();
  blindCut = TString(blindCut.data()).ReplaceAll("??", "Nominal.").Data();

  //bool doBlind = true;
  //string blind_tmp = "";
  //readConfigFile(config.data(), "doBlind", blind_tmp);
  //if(blind_tmp.find("TRUE") != std::string::npos) doBlind = true;
  //else doBlind = false;

  TString dirpath = "/publicfs/atlas/atlasnew/higgs/hgg/chenhr/vbfcp/ntuple/h026/sys/";
  std::string path_str = dirpath.Data();
  std::vector<std::string> sub_dirs = getDirBinsSortedPath(path_str);

  for(auto id : samples_id){
    TString ID = id.second;
    TString sample = id.first;
    cout<<"processing "<<sample<<".."<<endl;
  for(auto camp : lumi){
    TString mcCamp = camp.first;
    float luminosity = camp.second;

    std::vector<std::string> files(0);

    for(auto d : sub_dirs){
      if(d==".") continue;
      if(d.find(mcCamp.Data()) == std::string::npos) continue;
      if(d.find(".root") != std::string::npos) continue;
      cout<<"d: "<<path_str+d<<endl;
      std::vector<std::string> fs = getDirBinsSortedPath(path_str+d+"/");
      for(auto f : fs){
        if(f==".") continue;
        if(f.find(sample.Data()) == std::string::npos) continue;
        if(f.find("photonsys") == std::string::npos) continue;
        if(f.find(".root") == std::string::npos) continue;
        cout<<"f: "<<path_str+"/"+d+"/"+f<<endl;
        files.push_back(path_str+d+"/"+f);
        
      }
    }

    TFile *fSumWeight = new TFile(files.at(0).data(), "read");
    double sumOfWeights = getSumOfWeights(std::atoi(ID.Data()), fSumWeight); cout<<sumOfWeights<<endl;
    delete fSumWeight;

    cout<<"\n"<<endl;

    for(auto cat : catCuts){
      TString catName = cat.first; cout<<"=== "<<catName<<" ==="<<endl;
      string catCut = cat.second;

      catCut = TString(catCut.data()).ReplaceAll("??", "Nominal.").Data();

      TH1F *h_myy = new TH1F("m_yy_"+sample, "", 110, 105, 160);
      TH1F *h_oo = new TH1F("oo_"+sample, "", 40, -20, 20);

      TChain ch("output", "output");
    
      for(auto f : files){
        TString filepath = f.data();
    
        cout<<filepath<<endl;

        ch.Add(filepath);
 

      }

      ROOT::RDataFrame df(ch, {"Nominal.m_yy"});

      string allCuts = "";
      string allCuts_noBlind = Form("%s && %s", baseCuts.data(), catCut.data()); cout<<"all cuts: "<<allCuts<<endl;
      //if(doBlind) allCuts = Form("%s && %s && %s", baseCuts.data(), catCut.data(), blindCut.data());
      //else allCuts = allCuts_noBlind;
      allCuts = allCuts_noBlind;

      auto df_cut = df.Filter(allCuts);
      auto df_alias = df_cut.Alias("m_yy", "Nominal.m_yy")
                            .Alias("oo1", "Nominal.oo1")
                            .Alias("BDTout_ggH", "Nominal.BDTout_ggH")
                            .Alias("BDTout_yy", "Nominal.BDTout_yy")
                            .Alias("cat_BDTggH_BDTyy", "Nominal.cat_BDTggH_BDTyy");
      auto df_wt = df_alias.Define("weight", [&sumOfWeights, &luminosity](float xsec, float weight, float wjvt, float wfjvt){
        return (float) (luminosity*xsec*weight*wjvt*wfjvt/sumOfWeights);
      }, {"xsec_kF_eff", "Nominal.weight", "Nominal.weightJvt_30", "Nominal.weightFJvt_30"}); cout<<df_wt.Sum("weight").GetValue()<<endl;

      df_wt.Foreach([&h_myy, &h_oo](float myy, float oo, float weight){
        h_myy->Fill(myy/1000, weight);
        h_oo->Fill(oo, weight);
      }, {"m_yy", "oo1", "weight"});

      //df_cut.Snapshot("CollectionTree", Form("tree_%s_OO_%s_%s.root", sample.Data(), catName.Data(), binName.Data()), {"weight", "m_yy", "oo1"});
      df_wt.Snapshot("CollectionTree", Form("MC16Xs/%s/hists_%s_%s.root", mcCamp.Data(), sample.Data(), catName.Data()), {"weight", "m_yy", "oo1", "BDTout_ggH", "BDTout_yy", "cat_BDTggH_BDTyy"});
    
      TFile *f_out = new TFile(Form("MC16Xs/%s/hists_%s_%s.root", mcCamp.Data(), sample.Data(), catName.Data()),"update");

      f_out->cd();
      h_myy->Write();
      h_oo->Write();
      //f_out->Close();
    
      delete h_myy;
      delete h_oo;

    }
  }
  }
}
