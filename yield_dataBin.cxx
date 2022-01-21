void yield_dataBin(TString file){
  TFile *f = new TFile(file, "read");
  TTree *t = (TTree*) f->Get("CollectionTree");
  cout<<t->GetEntries()<<endl;
}
