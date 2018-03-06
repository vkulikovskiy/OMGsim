void calcK40rates(int level = 2, string fname = "coinc.root"){
	double scalek40 = 13750; //in Hz


	ostringstream oss_fname(fname);

	TFile *fK40 = new TFile(oss_fname.str().c_str());
	TH1F *hcoincK40;
	if (level == 2) hcoincK40 = (TH1F *)fK40->Get("hcoincidences");
	else if (level == 3) hcoincK40 = (TH1F *)fK40->Get("hcoincidencesMax");
	else {cout << "wrong level specified! (2 or 3 should be used)" << endl; return; }

	hcoincK40->Scale(scalek40);  //k40 bq = 13750 m-3 sec-1
	hcoincK40->SetLineWidth(2);
	hcoincK40->SetMarkerStyle(20);
	hcoincK40->SetMarkerStyle(4);
	hcoincK40->SetTitle("Rates; coincidence level; rate[Hz]");

	for (int ii = 2; ii < 32; ii++) {
		cout << "level " << hcoincK40->GetBinLowEdge(ii)  << " rate per DOM " << hcoincK40->GetBinContent(ii) << "+-" << hcoincK40->GetBinError(ii) << "Hz" << endl;
	}




}
