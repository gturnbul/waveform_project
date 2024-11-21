void sndisplay_calorimeter_test_values (bool with_palette = true)
{
  sncalo = new sndisplay::calorimeter ("sndiplay_test", with_palette);
  sncalo->draw_content_label("%.1f");

  // open tree
  TFile *File = TFile::Open("treebuild.root", "READ");
  TTree *tree = (TTree*)File->Get("treebuild");

  int max_om  = 520;

  // Use vectors to store the sum of chi2ndf values and the count of entries for each OM number
  std::vector<double> chi2ndf_sum(max_om, 0.0); // Initialize all elements to 0.0
  std::vector<int> chi2ndf_count(max_om, 0);    // Initialize all elements to 0

  // set up variables to hold the branch data
  int e_om_number;
  double e_chi2ndf;
  double e_event;

  //set branch address
  tree->SetBranchAddress("e_om_number", &e_om_number);
  tree->SetBranchAddress("e_chi2ndf", &e_chi2ndf);
  tree->SetBranchAddress("e_event", &e_event);

  // Loop over the entries in the tree to fill the arrays
  int nentries = tree->GetEntries();
  for (int i = 0; i < nentries; i++) {
      tree->GetEntry(i);

    // Get sums for each om number 
      if (e_event > -0.1 && e_om_number >= 0 && e_om_number < max_om) {
          chi2ndf_sum[e_om_number] += e_chi2ndf;
          chi2ndf_count[e_om_number]++;
      }
  }

  for (int omnum=0; omnum<520; ++omnum) // MW
    if(chi2ndf_count[omnum]>0){
      //double average = chi2ndf_sum[omnum]/chi2ndf_count[omnum];
      double counts = chi2ndf_count[omnum];
      sncalo->setcontent(omnum, counts );
    } else {
      sncalo->setcontent(omnum,0);
    }

  for (int omnum=520; omnum<648; ++omnum) // XW
    sncalo->setcontent(omnum, 0);

  for (int omnum=648; omnum<712; ++omnum) // GV
    sncalo->setcontent(omnum, 0);

  sncalo->setrange(0, 100); // to force z axis range

  sncalo->draw();

  sncalo->canvas_it->SaveAs("sndisplay-calorimeter-test-it.png");
  sncalo->canvas_fr->SaveAs("sndisplay-calorimeter-test-fr.png");

  // merge IT and FR canvas side by side using image magick (if installed)
  gSystem->Exec("which convert > /dev/null && convert sndisplay-calorimeter-test-it.png sndisplay-calorimeter-test-fr.png +append sndisplay-calorimeter-test.png");
}
