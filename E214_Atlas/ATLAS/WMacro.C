{
double xmin = 32;
double xmax = 50;
w.SetCutSelection("njet==0 && etmis>35 && el_etiso>1 && el_energy > 40");
w.SetQCDScaleFactor(1.85/3.2*192/225);
TCanvas* c1 = new TCanvas("c1");
c1->Clear();
c1->Divide(3,3);
c1->cd(1);w.GetHalfMaximum("MCW75",xmin,xmax);
c1->cd(2);w.GetHalfMaximum("MCW78",xmin,xmax);
c1->cd(3);w.GetHalfMaximum("MCW79",xmin,xmax);
c1->cd(4);w.GetHalfMaximum("MCW80",xmin,xmax);
c1->cd(5);w.GetHalfMaximum("MCW81",xmin,xmax);
c1->cd(6);w.GetHalfMaximum("MCW82",xmin,xmax);
c1->cd(7);w.GetHalfMaximum("MCW85",xmin,xmax);
c1->cd(8);w.GetHalfMaximum("Wenu",xmin,xmax);
c1->cd(9);w.GetHalfMaximum("Zee",xmin,xmax);
}
