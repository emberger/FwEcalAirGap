{
//=========Macro generated from canvas: Comparison/Comparison
//=========  (Fri Dec  8 16:44:02 2017) by ROOT version5.34/36
   TCanvas *Comparison = new TCanvas("Comparison", "Comparison",0,0,700,500);
   Comparison->Range(0,0,1,1);
   Comparison->SetFillColor(0);
   Comparison->SetBorderMode(0);
   Comparison->SetBorderSize(2);
   Comparison->SetFrameBorderMode(0);
   
   TGraph *graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("Res_PCA_X");
   graph->SetFillColor(1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.3066009539);
   graph->SetPoint(1,150,0.1576959847);
   graph->SetPoint(2,250,0.1238425196);
   graph->SetPoint(3,350,0.1058719917);
   graph->SetPoint(4,450,0.09681352711);
   graph->SetPoint(5,550,0.09374884517);
   graph->SetPoint(6,650,0.07918247035);
   graph->SetPoint(7,750,0.07755779967);
   graph->SetPoint(8,850,0.07715525799);
   graph->SetPoint(9,950,0.06937840799);
   graph->SetPoint(10,1500,0.06025170665);
   
   TH1F *Graph_Graph115 = new TH1F("Graph_Graph115","Res_PCA_X",100,0,1645);
   Graph_Graph115->SetMinimum(0);
   Graph_Graph115->SetMaximum(0.5);
   Graph_Graph115->SetDirectory(0);
   Graph_Graph115->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph115->SetLineColor(ci);
   Graph_Graph115->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph115->GetXaxis()->SetLabelFont(42);
   Graph_Graph115->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph115->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph115->GetXaxis()->SetTitleFont(42);
   Graph_Graph115->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");
   Graph_Graph115->GetYaxis()->SetLabelFont(42);
   Graph_Graph115->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph115->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph115->GetYaxis()->SetTitleOffset(1.23);
   Graph_Graph115->GetYaxis()->SetTitleFont(42);
   Graph_Graph115->GetZaxis()->SetLabelFont(42);
   Graph_Graph115->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph115->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph115->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph115);
   
   graph->Draw("ap");
   
   graph = new TGraph(11);
   graph->SetName("Graph1");
   graph->SetTitle("Res_PCA_X");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.4023719196);
   graph->SetPoint(1,150,0.1704925996);
   graph->SetPoint(2,250,0.1386163224);
   graph->SetPoint(3,350,0.1197869949);
   graph->SetPoint(4,450,0.1141012707);
   graph->SetPoint(5,550,0.108959119);
   graph->SetPoint(6,650,0.0953261878);
   graph->SetPoint(7,750,0.09067575531);
   graph->SetPoint(8,850,0.08647029635);
   graph->SetPoint(9,950,0.07887780217);
   graph->SetPoint(10,1500,0.06553743437);
   
   TH1F *Graph_Graph116 = new TH1F("Graph_Graph116","Res_PCA_X",100,0,1645);
   Graph_Graph116->SetMinimum(0.03185399);
   Graph_Graph116->SetMaximum(0.4360554);
   Graph_Graph116->SetDirectory(0);
   Graph_Graph116->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph116->SetLineColor(ci);
   Graph_Graph116->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph116->GetXaxis()->SetLabelFont(42);
   Graph_Graph116->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph116->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph116->GetXaxis()->SetTitleFont(42);
   Graph_Graph116->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");
   Graph_Graph116->GetYaxis()->SetLabelFont(42);
   Graph_Graph116->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph116->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph116->GetYaxis()->SetTitleFont(42);
   Graph_Graph116->GetZaxis()->SetLabelFont(42);
   Graph_Graph116->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph116->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph116->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph116);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph2");
   graph->SetTitle("Res_PCA_X");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.3099820229);
   graph->SetPoint(1,150,0.1536431945);
   graph->SetPoint(2,250,0.1248846198);
   graph->SetPoint(3,350,0.1053524671);
   graph->SetPoint(4,450,0.09725702241);
   graph->SetPoint(5,550,0.09481908917);
   graph->SetPoint(6,650,0.08860632242);
   graph->SetPoint(7,750,0.08229247296);
   graph->SetPoint(8,850,0.07646919615);
   graph->SetPoint(9,950,0.07413633305);
   graph->SetPoint(10,1500,0.06228128951);
   
   TH1F *Graph_Graph117 = new TH1F("Graph_Graph117","Res_PCA_X",100,0,1645);
   Graph_Graph117->SetMinimum(0.03751122);
   Graph_Graph117->SetMaximum(0.3347521);
   Graph_Graph117->SetDirectory(0);
   Graph_Graph117->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph117->SetLineColor(ci);
   Graph_Graph117->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph117->GetXaxis()->SetLabelFont(42);
   Graph_Graph117->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph117->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph117->GetXaxis()->SetTitleFont(42);
   Graph_Graph117->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");
   Graph_Graph117->GetYaxis()->SetLabelFont(42);
   Graph_Graph117->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph117->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph117->GetYaxis()->SetTitleFont(42);
   Graph_Graph117->GetZaxis()->SetLabelFont(42);
   Graph_Graph117->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph117->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph117->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph117);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph3");
   graph->SetTitle("Res_PCA_X");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.4890764494);
   graph->SetPoint(1,150,0.1782306377);
   graph->SetPoint(2,250,0.1472157976);
   graph->SetPoint(3,350,0.1165703362);
   graph->SetPoint(4,450,0.1115127118);
   graph->SetPoint(5,550,0.1050886375);
   graph->SetPoint(6,650,0.09545434827);
   graph->SetPoint(7,750,0.0956608737);
   graph->SetPoint(8,850,0.08632719159);
   graph->SetPoint(9,950,0.08588714829);
   graph->SetPoint(10,1500,0.06944928491);
   
   TH1F *Graph_Graph118 = new TH1F("Graph_Graph118","Res_PCA_X",100,0,1645);
   Graph_Graph118->SetMinimum(0.02748657);
   Graph_Graph118->SetMaximum(0.5310392);
   Graph_Graph118->SetDirectory(0);
   Graph_Graph118->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph118->SetLineColor(ci);
   Graph_Graph118->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph118->GetXaxis()->SetLabelFont(42);
   Graph_Graph118->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph118->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph118->GetXaxis()->SetTitleFont(42);
   Graph_Graph118->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");
   Graph_Graph118->GetYaxis()->SetLabelFont(42);
   Graph_Graph118->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph118->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph118->GetYaxis()->SetTitleFont(42);
   Graph_Graph118->GetZaxis()->SetLabelFont(42);
   Graph_Graph118->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph118->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph118->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph118);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph4");
   graph->SetTitle("Res_PCA_X");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(25);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.2858400207);
   graph->SetPoint(1,150,0.151471972);
   graph->SetPoint(2,250,0.1208378748);
   graph->SetPoint(3,350,0.1094803344);
   graph->SetPoint(4,450,0.1003672397);
   graph->SetPoint(5,550,0.08462008387);
   graph->SetPoint(6,650,0.08670357601);
   graph->SetPoint(7,750,0.0797003686);
   graph->SetPoint(8,850,0.07485705641);
   graph->SetPoint(9,950,0.07555611064);
   graph->SetPoint(10,1500,0.05681965661);
   
   TH1F *Graph_Graph119 = new TH1F("Graph_Graph119","Res_PCA_X",100,0,1645);
   Graph_Graph119->SetMinimum(0.03391762);
   Graph_Graph119->SetMaximum(0.3087421);
   Graph_Graph119->SetDirectory(0);
   Graph_Graph119->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph119->SetLineColor(ci);
   Graph_Graph119->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph119->GetXaxis()->SetLabelFont(42);
   Graph_Graph119->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph119->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph119->GetXaxis()->SetTitleFont(42);
   Graph_Graph119->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");
   Graph_Graph119->GetYaxis()->SetLabelFont(42);
   Graph_Graph119->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph119->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph119->GetYaxis()->SetTitleFont(42);
   Graph_Graph119->GetZaxis()->SetLabelFont(42);
   Graph_Graph119->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph119->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph119->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph119);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph5");
   graph->SetTitle("Res_PCA_X");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.4893560877);
   graph->SetPoint(1,150,0.1643010135);
   graph->SetPoint(2,250,0.1424274717);
   graph->SetPoint(3,350,0.1216583118);
   graph->SetPoint(4,450,0.1075570539);
   graph->SetPoint(5,550,0.1038292376);
   graph->SetPoint(6,650,0.09220496767);
   graph->SetPoint(7,750,0.09021965264);
   graph->SetPoint(8,850,0.08703734938);
   graph->SetPoint(9,950,0.08570877474);
   graph->SetPoint(10,1500,0.07016498413);
   
   TH1F *Graph_Graph120 = new TH1F("Graph_Graph120","Res_PCA_X",100,0,1645);
   Graph_Graph120->SetMinimum(0.02824587);
   Graph_Graph120->SetMaximum(0.5312752);
   Graph_Graph120->SetDirectory(0);
   Graph_Graph120->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph120->SetLineColor(ci);
   Graph_Graph120->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph120->GetXaxis()->SetLabelFont(42);
   Graph_Graph120->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph120->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph120->GetXaxis()->SetTitleFont(42);
   Graph_Graph120->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");
   Graph_Graph120->GetYaxis()->SetLabelFont(42);
   Graph_Graph120->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph120->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph120->GetYaxis()->SetTitleFont(42);
   Graph_Graph120->GetZaxis()->SetLabelFont(42);
   Graph_Graph120->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph120->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph120->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph120);
   
   graph->Draw("p");
   
   TLegend *leg = new TLegend(-2.353437e-185,-2.353437e-185,-2.353437e-185,-2.353437e-185,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph0","1mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph1","2mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph2","1mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph3","2mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph4","1mmLeadAbsorber_NOVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph5","2mmLeadAbsorber_NOVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   Comparison->Modified();
   Comparison->cd();
   Comparison->SetSelected(Comparison);
}
