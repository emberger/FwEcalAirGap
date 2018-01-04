{
//=========Macro generated from canvas: Comparison/Comparison
//=========  (Fri Dec  8 23:48:53 2017) by ROOT version5.34/36
   TCanvas *Comparison = new TCanvas("Comparison", "Comparison",0,0,700,500);
   Comparison->Range(0,0,1,1);
   Comparison->SetFillColor(0);
   Comparison->SetBorderMode(0);
   Comparison->SetBorderSize(2);
   Comparison->SetFrameBorderMode(0);
   
   TGraph *graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("GeometryComparison_PCA");
   graph->SetFillColor(1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.5356740741);
   graph->SetPoint(1,150,0.2668055556);
   graph->SetPoint(2,250,0.2038815789);
   graph->SetPoint(3,350,0.1713333333);
   graph->SetPoint(4,450,0.1494444444);
   graph->SetPoint(5,550,0.135);
   graph->SetPoint(6,650,0.125);
   graph->SetPoint(7,750,0.1137588652);
   graph->SetPoint(8,850,0.108649635);
   graph->SetPoint(9,950,0.1024137931);
   graph->SetPoint(10,1500,0.08259794872);
   
   TH1F *Graph_Graph7 = new TH1F("Graph_Graph7","GeometryComparison_PCA",100,0,1645);
   Graph_Graph7->SetMinimum(0);
   Graph_Graph7->SetMaximum(0.5);
   Graph_Graph7->SetDirectory(0);
   Graph_Graph7->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph7->SetLineColor(ci);
   Graph_Graph7->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph7->GetXaxis()->SetLabelFont(42);
   Graph_Graph7->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph7->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph7->GetXaxis()->SetTitleFont(42);
   Graph_Graph7->GetYaxis()->SetTitle("AngularResolution[rad]");
   Graph_Graph7->GetYaxis()->SetLabelFont(42);
   Graph_Graph7->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph7->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph7->GetYaxis()->SetTitleOffset(1.23);
   Graph_Graph7->GetYaxis()->SetTitleFont(42);
   Graph_Graph7->GetZaxis()->SetLabelFont(42);
   Graph_Graph7->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph7->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph7->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph7);
   
   graph->Draw("ap");
   
   graph = new TGraph(11);
   graph->SetName("Graph1");
   graph->SetTitle("GeometryComparison_PCA");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.62264);
   graph->SetPoint(1,150,0.316);
   graph->SetPoint(2,250,0.2359322034);
   graph->SetPoint(3,350,0.2045762712);
   graph->SetPoint(4,450,0.1804477612);
   graph->SetPoint(5,550,0.1643548387);
   graph->SetPoint(6,650,0.1510416667);
   graph->SetPoint(7,750,0.1391666667);
   graph->SetPoint(8,850,0.1325373134);
   graph->SetPoint(9,950,0.1255045872);
   graph->SetPoint(10,1500,0.09872670807);
   
   TH1F *Graph_Graph8 = new TH1F("Graph_Graph8","Res_68Quantile_PCA",100,0,1645);
   Graph_Graph8->SetMinimum(0.04633538);
   Graph_Graph8->SetMaximum(0.6750313);
   Graph_Graph8->SetDirectory(0);
   Graph_Graph8->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph8->SetLineColor(ci);
   Graph_Graph8->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph8->GetXaxis()->SetLabelFont(42);
   Graph_Graph8->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph8->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph8->GetXaxis()->SetTitleFont(42);
   Graph_Graph8->GetYaxis()->SetTitle("AngularResolution");
   Graph_Graph8->GetYaxis()->SetLabelFont(42);
   Graph_Graph8->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph8->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph8->GetYaxis()->SetTitleFont(42);
   Graph_Graph8->GetZaxis()->SetLabelFont(42);
   Graph_Graph8->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph8->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph8->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph8);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph2");
   graph->SetTitle("Res_68Quantile_PCA");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.5614631579);
   graph->SetPoint(1,150,0.2731466667);
   graph->SetPoint(2,250,0.2086619718);
   graph->SetPoint(3,350,0.1778846154);
   graph->SetPoint(4,450,0.1505333333);
   graph->SetPoint(5,550,0.1403846154);
   graph->SetPoint(6,650,0.1276666667);
   graph->SetPoint(7,750,0.1214473684);
   graph->SetPoint(8,850,0.1122764228);
   graph->SetPoint(9,950,0.1043382353);
   graph->SetPoint(10,1500,0.08428977273);
   
   TH1F *Graph_Graph9 = new TH1F("Graph_Graph9","Res_68Quantile_PCA",100,0,1645);
   Graph_Graph9->SetMinimum(0.03657243);
   Graph_Graph9->SetMaximum(0.6091805);
   Graph_Graph9->SetDirectory(0);
   Graph_Graph9->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph9->SetLineColor(ci);
   Graph_Graph9->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph9->GetXaxis()->SetLabelFont(42);
   Graph_Graph9->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph9->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph9->GetXaxis()->SetTitleFont(42);
   Graph_Graph9->GetYaxis()->SetTitle("AngularResolution");
   Graph_Graph9->GetYaxis()->SetLabelFont(42);
   Graph_Graph9->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph9->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph9->GetYaxis()->SetTitleFont(42);
   Graph_Graph9->GetZaxis()->SetLabelFont(42);
   Graph_Graph9->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph9->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph9->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph9);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph3");
   graph->SetTitle("Res_68Quantile_PCA");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.6135428571);
   graph->SetPoint(1,150,0.3097916667);
   graph->SetPoint(2,250,0.2381730769);
   graph->SetPoint(3,350,0.202260274);
   graph->SetPoint(4,450,0.1831666667);
   graph->SetPoint(5,550,0.1641860465);
   graph->SetPoint(6,650,0.1483928571);
   graph->SetPoint(7,750,0.1428217822);
   graph->SetPoint(8,850,0.1315384615);
   graph->SetPoint(9,950,0.1230555556);
   graph->SetPoint(10,1500,0.1020075758);
   
   TH1F *Graph_Graph10 = new TH1F("Graph_Graph10","Res_68Quantile_PCA",100,0,1645);
   Graph_Graph10->SetMinimum(0.05085405);
   Graph_Graph10->SetMaximum(0.6646964);
   Graph_Graph10->SetDirectory(0);
   Graph_Graph10->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph10->SetLineColor(ci);
   Graph_Graph10->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph10->GetXaxis()->SetLabelFont(42);
   Graph_Graph10->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph10->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph10->GetXaxis()->SetTitleFont(42);
   Graph_Graph10->GetYaxis()->SetTitle("AngularResolution");
   Graph_Graph10->GetYaxis()->SetLabelFont(42);
   Graph_Graph10->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph10->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph10->GetYaxis()->SetTitleFont(42);
   Graph_Graph10->GetZaxis()->SetLabelFont(42);
   Graph_Graph10->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph10->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph10->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph10);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph4");
   graph->SetTitle("Res_68Quantile_PCA");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(25);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.5361904762);
   graph->SetPoint(1,150,0.2617);
   graph->SetPoint(2,250,0.2015384615);
   graph->SetPoint(3,350,0.1727848101);
   graph->SetPoint(4,450,0.1494148936);
   graph->SetPoint(5,550,0.136366);
   graph->SetPoint(6,650,0.1243939394);
   graph->SetPoint(7,750,0.1174806202);
   graph->SetPoint(8,850,0.1113076923);
   graph->SetPoint(9,950,0.1048728814);
   graph->SetPoint(10,1500,0.08234939759);
   
   TH1F *Graph_Graph11 = new TH1F("Graph_Graph11","Res_68Quantile_PCA",100,0,1645);
   Graph_Graph11->SetMinimum(0.03696529);
   Graph_Graph11->SetMaximum(0.5815746);
   Graph_Graph11->SetDirectory(0);
   Graph_Graph11->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph11->SetLineColor(ci);
   Graph_Graph11->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph11->GetXaxis()->SetLabelFont(42);
   Graph_Graph11->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph11->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph11->GetXaxis()->SetTitleFont(42);
   Graph_Graph11->GetYaxis()->SetTitle("AngularResolution");
   Graph_Graph11->GetYaxis()->SetLabelFont(42);
   Graph_Graph11->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph11->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph11->GetYaxis()->SetTitleFont(42);
   Graph_Graph11->GetZaxis()->SetLabelFont(42);
   Graph_Graph11->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph11->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph11->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph11);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph5");
   graph->SetTitle("Res_68Quantile_PCA");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.6192714286);
   graph->SetPoint(1,150,0.3061363636);
   graph->SetPoint(2,250,0.2333695652);
   graph->SetPoint(3,350,0.2008181818);
   graph->SetPoint(4,450,0.1751948052);
   graph->SetPoint(5,550,0.1588586957);
   graph->SetPoint(6,650,0.1492783505);
   graph->SetPoint(7,750,0.1371138211);
   graph->SetPoint(8,850,0.1298290598);
   graph->SetPoint(9,950,0.123287037);
   graph->SetPoint(10,1500,0.09968);
   
   TH1F *Graph_Graph12 = new TH1F("Graph_Graph12","Res_68Quantile_PCA",100,0,1645);
   Graph_Graph12->SetMinimum(0.04772086);
   Graph_Graph12->SetMaximum(0.6712306);
   Graph_Graph12->SetDirectory(0);
   Graph_Graph12->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph12->SetLineColor(ci);
   Graph_Graph12->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph12->GetXaxis()->SetLabelFont(42);
   Graph_Graph12->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetXaxis()->SetTitleFont(42);
   Graph_Graph12->GetYaxis()->SetTitle("AngularResolution");
   Graph_Graph12->GetYaxis()->SetLabelFont(42);
   Graph_Graph12->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetYaxis()->SetTitleFont(42);
   Graph_Graph12->GetZaxis()->SetLabelFont(42);
   Graph_Graph12->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph12);
   
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
