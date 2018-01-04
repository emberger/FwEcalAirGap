{
//=========Macro generated from canvas: Comparison/Comparison
//=========  (Fri Dec  8 23:53:40 2017) by ROOT version5.34/36
   TCanvas *Comparison = new TCanvas("Comparison", "Comparison",0,0,700,500);
   Comparison->Range(0,0,1,1);
   Comparison->SetFillColor(0);
   Comparison->SetBorderMode(0);
   Comparison->SetBorderSize(2);
   Comparison->SetFrameBorderMode(0);
   
   TGraph *graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("Res_68Quantile_PCA");
   graph->SetFillColor(1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.5166);
   graph->SetPoint(1,150,0.2716666667);
   graph->SetPoint(2,250,0.2055084746);
   graph->SetPoint(3,350,0.1745762712);
   graph->SetPoint(4,450,0.1509375);
   graph->SetPoint(5,550,0.1405670103);
   graph->SetPoint(6,650,0.1253763441);
   graph->SetPoint(7,750,0.1174489796);
   graph->SetPoint(8,850,0.1108955224);
   graph->SetPoint(9,950,0.1062601626);
   graph->SetPoint(10,1500,0.08481481481);
   
   TH1F *Graph_Graph7 = new TH1F("Graph_Graph7","Res_68Quantile_PCA",100,0,1645);
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
   Graph_Graph7->GetYaxis()->SetTitle("AngularResolution");
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
   graph->SetTitle("Res_68Quantile_PCA");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.6048583333);
   graph->SetPoint(1,150,0.3145348837);
   graph->SetPoint(2,250,0.2346875);
   graph->SetPoint(3,350,0.1971052632);
   graph->SetPoint(4,450,0.1750675676);
   graph->SetPoint(5,550,0.1616860465);
   graph->SetPoint(6,650,0.1494623656);
   graph->SetPoint(7,750,0.1389634146);
   graph->SetPoint(8,850,0.1315714286);
   graph->SetPoint(9,950,0.1232113821);
   graph->SetPoint(10,1500,0.09951219512);
   
   TH1F *Graph_Graph8 = new TH1F("Graph_Graph8","Res_68Quantile_PCA",100,0,1645);
   Graph_Graph8->SetMinimum(0.04897758);
   Graph_Graph8->SetMaximum(0.6553929);
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
   graph->SetPoint(0,50,0.5338888889);
   graph->SetPoint(1,150,0.2663636364);
   graph->SetPoint(2,250,0.2040677966);
   graph->SetPoint(3,350,0.1704320988);
   graph->SetPoint(4,450,0.1532967033);
   graph->SetPoint(5,550,0.1407058824);
   graph->SetPoint(6,650,0.1277717391);
   graph->SetPoint(7,750,0.1175438596);
   graph->SetPoint(8,850,0.1133333333);
   graph->SetPoint(9,950,0.1065725806);
   graph->SetPoint(10,1500,0.0855);
   
   TH1F *Graph_Graph9 = new TH1F("Graph_Graph9","Res_68Quantile_PCA",100,0,1645);
   Graph_Graph9->SetMinimum(0.04066111);
   Graph_Graph9->SetMaximum(0.5787278);
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
   graph->SetPoint(0,50,0.6069);
   graph->SetPoint(1,150,0.3073684211);
   graph->SetPoint(2,250,0.2383333333);
   graph->SetPoint(3,350,0.2025);
   graph->SetPoint(4,450,0.1727777778);
   graph->SetPoint(5,550,0.1587179487);
   graph->SetPoint(6,650,0.1474468468);
   graph->SetPoint(7,750,0.1375274725);
   graph->SetPoint(8,850,0.1298095238);
   graph->SetPoint(9,950,0.123490566);
   graph->SetPoint(10,1500,0.09770588235);
   
   TH1F *Graph_Graph10 = new TH1F("Graph_Graph10","Res_68Quantile_PCA",100,0,1645);
   Graph_Graph10->SetMinimum(0.04678647);
   Graph_Graph10->SetMaximum(0.6578194);
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
