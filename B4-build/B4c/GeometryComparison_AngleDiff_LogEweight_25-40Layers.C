{
//=========Macro generated from canvas: Comparison/Comparison
//=========  (Tue Dec 19 10:22:24 2017) by ROOT version5.34/36
   TCanvas *Comparison = new TCanvas("Comparison", "Comparison",0,0,700,500);
   Comparison->Range(0,0,1,1);
   Comparison->SetFillColor(0);
   Comparison->SetBorderMode(0);
   Comparison->SetBorderSize(2);
   Comparison->SetFrameBorderMode(0);
   
   TGraph *graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("LogE_Res_68Quantile_PCA");
   graph->SetFillColor(1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.5277586207);
   graph->SetPoint(1,150,0.2676851852);
   graph->SetPoint(2,250,0.2163636364);
   graph->SetPoint(3,350,0.1788732394);
   graph->SetPoint(4,450,0.1573888889);
   graph->SetPoint(5,550,0.1432777778);
   graph->SetPoint(6,650,0.1306557377);
   graph->SetPoint(7,750,0.1236507937);
   graph->SetPoint(8,850,0.114789916);
   graph->SetPoint(9,950,0.1071875);
   graph->SetPoint(10,1500,0.08418478261);
   
   TH1F *Graph_Graph9 = new TH1F("Graph_Graph9","LogE_Res_68Quantile_PCA",100,0,1645);
   Graph_Graph9->SetMinimum(0);
   Graph_Graph9->SetMaximum(0.7);
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
   Graph_Graph9->GetYaxis()->SetTitleOffset(1.23);
   Graph_Graph9->GetYaxis()->SetTitleFont(42);
   Graph_Graph9->GetZaxis()->SetLabelFont(42);
   Graph_Graph9->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph9->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph9->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph9);
   
   graph->Draw("ap");
   
   graph = new TGraph(11);
   graph->SetName("Graph1");
   graph->SetTitle("LogE_Res_68Quantile_PCA");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.63106);
   graph->SetPoint(1,150,0.3211363636);
   graph->SetPoint(2,250,0.2509302326);
   graph->SetPoint(3,350,0.2114655172);
   graph->SetPoint(4,450,0.1842073171);
   graph->SetPoint(5,550,0.1712307692);
   graph->SetPoint(6,650,0.1568817204);
   graph->SetPoint(7,750,0.1420058252);
   graph->SetPoint(8,850,0.1380337079);
   graph->SetPoint(9,950,0.1293181818);
   graph->SetPoint(10,1500,0.1032733813);
   
   TH1F *Graph_Graph10 = new TH1F("Graph_Graph10","LogE_Res_68Quantile_PCA",100,0,1645);
   Graph_Graph10->SetMinimum(0.05049472);
   Graph_Graph10->SetMaximum(0.6838387);
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
   graph->SetName("Graph2");
   graph->SetTitle("LogE_Res_68Quantile_PCA");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.5540166667);
   graph->SetPoint(1,150,0.2761111111);
   graph->SetPoint(2,250,0.2163768116);
   graph->SetPoint(3,350,0.1800684932);
   graph->SetPoint(4,450,0.1619277108);
   graph->SetPoint(5,550,0.1405208333);
   graph->SetPoint(6,650,0.1350952381);
   graph->SetPoint(7,750,0.1239795918);
   graph->SetPoint(8,850,0.1159489051);
   graph->SetPoint(9,950,0.108358209);
   graph->SetPoint(10,1500,0.086625);
   
   TH1F *Graph_Graph11 = new TH1F("Graph_Graph11","LogE_Res_68Quantile_PCA",100,0,1645);
   Graph_Graph11->SetMinimum(0.03988583);
   Graph_Graph11->SetMaximum(0.6007558);
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
   graph->SetName("Graph3");
   graph->SetTitle("LogE_Res_68Quantile_PCA");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(25);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.623625);
   graph->SetPoint(1,150,0.31828125);
   graph->SetPoint(2,250,0.2483783784);
   graph->SetPoint(3,350,0.2105357143);
   graph->SetPoint(4,450,0.180890411);
   graph->SetPoint(5,550,0.1655797101);
   graph->SetPoint(6,650,0.1533888889);
   graph->SetPoint(7,750,0.1427222222);
   graph->SetPoint(8,850,0.1375520833);
   graph->SetPoint(9,950,0.1262820513);
   graph->SetPoint(10,1500,0.103556338);
   
   TH1F *Graph_Graph12 = new TH1F("Graph_Graph12","LogE_Res_68Quantile_PCA",100,0,1645);
   Graph_Graph12->SetMinimum(0.05154947);
   Graph_Graph12->SetMaximum(0.6756319);
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
   
   graph = new TGraph(11);
   graph->SetName("Graph4");
   graph->SetTitle("LogE_Res_68Quantile_PCA");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.5476333333);
   graph->SetPoint(1,150,0.2777906977);
   graph->SetPoint(2,250,0.2141176471);
   graph->SetPoint(3,350,0.1776666667);
   graph->SetPoint(4,450,0.1538709677);
   graph->SetPoint(5,550,0.1369186047);
   graph->SetPoint(6,650,0.1283027523);
   graph->SetPoint(7,750,0.117421875);
   graph->SetPoint(8,850,0.1107142857);
   graph->SetPoint(9,950,0.1050393701);
   graph->SetPoint(10,1500,0.08417586207);
   
   TH1F *Graph_Graph13 = new TH1F("Graph_Graph13","LogE_Res_68Quantile_PCA",100,0,1645);
   Graph_Graph13->SetMinimum(0.03783011);
   Graph_Graph13->SetMaximum(0.5939791);
   Graph_Graph13->SetDirectory(0);
   Graph_Graph13->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph13->SetLineColor(ci);
   Graph_Graph13->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph13->GetXaxis()->SetLabelFont(42);
   Graph_Graph13->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph13->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph13->GetXaxis()->SetTitleFont(42);
   Graph_Graph13->GetYaxis()->SetTitle("AngularResolution");
   Graph_Graph13->GetYaxis()->SetLabelFont(42);
   Graph_Graph13->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph13->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph13->GetYaxis()->SetTitleFont(42);
   Graph_Graph13->GetZaxis()->SetLabelFont(42);
   Graph_Graph13->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph13->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph13->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph13);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph5");
   graph->SetTitle("LogE_Res_68Quantile_PCA");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.6352705882);
   graph->SetPoint(1,150,0.3327380952);
   graph->SetPoint(2,250,0.2440833333);
   graph->SetPoint(3,350,0.2127118644);
   graph->SetPoint(4,450,0.1854477612);
   graph->SetPoint(5,550,0.1718421053);
   graph->SetPoint(6,650,0.1559770115);
   graph->SetPoint(7,750,0.1443298969);
   graph->SetPoint(8,850,0.1369565217);
   graph->SetPoint(9,950,0.1304310345);
   graph->SetPoint(10,1500,0.1014506173);
   
   TH1F *Graph_Graph14 = new TH1F("Graph_Graph14","LogE_Res_68Quantile_PCA",100,0,1645);
   Graph_Graph14->SetMinimum(0.04806862);
   Graph_Graph14->SetMaximum(0.6886526);
   Graph_Graph14->SetDirectory(0);
   Graph_Graph14->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph14->SetLineColor(ci);
   Graph_Graph14->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph14->GetXaxis()->SetLabelFont(42);
   Graph_Graph14->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph14->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph14->GetXaxis()->SetTitleFont(42);
   Graph_Graph14->GetYaxis()->SetTitle("AngularResolution");
   Graph_Graph14->GetYaxis()->SetLabelFont(42);
   Graph_Graph14->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph14->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph14->GetYaxis()->SetTitleFont(42);
   Graph_Graph14->GetZaxis()->SetLabelFont(42);
   Graph_Graph14->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph14->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph14->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph14);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph6");
   graph->SetTitle("LogE_Res_68Quantile_PCA");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.5746);
   graph->SetPoint(1,150,0.2820181818);
   graph->SetPoint(2,250,0.2185087719);
   graph->SetPoint(3,350,0.1851388889);
   graph->SetPoint(4,450,0.1548484848);
   graph->SetPoint(5,550,0.145);
   graph->SetPoint(6,650,0.13140625);
   graph->SetPoint(7,750,0.1243867925);
   graph->SetPoint(8,850,0.1143307087);
   graph->SetPoint(9,950,0.1060606061);
   graph->SetPoint(10,1500,0.085);
   
   TH1F *Graph_Graph15 = new TH1F("Graph_Graph15","LogE_Res_68Quantile_PCA",100,0,1645);
   Graph_Graph15->SetMinimum(0.03604);
   Graph_Graph15->SetMaximum(0.62356);
   Graph_Graph15->SetDirectory(0);
   Graph_Graph15->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph15->SetLineColor(ci);
   Graph_Graph15->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph15->GetXaxis()->SetLabelFont(42);
   Graph_Graph15->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph15->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph15->GetXaxis()->SetTitleFont(42);
   Graph_Graph15->GetYaxis()->SetTitle("AngularResolution");
   Graph_Graph15->GetYaxis()->SetLabelFont(42);
   Graph_Graph15->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph15->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph15->GetYaxis()->SetTitleFont(42);
   Graph_Graph15->GetZaxis()->SetLabelFont(42);
   Graph_Graph15->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph15->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph15->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph15);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph7");
   graph->SetTitle("LogE_Res_68Quantile_PCA");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(25);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.6259);
   graph->SetPoint(1,150,0.3231707317);
   graph->SetPoint(2,250,0.2496938776);
   graph->SetPoint(3,350,0.2108163265);
   graph->SetPoint(4,450,0.190915493);
   graph->SetPoint(5,550,0.1731506849);
   graph->SetPoint(6,650,0.155125);
   graph->SetPoint(7,750,0.1479411765);
   graph->SetPoint(8,850,0.13775);
   graph->SetPoint(9,950,0.1276470588);
   graph->SetPoint(10,1500,0.1047101449);
   
   TH1F *Graph_Graph16 = new TH1F("Graph_Graph16","LogE_Res_68Quantile_PCA",100,0,1645);
   Graph_Graph16->SetMinimum(0.05259116);
   Graph_Graph16->SetMaximum(0.678019);
   Graph_Graph16->SetDirectory(0);
   Graph_Graph16->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph16->SetLineColor(ci);
   Graph_Graph16->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph16->GetXaxis()->SetLabelFont(42);
   Graph_Graph16->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph16->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph16->GetXaxis()->SetTitleFont(42);
   Graph_Graph16->GetYaxis()->SetTitle("AngularResolution");
   Graph_Graph16->GetYaxis()->SetLabelFont(42);
   Graph_Graph16->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph16->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph16->GetYaxis()->SetTitleFont(42);
   Graph_Graph16->GetZaxis()->SetLabelFont(42);
   Graph_Graph16->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph16->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph16->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph16);
   
   graph->Draw("p");
   
   TLegend *leg = new TLegend(-2.353437e-185,-2.353437e-185,-2.353437e-185,-2.353437e-185,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph0","40innerLayers_1mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph1","40innerLayers_2mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph2","40innerLayers_1mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph3","40innerLayers_2mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph4","25innerLayers_1mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph5","25innerLayers_2mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph6","25innerLayers_1mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph7","25innerLayers_2mmLeadAbsorber_20mmSteelVessel","p");
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
