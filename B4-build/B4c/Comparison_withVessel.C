{
//=========Macro generated from canvas: Comparison/Comparison
//=========  (Thu Dec  7 15:02:45 2017) by ROOT version5.34/36
   TCanvas *Comparison = new TCanvas("Comparison", "Comparison",65,24,1855,1176);
   Comparison->Range(-205.625,0.02264458,1850.625,0.3982118);
   Comparison->SetFillColor(0);
   Comparison->SetBorderMode(0);
   Comparison->SetBorderSize(2);
   Comparison->SetFrameBorderMode(0);
   Comparison->SetFrameBorderMode(0);
   
   TGraph *graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("GeometryComparison");
   graph->SetFillColor(1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#3333ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.3356172937);
   graph->SetPoint(1,150,0.2394236189);
   graph->SetPoint(2,250,0.2029440568);
   graph->SetPoint(3,350,0.1759398437);
   graph->SetPoint(4,450,0.1582281385);
   graph->SetPoint(5,550,0.1416048169);
   graph->SetPoint(6,650,0.1310871083);
   graph->SetPoint(7,750,0.1202862244);
   graph->SetPoint(8,850,0.1162779217);
   graph->SetPoint(9,950,0.1071979973);
   graph->SetPoint(10,1500,0.08523913105);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","GeometryComparison",100,0,1645);
   Graph_Graph1->SetMinimum(0.06020131);
   Graph_Graph1->SetMaximum(0.3606551);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   graph->Draw("ap");
   
   TPaveText *pt = new TPaveText(0.4470201,0.94,0.5463509,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *text = pt->AddText("ResX");
   pt->Draw();
   
   graph = new TGraph(11);
   graph->SetName("Graph1");
   graph->SetTitle("GeometryComparison");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#3333ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.3545808981);
   graph->SetPoint(1,150,0.2662277362);
   graph->SetPoint(2,250,0.2203233063);
   graph->SetPoint(3,350,0.1983167104);
   graph->SetPoint(4,450,0.1847892336);
   graph->SetPoint(5,550,0.1688906674);
   graph->SetPoint(6,650,0.1579251342);
   graph->SetPoint(7,750,0.1474929141);
   graph->SetPoint(8,850,0.1358279866);
   graph->SetPoint(9,950,0.1311398154);
   graph->SetPoint(10,1500,0.1049392643);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","ResX",100,0,1645);
   Graph_Graph2->SetMinimum(0.0799751);
   Graph_Graph2->SetMaximum(0.3795451);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph2->SetLineColor(ci);
   Graph_Graph2->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph2);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph2");
   graph->SetTitle("ResX");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#666666");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.3417256553);
   graph->SetPoint(1,150,0.2483132627);
   graph->SetPoint(2,250,0.2030209352);
   graph->SetPoint(3,350,0.1783962762);
   graph->SetPoint(4,450,0.1595869423);
   graph->SetPoint(5,550,0.1461653918);
   graph->SetPoint(6,650,0.1372038789);
   graph->SetPoint(7,750,0.125975252);
   graph->SetPoint(8,850,0.1173133618);
   graph->SetPoint(9,950,0.1089496742);
   graph->SetPoint(10,1500,0.09049972527);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","ResX",100,0,1645);
   Graph_Graph3->SetMinimum(0.06537713);
   Graph_Graph3->SetMaximum(0.3668482);
   Graph_Graph3->SetDirectory(0);
   Graph_Graph3->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph3->SetLineColor(ci);
   Graph_Graph3->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph3->GetXaxis()->SetLabelFont(42);
   Graph_Graph3->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetXaxis()->SetTitleFont(42);
   Graph_Graph3->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");
   Graph_Graph3->GetYaxis()->SetLabelFont(42);
   Graph_Graph3->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetYaxis()->SetTitleFont(42);
   Graph_Graph3->GetZaxis()->SetLabelFont(42);
   Graph_Graph3->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph3);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph3");
   graph->SetTitle("ResX");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#666666");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.3549209018);
   graph->SetPoint(1,150,0.2646672565);
   graph->SetPoint(2,250,0.2289022186);
   graph->SetPoint(3,350,0.2039783351);
   graph->SetPoint(4,450,0.1828218528);
   graph->SetPoint(5,550,0.1702202749);
   graph->SetPoint(6,650,0.160375872);
   graph->SetPoint(7,750,0.1478361454);
   graph->SetPoint(8,850,0.1419928172);
   graph->SetPoint(9,950,0.1320501768);
   graph->SetPoint(10,1500,0.1123572123);
   
   TH1F *Graph_Graph4 = new TH1F("Graph_Graph4","ResX",100,0,1645);
   Graph_Graph4->SetMinimum(0.08810084);
   Graph_Graph4->SetMaximum(0.3791773);
   Graph_Graph4->SetDirectory(0);
   Graph_Graph4->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph4->SetLineColor(ci);
   Graph_Graph4->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph4->GetXaxis()->SetLabelFont(42);
   Graph_Graph4->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph4->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph4->GetXaxis()->SetTitleFont(42);
   Graph_Graph4->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");
   Graph_Graph4->GetYaxis()->SetLabelFont(42);
   Graph_Graph4->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph4->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph4->GetYaxis()->SetTitleFont(42);
   Graph_Graph4->GetZaxis()->SetLabelFont(42);
   Graph_Graph4->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph4->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph4->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph4);
   
   graph->Draw("p");
   
   TLegend *leg = new TLegend(0.5196978,0.7,0.8996222,0.9,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph0","1mmleadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#3333ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","2mmleadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#3333ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","1mmleadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#666666");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph3","2mmleadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#666666");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   leg->Draw();
   Comparison->Modified();
   Comparison->cd();
   Comparison->SetSelected(Comparison);
}
