{
//=========Macro generated from canvas: GeometryRating/GeometryRating
//=========  (Tue Dec 19 10:22:24 2017) by ROOT version5.34/36
   TCanvas *GeometryRating = new TCanvas("GeometryRating", "GeometryRating",65,52,700,500);
   GeometryRating->Range(-205.625,-0.325,1850.625,0.925);
   GeometryRating->SetFillColor(0);
   GeometryRating->SetBorderMode(0);
   GeometryRating->SetBorderSize(2);
   GeometryRating->SetFrameBorderMode(0);
   GeometryRating->SetFrameBorderMode(0);
   
   TGraph *graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,-0.08151997792);
   graph->SetPoint(1,150,-0.05082295241);
   graph->SetPoint(2,250,-0.009817159958);
   graph->SetPoint(3,350,-0.03384296778);
   graph->SetPoint(4,450,0.01640574038);
   graph->SetPoint(5,550,-0.01187739464);
   graph->SetPoint(6,650,-0.005711389642);
   graph->SetPoint(7,750,-0.005917017294);
   graph->SetPoint(8,850,0.004016482626);
   graph->SetPoint(9,950,0.010625);
   graph->SetPoint(10,1500,-0.009590792839);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,0,1645);
   Graph_Graph1->SetMinimum(-0.2);
   Graph_Graph1->SetMaximum(0.8);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   graph->Draw("ap");
   
   graph = new TGraph(11);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.09825965889);
   graph->SetPoint(1,150,0.1387080137);
   graph->SetPoint(2,250,0.1483760141);
   graph->SetPoint(3,350,0.1421993429);
   graph->SetPoint(4,450,0.1895971553);
   graph->SetPoint(5,550,0.1809018568);
   graph->SetPoint(6,650,0.1938680271);
   graph->SetPoint(7,750,0.1416471351);
   graph->SetPoint(8,850,0.2073196211);
   graph->SetPoint(9,950,0.2192857143);
   graph->SetPoint(10,1500,0.2149809564);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0,1645);
   Graph_Graph2->SetMinimum(0.08615705);
   Graph_Graph2->SetMaximum(0.2313883);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph2->SetLineColor(ci);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
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
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,-0.03582202112);
   graph->SetPoint(1,150,-0.02094570878);
   graph->SetPoint(2,250,-0.00975686384);
   graph->SetPoint(3,350,-0.02738698373);
   graph->SetPoint(4,450,0.04571711504);
   graph->SetPoint(5,550,-0.0308908046);
   graph->SetPoint(6,650,0.02807315554);
   graph->SetPoint(7,750,-0.003273664415);
   graph->SetPoint(8,850,0.01415364662);
   graph->SetPoint(9,950,0.02166311301);
   graph->SetPoint(10,1500,0.01911764706);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph",100,0,1645);
   Graph_Graph3->SetMinimum(-0.04397593);
   Graph_Graph3->SetMaximum(0.05387103);
   Graph_Graph3->SetDirectory(0);
   Graph_Graph3->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph3->SetLineColor(ci);
   Graph_Graph3->GetXaxis()->SetLabelFont(42);
   Graph_Graph3->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetXaxis()->SetTitleFont(42);
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
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(25);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.08532022276);
   graph->SetPoint(1,150,0.12858415);
   graph->SetPoint(2,250,0.1366975165);
   graph->SetPoint(3,350,0.1371771514);
   graph->SetPoint(4,450,0.1681768222);
   graph->SetPoint(5,550,0.1419290355);
   graph->SetPoint(6,650,0.1672876206);
   graph->SetPoint(7,750,0.1474065647);
   graph->SetPoint(8,850,0.203107065);
   graph->SetPoint(9,950,0.1906593407);
   graph->SetPoint(10,1500,0.2183098592);
   
   TH1F *Graph_Graph4 = new TH1F("Graph_Graph4","Graph",100,0,1645);
   Graph_Graph4->SetMinimum(0.07202126);
   Graph_Graph4->SetMaximum(0.2316088);
   Graph_Graph4->SetDirectory(0);
   Graph_Graph4->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph4->SetLineColor(ci);
   Graph_Graph4->GetXaxis()->SetLabelFont(42);
   Graph_Graph4->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph4->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph4->GetXaxis()->SetTitleFont(42);
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
   
   graph = new TGraph(11);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,-0.04693119851);
   graph->SetPoint(1,150,-0.01499011204);
   graph->SetPoint(2,250,-0.02009587456);
   graph->SetPoint(3,350,-0.04036009002);
   graph->SetPoint(4,450,-0.006312732782);
   graph->SetPoint(5,550,-0.05573376103);
   graph->SetPoint(6,650,-0.02361758064);
   graph->SetPoint(7,750,-0.0559940273);
   graph->SetPoint(8,850,-0.03163124754);
   graph->SetPoint(9,950,-0.0096287964);
   graph->SetPoint(10,1500,-0.009695740365);
   
   TH1F *Graph_Graph5 = new TH1F("Graph_Graph5","Graph",100,0,1645);
   Graph_Graph5->SetMinimum(-0.06096216);
   Graph_Graph5->SetMaximum(-0.001344603);
   Graph_Graph5->SetDirectory(0);
   Graph_Graph5->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph5->SetLineColor(ci);
   Graph_Graph5->GetXaxis()->SetLabelFont(42);
   Graph_Graph5->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph5->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph5->GetXaxis()->SetTitleFont(42);
   Graph_Graph5->GetYaxis()->SetLabelFont(42);
   Graph_Graph5->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph5->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph5->GetYaxis()->SetTitleFont(42);
   Graph_Graph5->GetZaxis()->SetLabelFont(42);
   Graph_Graph5->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph5->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph5->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph5);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.1055875187);
   graph->SetPoint(1,150,0.1798462535);
   graph->SetPoint(2,250,0.1170413489);
   graph->SetPoint(3,350,0.1489313006);
   graph->SetPoint(4,450,0.1976078512);
   graph->SetPoint(5,550,0.1851179673);
   graph->SetPoint(6,650,0.1869832028);
   graph->SetPoint(7,750,0.1603313669);
   graph->SetPoint(8,850,0.1978979519);
   graph->SetPoint(9,950,0.2297783251);
   graph->SetPoint(10,1500,0.1935366739);
   
   TH1F *Graph_Graph6 = new TH1F("Graph_Graph6","Graph",100,0,1645);
   Graph_Graph6->SetMinimum(0.09316844);
   Graph_Graph6->SetMaximum(0.2421974);
   Graph_Graph6->SetDirectory(0);
   Graph_Graph6->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph6->SetLineColor(ci);
   Graph_Graph6->GetXaxis()->SetLabelFont(42);
   Graph_Graph6->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph6->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph6->GetXaxis()->SetTitleFont(42);
   Graph_Graph6->GetYaxis()->SetLabelFont(42);
   Graph_Graph6->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph6->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph6->GetYaxis()->SetTitleFont(42);
   Graph_Graph6->GetZaxis()->SetLabelFont(42);
   Graph_Graph6->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph6->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph6->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph6);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0);
   graph->SetPoint(1,150,0);
   graph->SetPoint(2,250,0);
   graph->SetPoint(3,350,0);
   graph->SetPoint(4,450,0);
   graph->SetPoint(5,550,0);
   graph->SetPoint(6,650,0);
   graph->SetPoint(7,750,0);
   graph->SetPoint(8,850,0);
   graph->SetPoint(9,950,0);
   graph->SetPoint(10,1500,0);
   
   TH1F *Graph_Graph7 = new TH1F("Graph_Graph7","Graph",100,0,1645);
   Graph_Graph7->SetMinimum(0);
   Graph_Graph7->SetMaximum(1.1);
   Graph_Graph7->SetDirectory(0);
   Graph_Graph7->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph7->SetLineColor(ci);
   Graph_Graph7->GetXaxis()->SetLabelFont(42);
   Graph_Graph7->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph7->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph7->GetXaxis()->SetTitleFont(42);
   Graph_Graph7->GetYaxis()->SetLabelFont(42);
   Graph_Graph7->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph7->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph7->GetYaxis()->SetTitleFont(42);
   Graph_Graph7->GetZaxis()->SetLabelFont(42);
   Graph_Graph7->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph7->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph7->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph7);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph7");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(25);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.08927949878);
   graph->SetPoint(1,150,0.1459216197);
   graph->SetPoint(2,250,0.142717866);
   graph->SetPoint(3,350,0.1386928365);
   graph->SetPoint(4,450,0.2329180563);
   graph->SetPoint(5,550,0.1941426547);
   graph->SetPoint(6,650,0.1804994055);
   graph->SetPoint(7,750,0.1893640277);
   graph->SetPoint(8,850,0.2048381543);
   graph->SetPoint(9,950,0.2035294118);
   graph->SetPoint(10,1500,0.231884058);
   
   TH1F *Graph_Graph8 = new TH1F("Graph_Graph8","Graph",100,0,1645);
   Graph_Graph8->SetMinimum(0.07491564);
   Graph_Graph8->SetMaximum(0.2472819);
   Graph_Graph8->SetDirectory(0);
   Graph_Graph8->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph8->SetLineColor(ci);
   Graph_Graph8->GetXaxis()->SetLabelFont(42);
   Graph_Graph8->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph8->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph8->GetXaxis()->SetTitleFont(42);
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
   
   TLegend *leg = new TLegend(0.5,0.7,0.9,0.9,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph0","40innerLayer_1mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","40innerLayer_2mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","40innerLayer_1mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph3","40innerLayer_2mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph4","25innerLayer_1mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph5","25innerLayer_2mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph6","25innerLayer_1mmLeadAbsorber_20mmSteelVessel_standard","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph7","25innerLayer_2mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   leg->Draw();
   
   TPaveText *pt = new TPaveText(0.4398567,0.9342405,0.5601433,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *text = pt->AddText("Graph");
   pt->Draw();
   GeometryRating->Modified();
   GeometryRating->cd();
   GeometryRating->SetSelected(GeometryRating);
}
