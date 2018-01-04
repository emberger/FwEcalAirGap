{
//=========Macro generated from canvas: GeometryRating/GeometryRating
//=========  (Fri Dec  8 16:44:02 2017) by ROOT version5.34/36
   TCanvas *GeometryRating = new TCanvas("GeometryRating", "GeometryRating",65,24,700,500);
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
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.07263130301);
   graph->SetPoint(1,150,0.0410901937);
   graph->SetPoint(2,250,0.02486509104);
   graph->SetPoint(3,350,-0.03295882048);
   graph->SetPoint(4,450,-0.03540709727);
   graph->SetPoint(5,550,0.1078793696);
   graph->SetPoint(6,650,-0.08674504561);
   graph->SetPoint(7,750,-0.02688279826);
   graph->SetPoint(8,850,0.03070120163);
   graph->SetPoint(9,950,-0.08176311077);
   graph->SetPoint(10,1500,0.060402513);
   
   TH1F *Graph_Graph109 = new TH1F("Graph_Graph109","Graph",100,0,1645);
   Graph_Graph109->SetMinimum(-0.2);
   Graph_Graph109->SetMaximum(0.8);
   Graph_Graph109->SetDirectory(0);
   Graph_Graph109->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph109->SetLineColor(ci);
   Graph_Graph109->GetXaxis()->SetLabelFont(42);
   Graph_Graph109->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph109->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph109->GetXaxis()->SetTitleFont(42);
   Graph_Graph109->GetYaxis()->SetLabelFont(42);
   Graph_Graph109->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph109->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph109->GetYaxis()->SetTitleFont(42);
   Graph_Graph109->GetZaxis()->SetLabelFont(42);
   Graph_Graph109->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph109->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph109->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph109);
   
   graph->Draw("ap");
   
   graph = new TGraph(11);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.4076822363);
   graph->SetPoint(1,150,0.1255719281);
   graph->SetPoint(2,250,0.14712645);
   graph->SetPoint(3,350,0.09414166067);
   graph->SetPoint(4,450,0.136837787);
   graph->SetPoint(5,550,0.2876271686);
   graph->SetPoint(6,650,0.09944932134);
   graph->SetPoint(7,750,0.1377081048);
   graph->SetPoint(8,850,0.1551388807);
   graph->SetPoint(9,950,0.04396324146);
   graph->SetPoint(10,1500,0.1534289062);
   
   TH1F *Graph_Graph110 = new TH1F("Graph_Graph110","Graph",100,0,1645);
   Graph_Graph110->SetMinimum(0.007591342);
   Graph_Graph110->SetMaximum(0.4440541);
   Graph_Graph110->SetDirectory(0);
   Graph_Graph110->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph110->SetLineColor(ci);
   Graph_Graph110->GetXaxis()->SetLabelFont(42);
   Graph_Graph110->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph110->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph110->GetXaxis()->SetTitleFont(42);
   Graph_Graph110->GetYaxis()->SetLabelFont(42);
   Graph_Graph110->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph110->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph110->GetYaxis()->SetTitleFont(42);
   Graph_Graph110->GetZaxis()->SetLabelFont(42);
   Graph_Graph110->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph110->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph110->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph110);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.08445983921);
   graph->SetPoint(1,150,0.01433415312);
   graph->SetPoint(2,250,0.03348904476);
   graph->SetPoint(3,350,-0.03770418966);
   graph->SetPoint(4,450,-0.03098837152);
   graph->SetPoint(5,550,0.1205270053);
   graph->SetPoint(6,650,0.02194542028);
   graph->SetPoint(7,750,0.03252311641);
   graph->SetPoint(8,850,0.02153624287);
   graph->SetPoint(9,950,-0.01879103588);
   graph->SetPoint(10,1500,0.09612224407);
   
   TH1F *Graph_Graph111 = new TH1F("Graph_Graph111","Graph",100,0,1645);
   Graph_Graph111->SetMinimum(-0.05352731);
   Graph_Graph111->SetMaximum(0.1363501);
   Graph_Graph111->SetDirectory(0);
   Graph_Graph111->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph111->SetLineColor(ci);
   Graph_Graph111->GetXaxis()->SetLabelFont(42);
   Graph_Graph111->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph111->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph111->GetXaxis()->SetTitleFont(42);
   Graph_Graph111->GetYaxis()->SetLabelFont(42);
   Graph_Graph111->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph111->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph111->GetYaxis()->SetTitleFont(42);
   Graph_Graph111->GetZaxis()->SetLabelFont(42);
   Graph_Graph111->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph111->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph111->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph111);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.7110146025);
   graph->SetPoint(1,150,0.1766575383);
   graph->SetPoint(2,250,0.2182918465);
   graph->SetPoint(3,350,0.06476050524);
   graph->SetPoint(4,450,0.1110469123);
   graph->SetPoint(5,550,0.2418876546);
   graph->SetPoint(6,650,0.1009274665);
   graph->SetPoint(7,750,0.2002563525);
   graph->SetPoint(8,850,0.1532271737);
   graph->SetPoint(9,950,0.1367333174);
   graph->SetPoint(10,1500,0.2222756886);
   
   TH1F *Graph_Graph112 = new TH1F("Graph_Graph112","Graph",100,0,1645);
   Graph_Graph112->SetMinimum(0.0001350955);
   Graph_Graph112->SetMaximum(0.77564);
   Graph_Graph112->SetDirectory(0);
   Graph_Graph112->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph112->SetLineColor(ci);
   Graph_Graph112->GetXaxis()->SetLabelFont(42);
   Graph_Graph112->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph112->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph112->GetXaxis()->SetTitleFont(42);
   Graph_Graph112->GetYaxis()->SetLabelFont(42);
   Graph_Graph112->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph112->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph112->GetYaxis()->SetTitleFont(42);
   Graph_Graph112->GetZaxis()->SetLabelFont(42);
   Graph_Graph112->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph112->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph112->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph112);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(25);
   graph->SetMarkerSize(1.5);
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
   
   TH1F *Graph_Graph113 = new TH1F("Graph_Graph113","Graph",100,0,1645);
   Graph_Graph113->SetMinimum(0);
   Graph_Graph113->SetMaximum(1.1);
   Graph_Graph113->SetDirectory(0);
   Graph_Graph113->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph113->SetLineColor(ci);
   Graph_Graph113->GetXaxis()->SetLabelFont(42);
   Graph_Graph113->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph113->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph113->GetXaxis()->SetTitleFont(42);
   Graph_Graph113->GetYaxis()->SetLabelFont(42);
   Graph_Graph113->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph113->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph113->GetYaxis()->SetTitleFont(42);
   Graph_Graph113->GetZaxis()->SetLabelFont(42);
   Graph_Graph113->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph113->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph113->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph113);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.711992906);
   graph->SetPoint(1,150,0.08469581096);
   graph->SetPoint(2,250,0.1786658107);
   graph->SetPoint(3,350,0.1112343826);
   graph->SetPoint(4,450,0.07163506927);
   graph->SetPoint(5,550,0.2270046647);
   graph->SetPoint(6,650,0.06345057399);
   graph->SetPoint(7,750,0.1319853877);
   graph->SetPoint(8,850,0.1627140253);
   graph->SetPoint(9,950,0.1343725084);
   graph->SetPoint(10,1500,0.2348716681);
   
   TH1F *Graph_Graph114 = new TH1F("Graph_Graph114","Graph",100,0,1645);
   Graph_Graph114->SetMinimum(0);
   Graph_Graph114->SetMaximum(0.7768471);
   Graph_Graph114->SetDirectory(0);
   Graph_Graph114->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph114->SetLineColor(ci);
   Graph_Graph114->GetXaxis()->SetLabelFont(42);
   Graph_Graph114->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph114->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph114->GetXaxis()->SetTitleFont(42);
   Graph_Graph114->GetYaxis()->SetLabelFont(42);
   Graph_Graph114->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph114->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph114->GetYaxis()->SetTitleFont(42);
   Graph_Graph114->GetZaxis()->SetLabelFont(42);
   Graph_Graph114->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph114->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph114->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph114);
   
   graph->Draw("p");
   
   TLegend *leg = new TLegend(0.5,0.7,0.9,0.9,NULL,"brNDC");
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

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","2mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","1mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph3","2mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph4","1mmLeadAbsorber_NOVessel_standard","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph5","2mmLeadAbsorber_NOVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.5);
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
