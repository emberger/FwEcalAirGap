{
//=========Macro generated from canvas: GeometryRating/GeometryRating
//=========  (Fri Dec  8 23:47:29 2017) by ROOT version5.34/36
   TCanvas *GeometryRating = new TCanvas("GeometryRating", "GeometryRating",0,0,700,500);
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
   graph->SetPoint(0,50,-0.01647284285);
   graph->SetPoint(1,150,0.0306542691);
   graph->SetPoint(2,250,-0.003167938931);
   graph->SetPoint(3,350,-0.005558855559);
   graph->SetPoint(4,450,0.05982512655);
   graph->SetPoint(5,550,-0.002131991719);
   graph->SetPoint(6,650,0.03726215957);
   graph->SetPoint(7,750,0.02144506763);
   graph->SetPoint(8,850,0.0001535744452);
   graph->SetPoint(9,950,0.02468376068);
   graph->SetPoint(10,1500,0.01532146631);
   
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
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,-0.1829439465);
   graph->SetPoint(1,150,-0.1818005081);
   graph->SetPoint(2,250,-0.1857506361);
   graph->SetPoint(3,350,-0.2007674865);
   graph->SetPoint(4,450,-0.2144959772);
   graph->SetPoint(5,550,-0.210332315);
   graph->SetPoint(6,650,-0.2060170524);
   graph->SetPoint(7,750,-0.2106410838);
   graph->SetPoint(8,850,-0.2144594228);
   graph->SetPoint(9,950,-0.2232027279);
   graph->SetPoint(10,1500,-0.1930724214);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0,1645);
   Graph_Graph2->SetMinimum(-0.2273429);
   Graph_Graph2->SetMaximum(-0.1776603);
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
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.0274178508);
   graph->SetPoint(1,150,0.04472429433);
   graph->SetPoint(2,250,0.03743638677);
   graph->SetPoint(3,350,0.01811897023);
   graph->SetPoint(4,450,0.04977387498);
   graph->SetPoint(5,550,0.03872709344);
   graph->SetPoint(6,650,0.0384160394);
   graph->SetPoint(7,750,0.02573685363);
   graph->SetPoint(8,850,0.005629568349);
   graph->SetPoint(9,950,0.02143862871);
   graph->SetPoint(10,1500,0.0321872714);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph",100,0,1645);
   Graph_Graph3->SetMinimum(0.001215138);
   Graph_Graph3->SetMaximum(0.05418831);
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
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.1769863825);
   graph->SetPoint(1,150,0.1922957949);
   graph->SetPoint(2,250,0.1908396947);
   graph->SetPoint(3,350,0.1392561285);
   graph->SetPoint(4,450,0.2031305625);
   graph->SetPoint(5,550,0.1879794084);
   graph->SetPoint(6,650,0.2172703378);
   graph->SetPoint(7,750,0.1996659518);
   graph->SetPoint(8,850,0.1852107809);
   graph->SetPoint(9,950,0.1968855219);
   graph->SetPoint(10,1500,0.2206462155);
   
   TH1F *Graph_Graph4 = new TH1F("Graph_Graph4","Graph",100,0,1645);
   Graph_Graph4->SetMinimum(0.1311171);
   Graph_Graph4->SetMaximum(0.2287852);
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
   
   TH1F *Graph_Graph5 = new TH1F("Graph_Graph5","Graph",100,0,1645);
   Graph_Graph5->SetMinimum(0);
   Graph_Graph5->SetMaximum(1.1);
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
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,50,0.154946714);
   graph->SetPoint(1,150,0.1697988675);
   graph->SetPoint(2,250,0.1579405908);
   graph->SetPoint(3,350,0.1622444222);
   graph->SetPoint(4,450,0.1725391021);
   graph->SetPoint(5,550,0.1649435758);
   graph->SetPoint(6,650,0.2000452051);
   graph->SetPoint(7,750,0.1671186359);
   graph->SetPoint(8,850,0.1663979114);
   graph->SetPoint(9,950,0.1755854845);
   graph->SetPoint(10,1500,0.2104520849);
   
   TH1F *Graph_Graph6 = new TH1F("Graph_Graph6","Graph",100,0,1645);
   Graph_Graph6->SetMinimum(0.1493962);
   Graph_Graph6->SetMaximum(0.2160026);
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
   
   TPaveText *pt = new TPaveText(0.4397126,0.9339831,0.5602874,0.995,"blNDC");
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
