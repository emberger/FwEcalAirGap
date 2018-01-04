{
//=========Macro generated from canvas: GeometryRating/GeometryRating
//=========  (Tue Dec 19 11:41:20 2017) by ROOT version5.34/36
   TCanvas *GeometryRating = new TCanvas("GeometryRating", "GeometryRating",65,52,700,500);
   GeometryRating->Range(-205.625,-0.6625,1850.625,0.9625);
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
   graph->SetPoint(0,50,-0.2265994372);
   graph->SetPoint(1,150,-0.3252127461);
   graph->SetPoint(2,250,-0.3543918815);
   graph->SetPoint(3,350,-0.3718229507);
   graph->SetPoint(4,450,-0.397509183);
   graph->SetPoint(5,550,-0.3966261372);
   graph->SetPoint(6,650,-0.429599171);
   graph->SetPoint(7,750,-0.4071763428);
   graph->SetPoint(8,850,-0.4471175999);
   graph->SetPoint(9,950,-0.4191604571);
   graph->SetPoint(10,1500,-0.4104683131);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,0,1645);
   Graph_Graph1->SetMinimum(-0.5);
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
   graph->SetPoint(0,50,0.2806065437);
   graph->SetPoint(1,150,0.04610446262);
   graph->SetPoint(2,250,-0.04895614442);
   graph->SetPoint(3,350,-0.09543544428);
   graph->SetPoint(4,450,-0.1168930808);
   graph->SetPoint(5,550,-0.1750393696);
   graph->SetPoint(6,650,-0.1765233054);
   graph->SetPoint(7,750,-0.2034624934);
   graph->SetPoint(8,850,-0.2354804617);
   graph->SetPoint(9,950,-0.2281586747);
   graph->SetPoint(10,1500,-0.2760703383);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0,1645);
   Graph_Graph2->SetMinimum(-0.331738);
   Graph_Graph2->SetMaximum(0.3362742);
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
   graph->SetPoint(0,50,-0.1674482979);
   graph->SetPoint(1,150,-0.2541606696);
   graph->SetPoint(2,250,-0.2722467511);
   graph->SetPoint(3,350,-0.2739669968);
   graph->SetPoint(4,450,-0.2784508358);
   graph->SetPoint(5,550,-0.2719452037);
   graph->SetPoint(6,650,-0.2619325549);
   graph->SetPoint(7,750,-0.2184825979);
   graph->SetPoint(8,850,-0.22376635);
   graph->SetPoint(9,950,-0.2153090835);
   graph->SetPoint(10,1500,-0.1711285126);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph",100,0,1645);
   Graph_Graph3->SetMinimum(-0.2895511);
   Graph_Graph3->SetMaximum(-0.156348);
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
   graph->SetPoint(0,50,0.276202278);
   graph->SetPoint(1,150,0.02826328866);
   graph->SetPoint(2,250,-0.004636848426);
   graph->SetPoint(3,350,-0.0866819401);
   graph->SetPoint(4,450,-0.1181588244);
   graph->SetPoint(5,550,-0.1451150925);
   graph->SetPoint(6,650,-0.1882058456);
   graph->SetPoint(7,750,-0.191238048);
   graph->SetPoint(8,850,-0.2259778877);
   graph->SetPoint(9,950,-0.242785533);
   graph->SetPoint(10,1500,-0.2609118192);
   
   TH1F *Graph_Graph4 = new TH1F("Graph_Graph4","Graph",100,0,1645);
   Graph_Graph4->SetMinimum(-0.3146232);
   Graph_Graph4->SetMaximum(0.3299137);
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
   graph->SetPoint(0,50,-0.1478822387);
   graph->SetPoint(1,150,-0.2623713075);
   graph->SetPoint(2,250,-0.2642677017);
   graph->SetPoint(3,350,-0.2732340173);
   graph->SetPoint(4,450,-0.339678352);
   graph->SetPoint(5,550,-0.3346504191);
   graph->SetPoint(6,650,-0.3402284888);
   graph->SetPoint(7,750,-0.3509235736);
   graph->SetPoint(8,850,-0.3517695957);
   graph->SetPoint(9,950,-0.3483926308);
   graph->SetPoint(10,1500,-0.3485729911);
   
   TH1F *Graph_Graph5 = new TH1F("Graph_Graph5","Graph",100,0,1645);
   Graph_Graph5->SetMinimum(-0.3721583);
   Graph_Graph5->SetMaximum(-0.1274935);
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
   graph->SetPoint(0,50,0.2594906729);
   graph->SetPoint(1,150,0.05867330675);
   graph->SetPoint(2,250,0.0007949970366);
   graph->SetPoint(3,350,-0.08527120018);
   graph->SetPoint(4,450,-0.1133597681);
   graph->SetPoint(5,550,-0.1204347761);
   graph->SetPoint(6,650,-0.1710701147);
   graph->SetPoint(7,750,-0.1689525796);
   graph->SetPoint(8,850,-0.1871868601);
   graph->SetPoint(9,950,-0.2025538702);
   graph->SetPoint(10,1500,-0.2320237861);
   
   TH1F *Graph_Graph6 = new TH1F("Graph_Graph6","Graph",100,0,1645);
   Graph_Graph6->SetMinimum(-0.2811752);
   Graph_Graph6->SetMaximum(0.3086421);
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
   graph->SetPoint(0,50,0.2866776004);
   graph->SetPoint(1,150,0.09875643903);
   graph->SetPoint(2,250,0.04141029);
   graph->SetPoint(3,350,-0.02780734157);
   graph->SetPoint(4,450,-0.05460340648);
   graph->SetPoint(5,550,-0.08328802094);
   graph->SetPoint(6,650,-0.09581631573);
   graph->SetPoint(7,750,-0.09233469631);
   graph->SetPoint(8,850,-0.09459454646);
   graph->SetPoint(9,950,-0.09102734283);
   graph->SetPoint(10,1500,-0.1101164192);
   
   TH1F *Graph_Graph8 = new TH1F("Graph_Graph8","Graph",100,0,1645);
   Graph_Graph8->SetMinimum(-0.1497958);
   Graph_Graph8->SetMaximum(0.326357);
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
