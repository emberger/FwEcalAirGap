{
//=========Macro generated from canvas: Comparison/Comparison
//=========  (Tue Dec 19 11:41:20 2017) by ROOT version5.34/36
   TCanvas *Comparison = new TCanvas("Comparison", "Comparison",10,10,700,500);
   Comparison->Range(0,0,1,1);
   Comparison->SetFillColor(0);
   Comparison->SetBorderMode(0);
   Comparison->SetBorderSize(2);
   Comparison->SetFrameBorderMode(0);
   
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
   graph->SetPoint(0,50,0.2156142377);
   graph->SetPoint(1,150,0.1282309883);
   graph->SetPoint(2,250,0.1022331226);
   graph->SetPoint(3,350,0.09004129242);
   graph->SetPoint(4,450,0.07993393437);
   graph->SetPoint(5,550,0.07479817437);
   graph->SetPoint(6,650,0.06752035046);
   graph->SetPoint(7,750,0.06595487808);
   graph->SetPoint(8,850,0.05980950654);
   graph->SetPoint(9,950,0.05936360204);
   graph->SetPoint(10,1500,0.05155677615);
   
   TH1F *Graph_Graph_Graph19 = new TH1F("Graph_Graph_Graph19","Graph",100,0,1645);
   Graph_Graph_Graph19->SetMinimum(0);
   Graph_Graph_Graph19->SetMaximum(0.4);
   Graph_Graph_Graph19->SetDirectory(0);
   Graph_Graph_Graph19->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph19->SetLineColor(ci);
   Graph_Graph_Graph19->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph_Graph19->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph19->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph19->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph19->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph19->GetYaxis()->SetTitle("#frac{#sigma}{E_{mean}}");
   Graph_Graph_Graph19->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph19->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph19->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph19->GetYaxis()->SetTitleOffset(1.23);
   Graph_Graph_Graph19->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph19->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph19->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph19->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph19->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph19);
   
   graph->Draw("ap");
   
   graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.3570168124);
   graph->SetPoint(1,150,0.1987930393);
   graph->SetPoint(2,250,0.150599381);
   graph->SetPoint(3,350,0.1296579711);
   graph->SetPoint(4,450,0.1171639609);
   graph->SetPoint(5,550,0.1022675208);
   graph->SetPoint(6,650,0.09747782995);
   graph->SetPoint(7,750,0.08861915935);
   graph->SetPoint(8,850,0.08270391011);
   graph->SetPoint(9,950,0.07888457636);
   graph->SetPoint(10,1500,0.06331038747);
   
   TH1F *Graph_Graph_Graph110 = new TH1F("Graph_Graph_Graph110","Graph",100,0,1645);
   Graph_Graph_Graph110->SetMinimum(0.03393974);
   Graph_Graph_Graph110->SetMaximum(0.3863875);
   Graph_Graph_Graph110->SetDirectory(0);
   Graph_Graph_Graph110->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph110->SetLineColor(ci);
   Graph_Graph_Graph110->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph_Graph110->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph110->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph110->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph110->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph110->GetYaxis()->SetTitle("#frac{#sigma}{E_{mean}}");
   Graph_Graph_Graph110->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph110->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph110->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph110->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_Graph110->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph110->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph110->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph110->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph110->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph110);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.2321048228);
   graph->SetPoint(1,150,0.1417331372);
   graph->SetPoint(2,250,0.115240941);
   graph->SetPoint(3,350,0.1040677147);
   graph->SetPoint(4,450,0.0957296973);
   graph->SetPoint(5,550,0.09025443918);
   graph->SetPoint(6,650,0.08736763697);
   graph->SetPoint(7,750,0.08694809047);
   graph->SetPoint(8,850,0.08397111494);
   graph->SetPoint(9,950,0.08019784442);
   graph->SetPoint(10,1500,0.07248794711);
   
   TH1F *Graph_Graph_Graph111 = new TH1F("Graph_Graph_Graph111","Graph",100,0,1645);
   Graph_Graph_Graph111->SetMinimum(0.05652626);
   Graph_Graph_Graph111->SetMaximum(0.2480665);
   Graph_Graph_Graph111->SetDirectory(0);
   Graph_Graph_Graph111->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph111->SetLineColor(ci);
   Graph_Graph_Graph111->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph_Graph111->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph111->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph111->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph111->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph111->GetYaxis()->SetTitle("#frac{#sigma}{E_{mean}}");
   Graph_Graph_Graph111->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph111->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph111->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph111->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_Graph111->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph111->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph111->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph111->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph111->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph111);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(25);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.3557889592);
   graph->SetPoint(1,150,0.1954026502);
   graph->SetPoint(2,250,0.1576174155);
   graph->SetPoint(3,350,0.130912676);
   graph->SetPoint(4,450,0.1169960316);
   graph->SetPoint(5,550,0.1059771301);
   graph->SetPoint(6,650,0.09609492662);
   graph->SetPoint(7,750,0.08997919584);
   graph->SetPoint(8,850,0.08373187603);
   graph->SetPoint(9,950,0.0773896661);
   graph->SetPoint(10,1500,0.0646360573);
   
   TH1F *Graph_Graph_Graph112 = new TH1F("Graph_Graph_Graph112","Graph",100,0,1645);
   Graph_Graph_Graph112->SetMinimum(0.03552077);
   Graph_Graph_Graph112->SetMaximum(0.3849042);
   Graph_Graph_Graph112->SetDirectory(0);
   Graph_Graph_Graph112->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph112->SetLineColor(ci);
   Graph_Graph_Graph112->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph_Graph112->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph112->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph112->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph112->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph112->GetYaxis()->SetTitle("#frac{#sigma}{E_{mean}}");
   Graph_Graph_Graph112->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph112->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph112->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph112->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_Graph112->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph112->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph112->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph112->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph112->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph112);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.2375595912);
   graph->SetPoint(1,150,0.1401728555);
   graph->SetPoint(2,250,0.1165044368);
   graph->SetPoint(3,350,0.1041727781);
   graph->SetPoint(4,450,0.08760649255);
   graph->SetPoint(5,550,0.08248109014);
   graph->SetPoint(6,650,0.07809947215);
   graph->SetPoint(7,750,0.07221330669);
   graph->SetPoint(8,850,0.07012402745);
   graph->SetPoint(9,950,0.06659629329);
   graph->SetPoint(10,1500,0.05696975621);
   
   TH1F *Graph_Graph_Graph113 = new TH1F("Graph_Graph_Graph113","Graph",100,0,1645);
   Graph_Graph_Graph113->SetMinimum(0.03891077);
   Graph_Graph_Graph113->SetMaximum(0.2556186);
   Graph_Graph_Graph113->SetDirectory(0);
   Graph_Graph_Graph113->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph113->SetLineColor(ci);
   Graph_Graph_Graph113->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph_Graph113->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph113->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph113->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph113->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph113->GetYaxis()->SetTitle("#frac{#sigma}{E_{mean}}");
   Graph_Graph_Graph113->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph113->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph113->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph113->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_Graph113->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph113->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph113->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph113->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph113->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph113);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(24);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.3511299763);
   graph->SetPoint(1,150,0.2011815185);
   graph->SetPoint(2,250,0.1584775573);
   graph->SetPoint(3,350,0.1311148878);
   graph->SetPoint(4,450,0.1176327342);
   graph->SetPoint(5,550,0.1090366637);
   graph->SetPoint(6,650,0.0981233433);
   graph->SetPoint(7,750,0.09245857623);
   graph->SetPoint(8,850,0.08792819738);
   graph->SetPoint(9,950,0.08150146676);
   graph->SetPoint(10,1500,0.06716242509);
   
   TH1F *Graph_Graph_Graph114 = new TH1F("Graph_Graph_Graph114","Graph",100,0,1645);
   Graph_Graph_Graph114->SetMinimum(0.03876567);
   Graph_Graph_Graph114->SetMaximum(0.3795267);
   Graph_Graph_Graph114->SetDirectory(0);
   Graph_Graph_Graph114->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph114->SetLineColor(ci);
   Graph_Graph_Graph114->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph_Graph114->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph114->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph114->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph114->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph114->GetYaxis()->SetTitle("#frac{#sigma}{E_{mean}}");
   Graph_Graph_Graph114->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph114->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph114->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph114->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_Graph114->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph114->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph114->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph114->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph114->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph114);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.278787278);
   graph->SetPoint(1,150,0.1900317286);
   graph->SetPoint(2,250,0.1583516682);
   graph->SetPoint(3,350,0.1433374437);
   graph->SetPoint(4,450,0.1326724526);
   graph->SetPoint(5,550,0.1239665471);
   graph->SetPoint(6,650,0.1183735139);
   graph->SetPoint(7,750,0.1112554759);
   graph->SetPoint(8,850,0.1081776279);
   graph->SetPoint(9,950,0.1022031003);
   graph->SetPoint(10,1500,0.0874537829);
   
   TH1F *Graph_Graph_Graph115 = new TH1F("Graph_Graph_Graph115","Graph",100,0,1645);
   Graph_Graph_Graph115->SetMinimum(0.06832043);
   Graph_Graph_Graph115->SetMaximum(0.2979206);
   Graph_Graph_Graph115->SetDirectory(0);
   Graph_Graph_Graph115->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph115->SetLineColor(ci);
   Graph_Graph_Graph115->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph_Graph115->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph115->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph115->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph115->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph115->GetYaxis()->SetTitle("#frac{#sigma}{E_{mean}}");
   Graph_Graph_Graph115->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph115->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph115->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph115->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_Graph115->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph115->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph115->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph115->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph115->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph115);
   
   graph->Draw("p");
   
   graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(25);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,50,0.3587093458);
   graph->SetPoint(1,150,0.2087985854);
   graph->SetPoint(2,250,0.1649090567);
   graph->SetPoint(3,350,0.1393516104);
   graph->SetPoint(4,450,0.1254280848);
   graph->SetPoint(5,550,0.1136416187);
   graph->SetPoint(6,650,0.1070313999);
   graph->SetPoint(7,750,0.1009827353);
   graph->SetPoint(8,850,0.09794461423);
   graph->SetPoint(9,950,0.09289982362);
   graph->SetPoint(10,1500,0.07782368548);
   
   TH1F *Graph_Graph_Graph116 = new TH1F("Graph_Graph_Graph116","Graph",100,0,1645);
   Graph_Graph_Graph116->SetMinimum(0.04973512);
   Graph_Graph_Graph116->SetMaximum(0.3867979);
   Graph_Graph_Graph116->SetDirectory(0);
   Graph_Graph_Graph116->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph116->SetLineColor(ci);
   Graph_Graph_Graph116->GetXaxis()->SetTitle("Energy[MeV]");
   Graph_Graph_Graph116->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph116->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph116->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph116->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph116->GetYaxis()->SetTitle("#frac{#sigma}{E_{mean}}");
   Graph_Graph_Graph116->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph116->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph116->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph116->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_Graph116->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph116->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph116->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph116->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph116->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph116);
   
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
   entry=leg->AddEntry("Graph0","40innerLayers_2mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph0","40innerLayers_1mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph0","40innerLayers_2mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph0","25innerLayers_1mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph0","25innerLayers_2mmLeadAbsorber_14mmTitanVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph0","25innerLayers_1mmLeadAbsorber_20mmSteelVessel","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph0","25innerLayers_2mmLeadAbsorber_20mmSteelVessel","p");
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
