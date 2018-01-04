{
        TFile * f1= new TFile("~/Geant4/Data/ForwardEcalWithPressureVessel/AngularResolution/AngReso_IT20_OT20_40Inner_50Outer_lead1mm_Polystyrene5mm_14mmTitanVessel/EnergyResolution.root");
        TCanvas * c1=(TCanvas*)f1->Get("c2p");
        TGraph * g1=(TGraph*)c1->GetListOfPrimitives()->At(1);


        TFile * f2= new TFile("~/Geant4/Data/ForwardEcalWithPressureVessel/AngularResolution/AngReso_IT20_OT20_40Inner_50Outer_lead2mm_Polystyrene5mm_14mmTitanVessel/EnergyResolution.root");
        TCanvas * c2=(TCanvas*)f2->Get("c2p");
        TGraph * g2=(TGraph*)c2->GetListOfPrimitives()->At(1);

        TFile * f3= new TFile("~/Geant4/Data/ForwardEcalWithPressureVessel/AngularResolution/AngReso_IT20_OT20_40Inner_50Outer_lead1mm_Polystyrene5mm_20mmSteelVessel/EnergyResolution.root");
        TCanvas * c3=(TCanvas*)f3->Get("c2p");
        TGraph * g3=(TGraph*)c3->GetListOfPrimitives()->At(1);

        TFile * f4= new TFile("~/Geant4/Data/ForwardEcalWithPressureVessel/AngularResolution/AngReso_IT20_OT20_40Inner_50Outer_lead2mm_Polystyrene5mm_20mmSteelVessel/EnergyResolution.root");
        TCanvas * c4=(TCanvas*)f4->Get("c2p");
        TGraph * g4=(TGraph*)c4->GetListOfPrimitives()->At(1);
/*
        TFile * f5= new TFile("AngReso_IT20_OT20_40Inner_50Outer_lead2mm_Polystyrene5mm_20mmSteelVessel/EnergyResolution.root");
        TCanvas * c5=(TCanvas*)f5->Get("c2p");
        TGraph * g5=(TGraph*)c5->GetListOfPrimitives()->At(1);

        TFile * f6= new TFile("AngReso_IT20_OT20_noVessel_25Inner_50Outer_lead2mm_Polystyrene5mm/EnergyResolution.root");
        TCanvas * c6=(TCanvas*)f6->Get("c2p");
        TGraph * g6=(TGraph*)c6->GetListOfPrimitives()->At(1);
 */

        TFile * f5= new TFile("~/Geant4/Data/ForwardEcalWithPressureVessel/AngularResolution/AngReso_IT20_OT20_25Inner_50Outer_lead1mm_Polystyrene5mm_14mmTitanVessel/EnergyResolution.root");
        TCanvas * c5=(TCanvas*)f5->Get("c2p");
        TGraph * g5=(TGraph*)c5->GetListOfPrimitives()->At(1);


        TFile * f6= new TFile("~/Geant4/Data/ForwardEcalWithPressureVessel/AngularResolution/AngReso_IT20_OT20_25Inner_50Outer_lead2mm_Polystyrene5mm_14mmTitanVessel/EnergyResolution.root");
        TCanvas * c6=(TCanvas*)f6->Get("c2p");
        TGraph * g6=(TGraph*)c6->GetListOfPrimitives()->At(1);

        TFile * f7= new TFile("~/Geant4/Data/ForwardEcalWithPressureVessel/AngularResolution/AngReso_IT20_OT20_25Inner_50Outer_lead1mm_Polystyrene5mm_20mmSteelVessel/EnergyResolution.root");
        TCanvas * c7=(TCanvas*)f7->Get("c2p");
        TGraph * g7=(TGraph*)c7->GetListOfPrimitives()->At(1);

        TFile * f8= new TFile("~/Geant4/Data/ForwardEcalWithPressureVessel/AngularResolution/AngReso_IT20_OT20_25Inner_50Outer_lead2mm_Polystyrene5mm_20mmSteelVessel/EnergyResolution.root");
        TCanvas * c8=(TCanvas*)f8->Get("c2p");
        TGraph * g8=(TGraph*)c8->GetListOfPrimitives()->At(1);
//Standard is 1cm Lead, NO Vessel

        Int_t nofPoints=g1->GetN();
        std::cout<<nofPoints<<std::endl;

        Double_t energy[nofPoints];

        Double_t r1[nofPoints];
        Double_t r2[nofPoints];
        Double_t r3[nofPoints];
        Double_t r4[nofPoints];
        Double_t r5[nofPoints];
        Double_t r6[nofPoints];
        Double_t r7[nofPoints];
        Double_t r8[nofPoints];

        for(Int_t i=0; i < nofPoints; i++) {

                Double_t stax,stay;
                g7->GetPoint(i, stax, stay);

                energy[i]=stax;

                Double_t x1,x2,x3,x4,x5,x6,x7, x8,y1,y2,y3,y4,y5,y6, y7,y8;

                g1->GetPoint(i,x1, y1);
                r1[i]=(y1/stay)-1;

                g2->GetPoint(i,x2, y2);
                r2[i]=(y2/stay) -1;

                g3->GetPoint(i,x3, y3);
                r3[i]=(y3/stay) -1;

                g4->GetPoint(i,x4, y4);
                r4[i]=(y4/stay) -1;

                g5->GetPoint(i,x5, y5);
                r5[i]=(y5/stay) -1;

                g6->GetPoint(i,x6, y6);
                r6[i]=(y6/stay) -1;

                g7->GetPoint(i,x7, y7);
                r7[i]=(y7/stay) -1;

                g8->GetPoint(i,x8, y8);
                r8[i]=(y8/stay) -1;
        }

        TGraph * rat1=new TGraph(nofPoints,energy, r1);
        TGraph * rat2=new TGraph(nofPoints,energy, r2);
        TGraph * rat3=new TGraph(nofPoints,energy, r3);
        TGraph * rat4=new TGraph(nofPoints,energy, r4);
        TGraph * rat5=new TGraph(nofPoints,energy, r5);
        TGraph * rat6=new TGraph(nofPoints,energy, r6);
        TGraph * rat7=new TGraph(nofPoints,energy, r7);
        TGraph * rat8=new TGraph(nofPoints,energy, r8);

        TCanvas * ratio1=new TCanvas("GeometryRating", "GeometryRating");

        rat1->SetMarkerSize(1.2);
        rat1->SetMarkerStyle(20);
        rat1->SetMarkerColor(kRed);

        rat2->SetMarkerSize(1.2);
        rat2->SetMarkerStyle(24);
        rat2->SetMarkerColor(kRed);

        rat3->SetMarkerSize(1.2);
        rat3->SetMarkerStyle(21);
        rat3->SetMarkerColor(kRed);

        rat4->SetMarkerSize(1.2);
        rat4->SetMarkerStyle(25);
        rat4->SetMarkerColor(kRed);

        rat5->SetMarkerSize(1.2);
        rat5->SetMarkerStyle(20);
        rat5->SetMarkerColor(kBlue);

        rat6->SetMarkerSize(1.2);
        rat6->SetMarkerStyle(24);
        rat6->SetMarkerColor(kBlue);

        rat7->SetMarkerSize(1.2);
        rat7->SetMarkerStyle(21);
        rat7->SetMarkerColor(kBlue);

        rat8->SetMarkerSize(1.2);
        rat8->SetMarkerStyle(25);
        rat8->SetMarkerColor(kBlue);
        TLegend * l1= new TLegend(0.5,0.7,0.9,0.9);
        l1->AddEntry(rat1, "40innerLayer_1mmLeadAbsorber_14mmTitanVessel", "p");
        l1->AddEntry(rat2, "40innerLayer_2mmLeadAbsorber_14mmTitanVessel", "p");
        l1->AddEntry(rat3, "40innerLayer_1mmLeadAbsorber_20mmSteelVessel", "p");
        l1->AddEntry(rat4, "40innerLayer_2mmLeadAbsorber_20mmSteelVessel", "p");
        //l1->AddEntry(rat5, "1mmLeadAbsorber_NOVessel_standard", "p");
        //l1->AddEntry(rat6, "2mmLeadAbsorber_NOVessel", "p");
        l1->AddEntry(rat5, "25innerLayer_1mmLeadAbsorber_14mmTitanVessel", "p");
        l1->AddEntry(rat6, "25innerLayer_2mmLeadAbsorber_14mmTitanVessel", "p");
        l1->AddEntry(rat7, "25innerLayer_1mmLeadAbsorber_20mmSteelVessel_standard", "p");
        l1->AddEntry(rat8, "25innerLayer_2mmLeadAbsorber_20mmSteelVessel", "p");
        ratio1->cd(0);
        rat1->GetYaxis()->SetRangeUser(-0.5,0.8);

        rat1->Draw("AP");
        rat2->Draw("P");
        rat3->Draw("P");
        rat4->Draw("P");
        rat5->Draw("P");
        rat6->Draw("P");
        rat7->Draw("P");
        rat8->Draw("P");
        l1->Draw();

        ratio1->Print("GeometryRating_EnergyResolution_25-40Layers.pdf");
        ratio1->Print("GeometryRating_EnergyResolution_25-40Layers.C");


        g1->SetMarkerSize(1.2);
        g1->SetMarkerStyle(20);
        g1->SetMarkerColor(kRed);

        g2->SetMarkerSize(1.2);
        g2->SetMarkerStyle(24);
        g2->SetMarkerColor(kRed);

        g3->SetMarkerSize(1.2);
        g3->SetMarkerStyle(21);
        g3->SetMarkerColor(kRed);

        g4->SetMarkerSize(1.2);
        g4->SetMarkerStyle(25);
        g4->SetMarkerColor(kRed);

        g5->SetMarkerSize(1.2);
        g5->SetMarkerStyle(20);
        g5->SetMarkerColor(kBlue);

        g6->SetMarkerSize(1.2);
        g6->SetMarkerStyle(24);
        g6->SetMarkerColor(kBlue);

        g7->SetMarkerSize(1.2);
        g7->SetMarkerStyle(21);
        g7->SetMarkerColor(kBlue);

        g8->SetMarkerSize(1.2);
        g8->SetMarkerStyle(25);
        g8->SetMarkerColor(kBlue);

        TLegend * l2= new TLegend(0.5,0.7,0.9,0.9);
        l2->AddEntry(g1, "40innerLayers_1mmLeadAbsorber_14mmTitanVessel", "p");
        l2->AddEntry(g2, "40innerLayers_2mmLeadAbsorber_14mmTitanVessel", "p");
        l2->AddEntry(g3, "40innerLayers_1mmLeadAbsorber_20mmSteelVessel", "p");
        l2->AddEntry(g4, "40innerLayers_2mmLeadAbsorber_20mmSteelVessel", "p");
        //l2->AddEntry(g5, "1mmLeadAbsorber_NOVessel", "p");
        //l2->AddEntry(g6, "2mmLeadAbsorber_NOVessel", "p");
        l2->AddEntry(g5, "25innerLayers_1mmLeadAbsorber_14mmTitanVessel", "p");
        l2->AddEntry(g6, "25innerLayers_2mmLeadAbsorber_14mmTitanVessel", "p");
        l2->AddEntry(g7, "25innerLayers_1mmLeadAbsorber_20mmSteelVessel", "p");
        l2->AddEntry(g8, "25innerLayers_2mmLeadAbsorber_20mmSteelVessel", "p");

        TCanvas * comp1 = new TCanvas("Comparison", "Comparison");
        g1->GetYaxis()->SetRangeUser(0.0,0.4);
        g1->GetYaxis()->SetTitleOffset(1.23);
        g1->Draw("AP");
        g2->Draw("P");
        g3->Draw("P");
        g4->Draw("P");
        g5->Draw("P");
        g6->Draw("P");
        g7->Draw("P");
        g8->Draw("P");
        l2->Draw();

        comp1->Print("GeometryComparison_EnergyResolution_25-40Layers.C");
        comp1->Print("GeometryComparison_EnergyResolution_25-40Layers.pdf");
}
