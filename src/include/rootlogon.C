void rootlogon()
{
	gROOT->SetStyle("Default");
	gStyle->SetDrawBorder(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasBorderSize(0);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetPadColor(kWhite);
	gStyle->SetFrameFillColor(kWhite);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(0);
	gStyle->SetOptFit(0);
	gStyle->SetStatColor(kWhite);
	gStyle->SetStatBorderSize(1);
	gStyle->SetStatFont(42);
	gStyle->SetOptStat("iou");
	gStyle->SetHistFillColor(kWhite);
	gStyle->SetTextFont(42);
	gStyle->SetTitleFont(42,"xyz");//Helvetica
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetTitleW(0.7);
	gStyle->SetTitleH(0.075);
	gStyle->SetTitleX(0.15);
	gStyle->SetTitleY(1.00);
	gStyle->SetTitleTextColor(kBlue);
	gStyle->SetFuncColor(kRed);
	gStyle->SetFuncWidth(3);
	gStyle->SetTitleXSize(0.06);
	gStyle->SetTitleXOffset(1.2);
	gStyle->SetTitleYSize(0.06);
	gStyle->SetTitleYOffset(1.15);
	gStyle->SetLabelFont(42,"xyz");
	gStyle->SetLabelSize(0.04,"XYZ");
	gStyle->SetLabelOffset(0.015,"x");
	gStyle->SetLabelOffset(0.03,"y");
	//gStyle->SetPadGridX(1);
	//gStyle->SetPadGridY(1);
	gStyle->SetLegendBorderSize(1);
	//gStyle->SetLegendFont(42);
	gStyle->SetPadColor(kWhite);
	gStyle->SetPadBorderSize(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetPadTopMargin(0.08);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.10);
	gStyle->SetTickLength(-0.02,"xy");
	//gStyle->SetPadTickX(1);
	//gStyle->SetPadTickY(1);

	//Color Palette
	const Int_t NRGBs = 5;
	const Int_t NCont = 99;

	//black, blue, magenta,red, yellow, white
	Double_t stops[NRGBs] = { 0.00, 0.15,  0.3,  0.7,   1.00};
	Double_t red[NRGBs]   = { 0.00, 0.00, 1.00, 1.00,  1.00}; 
	Double_t green[NRGBs] = { 0.00, 0.24, 0.00, 0.89,  1.00};
	Double_t blue[NRGBs]  = { 0.00, 1.00, 0.00, 0.01,  1.00};
	
        Int_t NewPalette = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
        gStyle->SetNumberContours(NCont);

        gStyle->UseCurrentStyle();
        gROOT->ForceStyle();

        printf("\n-- Custom gStyle is configured\n\n");
        return;
}

