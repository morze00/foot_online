#include "ext_data_clnt.hh"
#include "ext_data_struct_info.hh"
#include "ext_data_client.h"
#include "include/ext_h101_unpack.h"
#include "include/ext_h101_foot.h"
#include "include/libs.hh"
#define LENGTH(x) (sizeof x / sizeof *x)
#define NDETS LENGTH(foot_id)
using namespace std;
typedef struct EXT_STR_h101_t
{
	EXT_STR_h101_unpack_t unpack;
	EXT_STR_h101_FOOT_onion foot;
	//EXT_STR_h101_FOOT foot;
} EXT_STR_h101;

namespace{//don't change this!!!
	int foot_id[] = {2, 13, 4, 11, 6, 7, 12, 9};
}

Double_t pedestal[16][640];
Double_t sigma[16][640];
Double_t fine_sigma[16][640];
auto pedfilename = "/u/land/r3broot/202402_s091_s118/R3BParams_s091_s118/macros/exp/online/foot/event_display/pedestal.dat";

//Clustering data structures
struct cluster_data
{
	double position;
	double energy;
	int mul;
};
typedef std::pair<int, double> strip_data;
typedef std::vector<strip_data> cluster; 
typedef std::vector<cluster> foot_data; 
typedef std::vector<cluster_data> hits;

const double NSIGMA = 5;
const double FOOT_LENGTH   = 96.;//mm

double get_cog(cluster c)//calculate center of gravity of a cluster
{
	double value = 0;  double esum = 0.;
	for(auto & s: c){
		value  += (s.first * s.second);
		esum += s.second;
	}
	value /= esum;//center of gravity
	return value;
}

double get_esum(cluster c)//calculate cluster sum
{
	double esum = 0.;
	for(auto & s: c){
		esum += s.second;
	}
	return esum;
}

void Check_Strip(UInt_t number, double energy, foot_data& fdata)
{
	//cout << "\n\n---------------------------------------";
	//cout << "\nChecking strip: " << number << "\t Energy: " << energy << endl;
	strip_data strip = std::make_pair(number,energy);
	cluster clust;
	if(fdata.size()==0)//no cluster yet, create new
	{
		clust.push_back(strip);
		fdata.push_back(clust);
		//cout << "\n\t New cluster is created for this strip";
		return;
	}
	cluster    this_clust = fdata.back();
	strip_data this_strip = this_clust.back();
	if(abs(strip.first-this_strip.first)<2)//neighbour found 
	{
		//cout << "\n\tStrip belong to exisitng cluster! Adding it...";
		fdata.back().push_back(strip);
		return;
	}
	else
	{
		//cout << "\n\tStrip is a new cluster! Making it...";
		clust.clear();
		clust.push_back(strip);
		fdata.push_back(clust);
		return;
	}
}

bool Read_Pedestals()
{
	ifstream pedfile(pedfilename,ifstream::in);
	if ( !pedfile.is_open() ){
		cout << "\n\nERROR: Cannot open pedestal file! Try to run first with the option --pedestal \n\n";
		return false;
	}
	std::cout << "\n-- Reading pedestals from the file " << pedfilename;
	int Nlines = 0;
	int strip_no, det_no;
	std::string line;
	while( std::getline(pedfile, line) ) Nlines++;
	cout << "\n-- Number of lines in the pedestal file: " << Nlines << "\n\n";
	pedfile.clear();
	pedfile.seekg( 0, std::ios::beg );
	for(int d=0; d<Nlines/641; d++){
		pedfile >>det_no;
		for(Int_t i=0 ; i<640 ; i++){
			pedfile >> strip_no >>  pedestal[det_no][i] >> sigma[det_no][i];
		}
	}
	pedfile.close();
	return true;
}

bool is_asic_edge_strip(UInt_t strip)
{
	if((strip%64) == 1)
		return true;//ASIC edge strip with small sigma
	else 
		return false;
}

bool is_good_strip(UInt_t det, UInt_t strip)
{
	switch(det){
		case 4:
			if(strip==205 || strip==204 || strip==203) return false;
		case 9:
			if(strip==20 || strip == 229 || strip==230) return false;
		case 12:
			if(strip==31 || strip == 322) return false;
	}
	if(is_asic_edge_strip(strip)) return false;//ASIC edge
	if((strip%64)>61 || (strip%64)<4) return false;//ASIC edge
	return true;
}
void Draw_Asic_Borders(double ymin, double ymax)
{
	for(int i_asic=1; i_asic<10; i_asic++)
	{
		TLine* l = new TLine(64*i_asic,ymin,64*i_asic,ymax);
		l->SetLineStyle(3);
		l->SetLineWidth(1);
		l->SetLineColor(14);
		l->Draw("same");
	}
}
class FootOnline : public TNamed
{
	public:
		FootOnline();
		void Clear(Option_t *option);
		TH2F Get_Clean_TH2F(TH2F* h2_in, int thr);
		void Calculate_Fine_Sigmas(int threshold);
		void Calculate_Pedestals(int threshold);
		void Draw_Pedestals(Bool_t is_fine_sigma);
		TCanvas * c_ped;
		TCanvas * c_fine;
		TCanvas * c_raw;
		TCanvas * c_calib;
		TCanvas * c_corr_xy;
		static constexpr size_t ARRAY_SIZE = 16;
		TH2F * h2_raw[ARRAY_SIZE] = {nullptr};
		TH2F * h2_calib[ARRAY_SIZE] = {nullptr};
		TH2F * h2_fine[ARRAY_SIZE] = {nullptr};
		TH1D * h1_peds[ARRAY_SIZE] = {nullptr};
		TH1D * h1_sigma[ARRAY_SIZE] = {nullptr};
		TH1D * h1_sigma_fine[ARRAY_SIZE] = {nullptr};

		TH2F * h2_left_X_corr = nullptr;
		TH2F * h2_left_Y_corr = nullptr;
		TH2F * h2_right_X_corr = nullptr;
		TH2F * h2_right_Y_corr = nullptr;

		bool is_fine_sigma;
};

void FootOnline::Clear(Option_t *option ="")
{
	for (int i=0; i<16; ++i){
		h2_raw[i]->Reset();
		h2_calib[i]->Reset();
		h2_fine[i]->Reset();
		if(h1_peds[i]=nullptr) h1_peds[i]->Reset();
		if(h1_sigma[i]!=nullptr) h1_sigma[i]->Reset();
		if(h1_sigma_fine[i]!=nullptr) h1_sigma_fine[i]->Reset();
		if(h2_left_X_corr!=nullptr) h2_left_X_corr->Reset();
		if(h2_left_Y_corr!=nullptr) h2_left_Y_corr->Reset();
		if(h2_right_X_corr!=nullptr) h2_right_X_corr->Reset();
		if(h2_right_Y_corr!=nullptr) h2_right_Y_corr->Reset();
	}
	gStyle->SetPalette(kRainBow);
}

FootOnline::FootOnline() : TNamed("FootAna","FootAna")
{
	c_ped = new TCanvas("Pedestals","Pedestals",1000, 1000);
	c_ped->Divide(4,4);
	c_raw = new TCanvas("Raw signals","Raw signals",1000,500);
	c_raw->Divide(4,2);
	c_calib = new TCanvas("Subtracted pedestals","Subtracted pedestals",1000,500);
	c_calib->Divide(4,2);
	c_fine = new TCanvas("Fine baseline","Fine baseline",1000,500);
	c_fine->Divide(4,2);
	c_corr_xy = new TCanvas("Position correlations","Position correlations", 1000, 1000);
	c_corr_xy->Divide(2,2);

	for(int i=0; i<16; ++i){
		TString hname = TString::Format("FOOT%d raw",i+1);
		h2_raw[i] = new TH2F(hname.Data(),hname.Data(),640,1,641,700,0,700);
		hname = TString::Format("FOOT%d calib",i+1);
		h2_calib[i] = new TH2F(hname.Data(),hname.Data(),640,1,641,200,-100,100);
		hname = TString::Format("FOOT%d calib fine",i+1);
		h2_fine[i] = new TH2F(hname.Data(),hname.Data(),640,1,641,200,-100,100);
		h2_raw[i]->SetDirectory(0);//do not store histo for http server
		h2_calib[i]->SetDirectory(0);
		h2_fine[i]->SetDirectory(0);
	}
	for(int i=0;i<NDETS;++i)
	{
		int d = foot_id[i]-1;
		c_raw->cd(i+1);
		h2_raw[d]->Draw("colz");
		Draw_Asic_Borders(h2_raw[d]->GetYaxis()->GetXmin(),
				h2_raw[d]->GetYaxis()->GetXmax());
		c_calib->cd(i+1);
		h2_calib[d]->Draw("colz");
		Draw_Asic_Borders(h2_calib[d]->GetYaxis()->GetXmin(),
				h2_calib[d]->GetYaxis()->GetXmax());
		c_fine->cd(i+1);
		h2_fine[d]->Draw("colz");
		Draw_Asic_Borders(h2_fine[d]->GetYaxis()->GetXmin(),
				h2_fine[d]->GetYaxis()->GetXmax());
	}

	int NBINS=400;
	double min = -50;
	double max = 50;

	c_corr_xy->cd(1);
	h2_left_X_corr = new TH2F("X_corr_left","X corr (left arm)",NBINS,min,max,NBINS,min,max);
	h2_left_X_corr->Draw("colz");
	h2_left_X_corr->GetXaxis()->SetTitle("Left FOOT 1, X [mm]");
	h2_left_X_corr->GetYaxis()->SetTitle("Left FOOT 3, X [mm]");
	h2_left_X_corr->SetDirectory(0);//do not store histo for http server

	c_corr_xy->cd(2);
	h2_left_Y_corr = new TH2F("Y_corr_left","Y corr (left arm)",NBINS,min,max,NBINS,min,max);
	h2_left_Y_corr->Draw("colz");
	h2_left_Y_corr->GetXaxis()->SetTitle("Left FOOT 2, Y [mm]");
	h2_left_Y_corr->GetYaxis()->SetTitle("Left FOOT 4, Y [mm]");
	h2_left_Y_corr->SetDirectory(0);//do not store histo for http server

	c_corr_xy->cd(3);
	h2_right_X_corr = new TH2F("X_corr_right","X corr (right arm)",NBINS,min,max,NBINS,min,max);
	h2_right_X_corr->Draw("colz");
	h2_right_X_corr->GetXaxis()->SetTitle("Right FOOT 1, Y [mm]");
	h2_right_X_corr->GetYaxis()->SetTitle("Right FOOT 3, Y [mm]");
	h2_right_X_corr->SetDirectory(0);//do not store histo for http server

	c_corr_xy->cd(4);
	h2_right_Y_corr = new TH2F("Y_corr_right","Y corr (right arm)",NBINS,min,max,NBINS,min,max);
	h2_right_Y_corr->Draw("colz");
	h2_right_Y_corr->GetXaxis()->SetTitle("Right FOOT 2, Y [mm]");
	h2_right_Y_corr->GetYaxis()->SetTitle("Right FOOT 4, Y [mm]");
	h2_right_Y_corr->SetDirectory(0);//do not store histo for http server

	is_fine_sigma = false;
}

TH2F FootOnline::Get_Clean_TH2F(TH2F* h2_in, int thr)
{
	TString hname = TString(h2_in->GetTitle()) + "_clean";
	TH2F h2_out_clean(hname.Data(),hname.Data(),
			h2_in->GetNbinsX(),
			h2_in->GetXaxis()->GetXmin(), h2_in->GetXaxis()->GetXmax(),
			h2_in->GetNbinsY(),
			h2_in->GetYaxis()->GetXmin(), h2_in->GetYaxis()->GetXmax());
	h2_out_clean.SetDirectory(0);
	for(int xbin=1; xbin<=h2_in->GetNbinsX(); ++xbin){//clean up low stat bins
		for(int ybin=1; ybin<=h2_in->GetNbinsY(); ++ybin){
			double binContent = h2_in->GetBinContent(xbin, ybin);
			h2_out_clean.SetBinContent(xbin, ybin,(binContent < thr) ? 0.0 : binContent); 
		}
	}
	return h2_out_clean;
}

void FootOnline::Calculate_Pedestals(int threshold)
{
	cout << "\n-- Calculating pedestals";
	std::ofstream fout(pedfilename);
	if (!fout.is_open()){
		std::cerr << "Error opening file: " << pedfilename << std::endl;
		return;
	}
	TF1 foo("foo","gaus",50,700);
	for(int i=0; i<16; i++)
	{
		if(h2_raw[i]->Integral()<1) continue;//do not remove this cond
		auto h2_raw_clean = Get_Clean_TH2F(h2_raw[i], threshold);
		h2_raw_clean.SetDirectory(0);
		gErrorIgnoreLevel = kFatal;//to suppress printout info on the terminal
		h2_raw_clean.FitSlicesY(&foo,1,h2_raw_clean.GetNbinsX(),0,"QNR",0);
		gErrorIgnoreLevel = kPrint; // Reset to default level
		//Generated automatically, ignore:	
		auto h_ignore0 =  dynamic_cast<TH1D*> (gDirectory->Get(TString(h2_raw_clean.GetTitle()) +"_0"));
		h_ignore0->SetDirectory(0);//do not store for online
		auto h_ignore_chi2 =  dynamic_cast<TH1D*> (gDirectory->Get(TString(h2_raw_clean.GetTitle()) +"_chi2"));
		h_ignore_chi2->SetDirectory(0);//do not store for online
		h1_peds[i] =  dynamic_cast<TH1D*> (gDirectory->Get(TString(h2_raw_clean.GetTitle()) +"_1"));
		h1_peds[i]->SetDirectory(0);//do not store for online
		h1_sigma[i] = dynamic_cast<TH1D*> (gDirectory->Get(TString(h2_raw_clean.GetTitle()) +"_2"));
		h1_sigma[i]->SetDirectory(0);//do not store for online
		fout << i << std::endl;
		for(int j=0; j<640; j++){
			pedestal[i][j] = h1_peds[i]->GetBinContent(j+1);
			sigma[i][j]    = h1_sigma[i]->GetBinContent(j+1);
			fout << j+1 << "  " <<  pedestal[i][j] << "  " <<  sigma[i][j] << std::endl;
		}
	}
	std::cout << "\n-- Writing new pedestals and sigmas into file " << pedfilename << std::endl;
	fout.close();
}

void FootOnline::Calculate_Fine_Sigmas(int threshold)
{
	cout << "\n-- Calculating fine sigmas";
	TF1 foo("foo","gaus",-200,400);
	for(int i=0; i<16; i++)
	{
		if(h2_fine[i]->Integral()<1) continue;//do not remove this cond
		auto h2_clean = Get_Clean_TH2F(h2_fine[i], threshold);
		h2_clean.SetDirectory(0);
		gErrorIgnoreLevel = kFatal;//to suppress printout info on the terminal
		h2_clean.FitSlicesY(&foo,1,h2_clean.GetNbinsX(),0,"QNR",0);
		gErrorIgnoreLevel = kPrint; // Reset to default level
		//Generated automatically, ignore:	
		auto h_ignore0 = dynamic_cast<TH1D*> (gDirectory->Get(TString(h2_clean.GetTitle()) +"_0"));
		h_ignore0->SetDirectory(0);
		auto h_ignore1 = dynamic_cast<TH1D*> (gDirectory->Get(TString(h2_clean.GetTitle()) +"_1"));
		h_ignore1->SetDirectory(0);
		auto h_ignore_chi2 = dynamic_cast<TH1D*> (gDirectory->Get(TString(h2_clean.GetTitle()) +"_chi2"));
		h_ignore_chi2->SetDirectory(0);
		h1_sigma_fine[i] = dynamic_cast<TH1D*> (gDirectory->Get(TString(h2_clean.GetTitle()) +"_2"));
		h1_sigma_fine[i]->SetDirectory(0);//do not store for online
		for(int j=0; j<640; j++){
			fine_sigma[i][j]    = h1_sigma_fine[i]->GetBinContent(j+1);
		}
	}
	is_fine_sigma=true;//set flag to use fine sigmas threshold
}

void FootOnline::Draw_Pedestals(Bool_t is_fine_sigma)
{
	std::cout << "\n-- Plotting pedestals for all FOOTs\n";
	for(int i=0;i<NDETS;++i)
	{
		int d = foot_id[i]-1;
		if(h2_raw[d]->Integral()<1) continue;//don't remove this cond!
		c_ped->cd(i*2+1);
		h2_raw[d]->Draw("colz");
		Draw_Asic_Borders(h2_raw[d]->GetYaxis()->GetXmin(), 
				h2_raw[d]->GetYaxis()->GetXmax());
		h1_peds[d]->SetMarkerStyle(kFullCircle);
		h1_peds[d]->SetMarkerSize(0.01);
		h1_peds[d]->SetMarkerColor(kRed);
		h1_peds[d]->SetLineColor(kRed);
		h1_peds[d]->Draw("same");		
		c_ped->cd(i*2+2);
		h1_sigma[d]->SetMarkerStyle(kFullCircle);
		h1_sigma[d]->SetMarkerSize(0.01);
		h1_sigma[d]->SetMarkerColor(kBlue);
		h1_sigma[d]->SetLineColor(kBlue);
		h1_sigma[d]->GetYaxis()->SetRangeUser(0,15);
		h1_sigma[d]->GetXaxis()->SetTitle("Strip No.");
		h1_sigma[d]->GetYaxis()->SetTitle("ADC sigma");
		TString htitle;
		htitle.Form("Sigmas for FOOT%d",foot_id[i]);
		h1_sigma[d]->SetTitle(htitle.Data());
		h1_sigma[d]->Draw();
		if(is_fine_sigma)
		{	
			h1_sigma_fine[d]->SetMarkerStyle(kFullCircle);
			h1_sigma_fine[d]->SetMarkerSize(0.01);
			h1_sigma_fine[d]->SetMarkerColor(kRed);
			h1_sigma_fine[d]->SetLineColor(kRed);
			h1_sigma_fine[d]->Draw("same");
		}
		Draw_Asic_Borders(0, 15);
	}
}

int main(Int_t argc, Char_t* argv[])
{
	cout << "\n\n--------------------------------";
	cout <<   "\n  FOOT  online monitor   ";
	cout <<   "\n-------------------------------- \n\n";
	gROOT->SetBatch(kTRUE);//only for online
	gROOT->Macro("include/rootlogon.C");
	gStyle->SetPalette(kRainBow);

	EXT_STR_h101 ucesb_struct;
	ext_data_clnt fClient;
	Bool_t status;
	std::ostringstream command;
	//command << "$UCESB_DIR/../upexps/202402_s091_s118/202402_s091_s118 --stream=lxlanddaq01:9000 --allow-errors --input-buffer=300Mi --ntuple=RAW,FOOT,STRUCT,-";
	command << "$UCESB_DIR/../upexps/202402_s091_s118/202402_s091_s118 /lustre/r3b/202402_s118/lmd_stitched/main0031_0001_stitched.lmd --allow-errors --input-buffer=300Mi --ntuple=RAW,FOOT,STRUCT,-";

	/* Fork off ucesb (calls fork() and pipe()) */
	FILE* fFd = popen(command.str().c_str(), "r");
	if (nullptr == fFd){
		cout << "\npopen() failed\n";
		return 0;
	}
	/* Connect to forked instance */
	status = fClient.connect(fileno(fFd));
	if (kFALSE == status){
		cout << "\next_data_clnt::connect() failed";
		cout << "\nucesb error: " << fClient.last_error();
		return 0;
	}
	uint32_t struct_map_success = 0;
	auto fStructInfo = new ext_data_struct_info();
	int ok;
	size_t fOffset = offsetof(EXT_STR_h101, foot);
	EXT_STR_h101_FOOT_ITEMS_INFO(ok, *fStructInfo, fOffset, EXT_STR_h101_FOOT, 0);

	if (!ok){
		perror("ext_data_struct_info_item");
		fprintf (stderr,"Failed to setup structure information.\n");
		exit(1);
	}
	status = fClient.setup(NULL, 0, fStructInfo, &struct_map_success, sizeof(ucesb_struct));
	if (status != 0){
		cout << "\next_data_clnt::setup() failed";
		cout << "\nUCESB error: " << fClient.last_error();
		return 0;
	}

	if(!Read_Pedestals()){
		std::cerr << "ERROR!  Cannot read pedestals!" <<  std::endl;
		return 0;
	}
	FootOnline* FootAna = new FootOnline();	
	auto serv = new THttpServer("http:8886");//setup web display
	//gROOT->SetWebDisplay("browser"); // or "cef" for embedded browser
	serv->Register("",dynamic_cast<FootOnline*> (FootAna));
	serv->RegisterCommand("Reset_Histograms", Form("/Objects/FootAna/->Clear()"));
	//serv->RegisterCommand("Reset_Histograms", Form("FootAna/->Clear()"));
	int ret;
	const void* raw;
	ssize_t raw_words;
	int Nevents=0;
	int asic_pitch = 64;//chunk of strips, 64 in 1 asic
	const int Nasic = 640/asic_pitch;
	double  asic_offset[Nasic];
	double  signal = 0;
	int counter_asic =0;
	int stat=0;
	double signal_raw;
	int strip;
	//clustering
	foot_data left_arm_fdata[4];
	foot_data right_arm_fdata[4];
	hits left_hits[4];
	hits right_hits[4];
	bool is_good_event;//at least one event
	while(1)//Eventloop
	{
		is_good_event = false;
		//cout << "\n\n########### NEW EVENT ##########";
		for(auto f=0; f<4; ++f){//clear all cluster data
			left_arm_fdata[f].clear(); 
			right_arm_fdata[f].clear(); 
			left_hits[f].clear();
			right_hits[f].clear();
		}
		ret = fClient.fetch_event(&ucesb_struct, sizeof(ucesb_struct));
		if (-1 == ret){
			cout << "\next_data_clnt::fetch_event() failed! \nUCESB error: " << fClient.last_error() << endl;
			return 0;
		}
		ret = fClient.get_raw_data(&raw, &raw_words);
		if (0 != ret){
			cout << "\nFailed to get raw data.\n";
			return 0;
		}
		for(int d=0;d<16;d++)
		{
			if(ucesb_struct.foot.FOOT[d]._ !=640) continue;
			stat=0; counter_asic=0;
			for(int i=0; i<Nasic; i++){  asic_offset[i]=0.; }//reset asic baselines
			for(int  j=0 ; j<ucesb_struct.foot.FOOT[d]._; j++)
			{
				strip = ucesb_struct.foot.FOOT[d].I[j];
				signal_raw = ucesb_struct.foot.FOOT[d].E[j];
				FootAna->h2_raw[d]->Fill(strip, signal_raw); 

				signal = signal_raw - pedestal[d][j];
				FootAna->h2_calib[d]->Fill(strip, signal);
				if(fabs(signal) < (6 * sigma[d][j]) && 
						is_good_strip(d+1,strip)){//careful with d+1 here
					stat++;
					asic_offset[counter_asic] += signal;
				}
				if((strip % asic_pitch)==0 && strip>1){//switch to next asic
					if(stat>1) asic_offset[counter_asic] /= stat;
					counter_asic++;  stat=0;
				}
			}
			//----- fine baseline correction
			counter_asic=0;  stat=0;
			for(int  j=0 ; j<ucesb_struct.foot.FOOT[d]._; j++)
			{
				strip = ucesb_struct.foot.FOOT[d].I[j];
				signal_raw = ucesb_struct.foot.FOOT[d].E[j];
				if((strip % asic_pitch) == 1 && strip>1) counter_asic++; 


				if(is_asic_edge_strip(strip))
				{
					signal = signal_raw - pedestal[d][j];//simple pedestal subtraction for those strips
				}
				else
				{
					signal = signal_raw - pedestal[d][j] - asic_offset[counter_asic];
				}
				FootAna->h2_fine[d]->Fill(strip, signal);

				//Clustering
				if(!is_good_strip(d+1,strip)) continue;
				if( FootAna->is_fine_sigma  && signal < (NSIGMA * fine_sigma[d][j]) ) continue;
				else if( (!FootAna->is_fine_sigma) && signal < (NSIGMA * sigma[d][j]) ) continue;
				is_good_event=true;
				for(auto f=0; f<NDETS; ++f)
				{
					if((d+1) == foot_id[f])
					{
						//cout << "\nGood signal in " << foot_id[f] << "detector, strip=" << strip;
						if(f<4)	Check_Strip(strip, signal, left_arm_fdata[f]);
						else Check_Strip(strip, signal, right_arm_fdata[f-4]);
					}
				}
			}//fine baseline correction
		}//16 detectors
		Nevents++;
		//Transform to hits 
		cluster_data this_clust;
		for(int k=0; k<4; ++k)//iterate 4 detectors on left/right
		{
			//if(is_good_event)
			//{
			//	cout << "\n\nDet: left = " << k+1 << ", Hit mul = " << left_arm_fdata[k].size();
			//	cout << "\nDet: right = " << k+1 << ", Hit mul = " <<  right_arm_fdata[k].size();
			//}
			for(auto & c: left_arm_fdata[k]){
				this_clust.position = get_cog(c) * FOOT_LENGTH/640 - 0.5*FOOT_LENGTH;
				this_clust.energy = get_esum(c);
				this_clust.mul = c.size();
				left_hits[k].push_back(this_clust);
				//cout << "\nPosition = " << this_clust.position << ", energy = " <<  this_clust.energy << " mul =" << this_clust.mul;
			}
			for(auto & c: right_arm_fdata[k]){
				this_clust.position = get_cog(c) * FOOT_LENGTH/640 - 0.5*FOOT_LENGTH;
				this_clust.energy = get_esum(c);
				this_clust.mul = c.size();
				//cout << "\nPosition = " << this_clust.position << ", energy = " <<  this_clust.energy << " mul =" << this_clust.mul;
				right_hits[k].push_back(this_clust);
			}
		}
		//Fill icorrelation plots
		for(auto & h1: left_hits[0]){
			for(auto & h2: left_hits[2]){
				//cout << "\nFill corr";
				FootAna->h2_left_X_corr->Fill(h1.position, h2.position);
			}
		}
		for(auto & h1: left_hits[1]){
			for(auto & h2: left_hits[3]){
				//cout << "\nFill corr";
				FootAna->h2_left_Y_corr->Fill(h1.position, h2.position);
			}
		}
		for(auto & h1: right_hits[0]){
			for(auto & h2: right_hits[2]){
				//cout << "\nFill corr";
				FootAna->h2_right_X_corr->Fill(h1.position, h2.position);
			}
		}
		for(auto & h1: right_hits[1]){
			for(auto & h2: right_hits[3]){
				//cout << "\nFill corr";
				FootAna->h2_right_Y_corr->Fill(h1.position, h2.position);
			}
		}

		if((Nevents % 5000)==0)//update every N events
		{
			FootAna->c_raw->Update();
			FootAna->c_calib->Update();
			FootAna->c_fine->Update();
			FootAna->c_corr_xy->Update();
			gSystem->ProcessEvents();

		}
		if((Nevents%500000) == 0)
		{
			FootAna->Calculate_Pedestals(10);
			FootAna->Calculate_Fine_Sigmas(10);
			FootAna->Draw_Pedestals(kTRUE);
		}
	}//while
	for(int i=0;i<NDETS;++i)
	{
		int d = foot_id[i]-1;
		if(FootAna->h2_fine[d]->Integral()<1) continue;//don't remove this cond!
		FootAna->c_fine->cd(i+1);//displayed in browser
		FootAna->h2_fine[d]->Draw("colz");
		Draw_Asic_Borders(FootAna->h2_fine[d]->GetYaxis()->GetXmin(),
				FootAna->h2_fine[d]->GetYaxis()->GetXmax());
	}
	return 0;
}
