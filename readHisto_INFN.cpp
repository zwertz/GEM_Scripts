/*
 * Read mpdLibtest Histo mode output file
 *
 * Mar2021: added analysis of the Sample mode output file to get optimal phase_space and FIR parameters 
 * 
 * Last update: Mar/2021
 * Author: E. Cisbani
 *
 */

#include <unistd.h>
#include <iostream>
#include <fstream>
//#include <Riostream.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TMath.h>
#include <TProfile.h>
#include <TString.h>
#include <TText.h>
#include <TCut.h>
#include <TLine.h>

/* valid before 2020 - Oct - 28
// APV flag = 1 present, 0 = empty ADC channel, added Aug/2019
Int_t flagon[16][16]={
{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0}, // 1st chamber
{0,1,1,1,1,0,1,1,1,1,1,1,1,1,0,0},
{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0},
{1,1,1,1,0,0,1,1,1,1,0,1,1,1,1,0},
{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0}, // 2nd chamber
{1,1,1,1,0,0,1,1,1,1,1,1,1,1,0,0}, 
{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0},
{1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0}, 
{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0}, // 3rd chamber
{1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0}, 
{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0},
{1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0}, 
{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0}, // 4rd chamber
{1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0}, 
{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0},
{1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0} 
};
*/

// APV flag = 1 present, 0 = empty ADC channel, since 24 / Mar /2022
Int_t flagon[16][16]={
  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0}, // J1 Layer
  {1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0},
  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0},
  {1,1,1,1,0,0,1,1,1,1,1,1,1,1,0,0},
  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0}, // J3 Layer
  {1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0}, 
  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0},
  {1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0}, 
 /* {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0}, // 3rd chamber
  {0,1,1,1,1,1,1,1,1,0,1,1,1,1,0,0}, 
  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0},
  {1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0}, 
  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0}, // 4rd chamber
  {1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0}, 
  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0},
  {1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0} 
  };*/

};


void setStyle() {

  gStyle->SetCanvasColor(kWhite);     // background is no longer mouse-dropping white
  gStyle->SetCanvasDefW(500);
  gStyle->SetCanvasDefH(500);
  gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
  gStyle->SetCanvasBorderMode(0);     // turn off canvas borders
  gStyle->SetPadBorderMode(0);
  gStyle->SetPaintTextFormat("5.2f");  // What precision to put numbers if plotted with "TEXT"

  gStyle->SetTextSize(0.08);
  gStyle->SetLabelSize(0.1,"xyz"); // size of axis value font
  gStyle->SetTitleSize(0.1,"xyz"); // size of axis title font
  gStyle->SetTitleFont(62,"xyz"); // font option
  gStyle->SetLabelFont(62,"xyz");
  gStyle->SetTitleOffset(0.8,"y"); 
  gStyle->SetTitleOffset(1.0,"x");

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.05);

  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

/*
 *
 * infile: output of the sbsvme20:/home/daq/ben/mpd/libsrc4.0/test/mpdLibTest out.txt 4
 *
 * oflag=0 : visualize one canvas for each MPD hists
 *       1 : one canvas for all histos, each MPD is a column
 *       2 : as 1 and save canvas in a file
 *      negative values change the order of the visualization
 */
int readHisto_INFN(TString infile, TString outpdf, Int_t oflag=0) {

  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.05);

  Int_t i;
   
  // clean input file and read relevant information

  /*
    gSystem->Exec(Form("awk '$1 ~ /^[#.0-9]+$/ { print $0 }' %s >temp.txt",infile.Data()));
    TString dummy = gSystem->GetFromPipe("awk 'BEGIN {flg=0;} $2==\"Time=\"&&flg==0 { print $3; flg=1; }' temp.txt");
    Int_t time_off = dummy.Atoi();
    printf("Start absolute time = %d sec\n",time_off);
    dummy =  gSystem->GetFromPipe("awk 'BEGIN {flg=0;} $1!=\"#\"&&flg==0 { print NF; flg=1; }' temp.txt");
    Int_t nch = (dummy.Atoi()-2)/3;
    printf("Number of channel(s) in input file = %d\n",nch);
    TString stime = gSystem->GetFromPipe("awk 'BEGIN {flg=0;} $2==\"Time_Start=\"&&flg==0 { $1=$2=\"\"; print $0; flg=1; }' temp.txt");
    stime = stime.Strip(TString::kBoth);
    printf("Start Time = %s\n",stime.Data());
  */
  
  TTree *thm = new TTree("thm","MPD Histo Mode data");

  TString rform="mpd/I:adc/I:gain/I:val[4096]/I";

  Int_t nrow = thm->ReadFile(infile.Data(),rform);

  printf("%d rows of data from input file %s\n", nrow, infile.Data());

  //  thm->Draw("val:Iteration$","(mpd==6)&&(adc==0)");

  thm->SetMarkerStyle(7);

  // count mpd
  thm->Draw("mpd","adc==0");
  Int_t nmpd = thm->GetSelectedRows();

  printf("Number of mpds in data = %d\n",nmpd);
  
  // count mpd x adc
  thm->Draw("mpd:adc");

  Int_t nscan = thm->GetSelectedRows();

  printf("Number of total scanned channels  = %d\n",nscan);

  Int_t *impd = new Int_t[nscan];
  Int_t *iadc = new Int_t[nscan];

  for (Int_t i=0;i<nscan;i++) {
    impd[i] = thm->GetV1()[i];
    iadc[i] = thm->GetV2()[i];
    //  printf("%d : %d %d\n",i,impd[i],iadc[i]);
  }

  // MPD index = VME slot where the MPD is inserted
  // Connector = Analog HDMI Connector on MPD, 0=upper, 3=lower
  // Backplane = 0, 1 and 3: is the backplane that host 4 or 5 cards
  // ADC: is the ADC index in the MPD corresponding to a single APV card (128 channels)

  TCanvas **cc;
  if (oflag==0) {
    cc = new TCanvas*[nmpd];
    for (Int_t i=0;i<nmpd;i++) {
      printf("Prepare canvas %d for mpd %d\n",i,impd[i]);
      cc[i] = new TCanvas(Form("mpd%d",impd[i]),Form("MPD %d : Connector # Backplane * ADC (anomaly)",impd[i]));
      cc[i]->Divide(4,4,0.005,0.005);
      cc[i]->Update();
    }
  } else {
    cc = new TCanvas*[1];
    cc[0] = new TCanvas("mpda","MPD : Connector # Backplane * ADC (anomaly)");
    cc[0]->Divide(16,16,0.001,0.001);
    cc[0]->Update();
  }
  
  // now plot

  TGraph **gr = new TGraph*[nscan];
  Int_t cabcol[4]={2,4,8,20};
  Int_t j=0; // mpd idx
  Int_t anomaly_flag=0;

  TGraphErrors *grsta[3]; // region statistics

  for (Int_t i=0;i<3;i++) {
    grsta[i]=new TGraphErrors(nscan);
  }
  
  for (Int_t i=0;i<nscan;i++) {
    // printf("%d %d : %d %d\n",i,j,impd[i],iadc[i]);
    Int_t cidx=j;
    Int_t pidx=1+iadc[i];
    
    if (oflag!=0) {
      cidx=0;
      pidx=1+iadc[i]+j*16;
      if (oflag<0) pidx=1+iadc[i]*16+j;
    }
    
    cc[cidx]->cd(pidx)->SetLogy(1);


    thm->Draw("val:Iteration$",Form("(mpd==%d)&&(adc==%d)",impd[i],iadc[i]));
    Int_t nnn = thm->GetSelectedRows();
    gr[i] = new TGraph(thm->GetSelectedRows(), thm->GetV2(), thm->GetV1());

    gr[i]->SetMinimum(1);
    gr[i]->SetMarkerStyle(6);
    gr[i]->SetMarkerColor(kBlack);
    gr[i]->SetLineColor(20);
    gr[i]->SetTitle("Empty");
    if (flagon[j][iadc[i]]) {  // card expected in this position
      // try to identify anomalies in the three main regions: "near zero", "intermediate", "around 1"
      Double_t rms[3]={0,0,0};
      Double_t mean[3]={0,0,0};
      Double_t count[3]={0,0,0};
      Int_t ir=0;
      for (Int_t ig=0;ig<nnn;ig++) {
	Double_t x,y;
	gr[i]->GetPoint(ig, x, y);
	if (x<700) ir=0; 
	if ((x>=700) && (x<2550)) ir=1;
	if (x>=2550) ir=2;
	mean[ir] += y*x;      
	rms[ir] += y*x*x;
	count[ir] +=y;
      }

      for (ir=0;ir<3;ir++) {
	if (count[ir]) {
	  rms[ir] = TMath::Sqrt((rms[ir] - mean[ir]*mean[ir]/count[ir])/count[ir]);
	  mean[ir] = mean[ir]/count[ir];
	}
	Int_t ix = impd[i]*16+iadc[i];
	grsta[ir]->SetPoint(i, ix, mean[ir]);
	grsta[ir]->SetPointError(i, 0.1, rms[ir]);
	//	printf("%2d %3d : %1d : %.2f %.1f \n", j,i, ir, mean[ir], rms[ir]);
      }
      anomaly_flag=0;
      if (mean[1]>0) anomaly_flag=1;
      if (rms[2]>10) anomaly_flag +=2;
      if (mean[2]<2550) anomaly_flag +=4;
      
      Int_t cacos=cabcol[iadc[i]/5];
      Int_t bckcol=19;
      if (anomaly_flag) bckcol=41;

      cc[cidx]->cd(pidx)->SetFillColorAlpha(bckcol,.22);

      gr[i]->SetLineColor(cacos);
      gr[i]->SetMarkerColor(cacos);
      gr[i]->SetTitle(Form("%d : %d # %d > %d (%x)",impd[i],iadc[i]/4,iadc[i]/5,iadc[i], anomaly_flag));
    }
    gr[i]->Draw("pawl");
    cc[cidx]->Update();
    j++;
    j = j % nmpd;
  }

  TCanvas *cs = new TCanvas("cs","Summary");
  cs->Divide(1,3);
  cs->Update();

  for (Int_t i=0;i<3;i++) {
    cs->cd(i+1);
    grsta[i]->SetMarkerStyle(20+i);
    grsta[i]->SetMarkerSize(2);
    
    grsta[i]->Draw("PAW");
  }
  cs->Update();

  if (TMath::Abs(oflag)==0) {
    for(int icanvas=0; icanvas < nmpd; icanvas++){
      if(icanvas == 0) cc[icanvas]->SaveAs("plots/" + outpdf + "(");
      else if(icanvas == nmpd-1) cc[icanvas]->SaveAs("plots/" + outpdf + ")");
      else cc[icanvas]->SaveAs("plots/" + outpdf);
      
    }
  }

  return 0;

}

/*****
 */

 /*
  * Porting of the FIR coefficients estimation code of the BELLE II note B2N/0007 (M. Friedl et al., 18/Mar/2011) 
  *
  * Find FIR filters by Matrix inversion
  *
  * number of parameters has generalized (from 8 to ncoeff)
  *
  * ncoeff : number of coefficients of the FIR filter
  * vsample: array of the impulse sampled values (number of elements must be >= ncoeff), the first element must be the maximum. 
  * hmin : baseline value (corresponds to vsample[34] or so)
  * hdel : normalization so that the sampled distribution is normalized to 1 (corresponds to the max of vsample)
  *
  * this method is called by DQsample - you generally do not need to call it directly 
  *
  * Return the string of FIR coefficients
  *
  */

TString estimateFIR(Int_t ncoeff, Float_t *vsample, Float_t hmin, Float_t hdel) {

  // calculate fir coefficients (from Belle II)
  Double_t h[ncoeff];
  for (int i=0;i<ncoeff;i++) {
    h[i] = (vsample[i]-hmin)/hdel;
    cout << h[i] << " ";
  }
  
  cout << endl;

  // initialize matrix H:
  Int_t j, k, l;
  Double_t H[ncoeff*ncoeff];
  Double_t hF[ncoeff], vdelta[ncoeff];
  Double_t firconst[ncoeff];
  Double_t sumcoeff=0.,sum=0.;
  
  for (j=0; j<ncoeff; j++) {
    for (k=0; k<=j; k++) {
      H[j*ncoeff+k]= h[j-k];
    }
    for (k=j+1; k<ncoeff; k++) {
      H[j*ncoeff+k]=0.;
    }
    hF[j]=0;
    vdelta[j]=0;
  }
  vdelta[0]=1.;
  
  // calculate filter coefficients (recursive inversion of
  // triangular matrix):
  sumcoeff=0.;
  for (j=0; j<ncoeff; j++) {
    sum=0;
    for (l=0; l<j; l++) {
      sum = sum+H[j*ncoeff+l]*hF[l];
    }
    hF[j] = 1./H[j*ncoeff+j]*(vdelta[j]-sum);
    sumcoeff += hF[j];
  }
  // normalize filter coefficients (sum of all coefficients := 1)
  cout << "## FIR Constants: ";
  
  for(j=0; j<ncoeff; j++) {
    firconst[j]=hF[j]/sumcoeff;
    cout << firconst[j] << " "; 
  }
  cout << endl;
  
  // scaled to 2^13
  TString sret="[ ";
  
  cout << " fir_coeff = [ ";
  for(j=0; j<ncoeff; j++) {
    //  cout << firconst[j]*pow(2,13) << " "; 
    printf("%d",(int) (firconst[j]*pow(2,13))); 
    sret=sret + " " + Form("%d", (int) (firconst[j]*pow(2,13)));
    if (j==(ncoeff-1)) { printf(" ];\n"); } else { printf(", "); }
  }
  sret = sret + " ]";
  return sret;

}

TString wrapSyncPulse(const double *v_adc, int nsample=1024) {

  const int cycle = 35; // sync pulse period [clock]
  
  float vsample[cycle];
  float vnormal[cycle];
  
  for (int i=0;i<cycle;i++) {
    vsample[i] = 0;
    vnormal[i] = 0;
  }
    
  for (int i=0;i<nsample;i++) {
    vsample[i % cycle] += v_adc[i];
    vnormal[i % cycle] += 1.;
  }

  float smax=-1;
  int imax = -1;

  for (int i=0;i<cycle;i++) { // find maximum sample
    float ratio = vsample[i]/vnormal[i];
    if (smax < ratio) {
      imax = i;
      smax = ratio;
    }
    vsample[i]=ratio;
  }

  float vonorm = smax;
  float vobase = 0;
  float vosam[cycle];
  float vclk[cycle];
  
  for (int i=0;i<cycle;i++) { // shift max at 0 and estimate baseline
    vosam[i] = vsample[(imax+i)%cycle];
    vclk[i] = i*25.;
    if (i>=(cycle-11)&&(i<(cycle-3))) { // assume the baseline is the tail, remove last two points
      vobase +=vosam[i];
    }
  }
  vobase = vobase/8.;

  TGraph *cgra = new TGraph(cycle,vclk,vosam);
  cgra->SetMarkerStyle(20);
  cgra->SetMarkerColor(2);
  cgra->Draw("PAW");
  cgra->SetTitle(Form("Sync Pulse Sampled"));
  cgra->GetXaxis()->SetTitle("Time [ns]");
  
  TString cfir;
  cfir = estimateFIR(16, vosam, vobase, vonorm);
  
  return cfir;
  
}

float optimalPhaseSearch(TH2F *h2ap, float midadc) {

  float optdist = -1;
  float optphase = -1;

  int nphase = h2ap->GetXaxis()->GetNbins();

  //  printf(" Number of bins: %d %f\n",nphase, midadc);
  
  for (int i=0;i<nphase;i++) {
    int iph = (int) h2ap->GetXaxis()->GetBinCenter(i);
    TH1D* hhd = (TH1D*)h2ap->ProjectionY("_py",iph,iph);
    hhd->GetXaxis()->SetRange(0,midadc);
    float m0 = hhd->GetMean();
    float s0 = hhd->GetRMS();
    hhd->GetXaxis()->SetRange(midadc,hhd->GetNbinsX());
    float m1 = hhd->GetMean();
    float s1 = hhd->GetRMS();
    if ((m0>0)&&(m1>0)) {
      float dd = (m1-m0)/(s1);
      
      if (dd>optdist) {
	optphase    = (float) h2ap->GetXaxis()->GetBinCenter(i+1);
	optdist     = dd;
      }
      //      printf("  %d : %f %f %f %f : %f (current: %f %f)\n",iph,m0,m1,s0,s1, dd, optphase, optdist);
    }

  }
  
  //  printf("Optimal phase: %f %f)\n", optphase, optdist);
  return optphase;
}

/*
 * read and process the mpdLibTest sample mode output file
 *
 * oflag = 0: all ADC samples are plotted
 * oflag = 1: select upper part of ADC samples and plot them
 */

int readSample(TString infile, TString outpdf, Int_t oflag=0) {

  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.05);

  // read config lines (beginning with "#") first and put in TTree the relevant ones

  Int_t clkph[3]; // first, last and delta clock_phase
  TString line, tok;
  ifstream ifile;
  ofstream ofile;
  ifile.open(infile);
  ofile.open("temp.txt",ofstream::out);
  int count=0;
  while(!ifile.eof()) 
    {
      line.ReadLine(ifile);
      if (line[0]=='#') {
	TString tok;
	Ssiz_t from = 0;
	from = line.Index(":")+2; // move after ":"
	if (line.Index("CLOCK_RANGE")>=0) {
	  int ii=0;
	  while (line.Tokenize(tok, from, "[ ]")) {
	    clkph[ii]=tok.Atoi();
	    ii++;
	  }
	}
	if (line.Index("MPD_ADC_COUNT_CLOCK_TOUT")>=0) {
	  while (line.Tokenize(tok, from, ":")) {
	    ofile<<tok.Data()<<endl;
	  }
	  count++;
	}
      }
      //      if (count == 10) { break; }
    }
  ifile.close();
  ofile.close();

  printf("Clock Phase Loop: %d %d %d\n",clkph[0],clkph[1],clkph[2]);

  TTree *tc = new TTree("tc","MPD Sample Mode Config Data");
  TString rform="mpd/I:ch/I:cnt/I:clk/I:tout/I";
  Int_t nrowc = tc->ReadFile("temp.txt",rform);
  
  TTree *tsm = new TTree("tsm","MPD Sample Mode ADC data");
  rform="adc[1024]/I";
  Int_t nrows = tsm->ReadFile(infile.Data(),rform);

  if (nrowc != nrows) {
    printf("ERROR, Config Rows %d != Data Rows %d\n",nrowc, nrows);
    return -1;
  }

  tc->Draw("mpd","","goff");
  int nn = tc->GetSelectedRows();
  int mpd0 = TMath::MinElement(nn,tc->GetV1());
  int mpd1 = TMath::MaxElement(nn,tc->GetV1());
  int dmpd = mpd1-mpd0+1;

  tc->Draw("cnt",Form("(mpd==%d)",mpd0));
  int nsample = tc->GetV1()[0];
  
  printf("Total Rows : %d, MPD min %d max %d (events = %d)\n",nrowc,mpd0,mpd1,nsample);

  tsm->SetMarkerStyle(7);

  tsm->AddFriend("tc");

  TCanvas *cc[dmpd];
  TCanvas *cpul[dmpd];
  
  TText *txt = new TText(05,0.5,"nan");
  txt->SetTextSize(.05);

  for (int i=0;i<dmpd;i++) {

    TCut c_mpd = Form("(tc.mpd==%d)",i+mpd0);

    tsm->Draw("tc.mpd",c_mpd,"goff");
    
    if (tc->GetSelectedRows()>0) {
      cc[i] = new TCanvas(Form("cc%d",i),Form("MPD %d Sync DataStream Samples",i+mpd0));
      cc[i]->Divide(4,4);
      cc[i]->Update();

      cpul[i] = new TCanvas(Form("cpul%d",i),Form("MPD %d Sync Pulse Reconstructed",i+mpd0));
      cpul[i]->Divide(4,4);
      cpul[i]->Update();
      
      float thr = -1; // disabled
      for (int ic=0;ic<16;ic++) {

	TCut c_ch = Form("(tc.ch==%d)",ic);
      
	cc[i]->cd(ic+1);
	tsm->Draw("adc[]",c_mpd&&c_ch);
	int ns = tsm->GetSelectedRows();
	if (ns>0) {
	  float_t adc0 = TMath::MinElement(ns,tsm->GetV1());
	  float_t adc1 = TMath::MaxElement(ns,tsm->GetV1());
	  thr = (adc1+adc0)/2.;
	  if (adc0==adc1) {
	    printf("WARNING: mpd/ch %d/%d cannot see sync pulse\n",i+mpd0,ic);
	  } else {
	    TCut c_adc="(1==1)"; // disabled
	    if (oflag==1) c_adc = Form("(adc>%f)",thr);
	    tsm->Draw("adc[]:tc.clk",c_mpd&&c_ch&&c_adc,"colz");
	    if ( tsm->GetSelectedRows()>0) {
	      TH2F *htemp2 = (TH2F*) gPad->GetPrimitive("htemp");
	      htemp2->SetXTitle("Clock Phase [0.5 ns]");
	      htemp2->SetYTitle("ADC");

	      // find optimal clock phase
	      TH2F *h2ap=new TH2F("h2ap",Form("ADC (sample) vs Clock Phase, ch %d",ic),64,-0.5,63.5,adc1-adc0+1,adc0-0.5,adc1-0.5);
	      tsm->Draw("adc[]:tc.clk>>h2ap",c_mpd&&c_ch,"goff");
	      float optcphase = optimalPhaseSearch(h2ap, thr);
	      htemp2->SetTitle(Form("ADc vs Phase, Opt. Phase: #sim %.1f #pm 2",optcphase));
	      delete h2ap;

	      // find sync pulse and wrap samples within one sync period (35 clocks)
	      tsm->Draw("adc:Iteration$",c_mpd&&c_ch,"goff");
	      int n_one = tsm->GetSelectedRows();
	      double *v_adc, *v_evt;
	      v_adc = tsm->GetV1();
	      v_evt = tsm->GetV2();

	      cpul[i]->cd(ic+1)->SetLogy();
	      TString stra = wrapSyncPulse(v_adc, nsample);
	      txt->DrawTextNDC(.12,.7,Form("Ch. %d, FIR coeffs",ic));
	      txt->DrawTextNDC(.12,.6,stra.Data());
	      
	      printf(" ### %s\n",stra.Data());
	      
	    }
	  }
	}
      }
      cc[i]->Update();

    }
    
  }

 if (TMath::Abs(oflag)==0) {
    for(int icanvas=0; icanvas < dmpd; icanvas++){
  if(icanvas == 0){ cc[icanvas]->SaveAs("plots/" + outpdf + "(");
      cpul[icanvas]->SaveAs("plots/" + outpdf);}
       else if(icanvas == dmpd-1){ cc[icanvas]->SaveAs("plots/" + outpdf);
       cpul[icanvas]->SaveAs("plots/" + outpdf + ")");}
       else{ cc[icanvas]->SaveAs("plots/" + outpdf);
       cpul[icanvas]->SaveAs("plots/" + outpdf);
       }

  }

 }
  
  return 0;
  
}
