/* basic sorter. 
First run trace graber(./tgetwf -ip your.ip.adress -f filename -c channel -r numruns), provided then convert binary file to hex using unix/linux commmand: hexdump -x < binaryfile.wf > hextextfile.txt
each line should follow this format. if it doesn't this string sorter won't work and need to make appropriate adjustments
0000000    4f00    5000    5100    5100    5300    5400    5100    5000
0000010    5100    5200    5300    5100    5000    5200    5000    5100
0000020    4f00    4f00    5000    4d00    5000    4e00    4f00    4d00
0000030    4f00    4f00    4d00    4d00    4f00    4d00    4c00    4e00

Each line should be 71 characters long, if not it will be dismissed (sometimes a line is emty or filled with ** garbage). The vxi-11/tekVISA binary transfer protocol (2byte 8 bit signed).
these are formated into a vector like this (first line from above): 4f    50    51    51    53    54    51    50
then there are converted to 8bit signed intergers, then adjusted based on parameters from the .wfi file provided by trace graber (./tgetwf)

For 1PE: average taken from part A1 results


-26V  -3.5mV
-26.5 -4.5mV
-27v  -6mV

*/
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <TLatex.h>
#include <TImage.h>
#include "TF1.h"
#include "TH2F.h"
#include <TRandom3.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cstdlib>
#include <TPad.h>
#include <TFrame.h>
#include <cmath>
#include <TGraph2D.h>
#include "TGraph.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFile.h"
#include "TPaveStats.h"
#include <TLegend.h>
#include <TLatex.h>

using namespace std;

template <typename T> string tostr(const T& t) { 
   ostringstream os; 
   os<<t; 
   return os.str();
}




int run() 
{  TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);



	 const char* c; 
	string titlefile2 = "6mm C-series Single PE Response px2; u/U; delta"; 
	c = titlefile2.c_str();	
	TGraph *gr = new TGraph();
	gr->SetTitle(c);
	
			TGraph *gr1 = new TGraph();
			gr1->SetMarkerColor(4);
			gr1->SetMarkerStyle(7);
			gr1->SetTitle(c);	

			TGraph *gr2 = new TGraph();
			gr2->SetMarkerColor(5);
			gr2->SetMarkerStyle(7);
			gr2->SetTitle(c);

			TGraph *gr3 = new TGraph();
			gr3->SetMarkerColor(6);
			gr3->SetMarkerStyle(7);
			gr3->SetTitle(c);

			TGraph *gr4 = new TGraph();
			gr4->SetMarkerColor(7);
			gr4->SetMarkerStyle(7);
			gr4->SetTitle(c);

			TGraph *gr5 = new TGraph();
			gr5->SetMarkerColor(9);
			gr5->SetMarkerStyle(7);	

float t = 0;
float omega = 1.0;
float vis = 1.0;
float u;
float y = 0;
float PI = 3.14159;		
			
			//Data filling *********************************
float h = 8.0*sqrt(vis/omega);
	

				
			for (int i = 0; y < h; ++i)
			{		
                   
                    y = i*(0.001)*sqrt(omega/(2.0*vis));
                    t = 0.0;
					u = exp(-y*sqrt(omega/(2.0*vis)))*cos(omega*t-y*sqrt(omega/(2.0*vis)));
					gr1->SetPoint(i,u,y);		
		
				
                    t = (PI/2.0)/omega;
					u = exp(-y*sqrt(omega/(2.0*vis)))*cos(omega*t-y*sqrt(omega/(2.0*vis)));
					gr2->SetPoint(i,u,y);		
		
		
                    t = (3.0*PI/2.0)/omega;
					u = exp(-y*sqrt(omega/(2.0*vis)))*cos(omega*t-y*sqrt(omega/(2.0*vis)));
					gr3->SetPoint(i,u,y);		
	
	
                    t = (PI)/omega;
					u = exp(-y*sqrt(omega/(2.0*vis)))*cos(omega*t-y*sqrt(omega/(2.0*vis)));
					gr4->SetPoint(i,u,y);		
	


                    t = (PI)/omega;
					u = exp(-y*sqrt(omega/(2.0*vis)))*cos(omega*t-y*sqrt(omega/(2.0*vis)));
					gr5->SetPoint(i,u,y);		
			}


c1->SetGridx();
c1->SetGridy();
c1->SetTickx(1);
c1->SetTicky(1);

gr->SetPoint(0,1,h);
gr->SetPoint(1,-1,0);


    gr->Draw("ap");

    gr1->Draw("samep");
	//gr1->GetYaxis()->SetRangeUser(0,yaxis);
   // gr5->Draw("samep");

	//gr2->GetYaxis()->SetRangeUser(0,yaxis);
    gr2->Draw("samep");

	//gr3->GetYaxis()->SetRangeUser(0,yaxis);
    gr3->Draw("samep");
 
	//gr4->GetYaxis()->SetRangeUser(0,yaxis);
    gr4->Draw("samep");
	

	TLegend *leg = new TLegend(0.6,0.3,0.7,0.55);

leg->AddEntry(gr1,"t = 0","AP");
leg->AddEntry(gr2,"t = PI/2","AP");
leg->AddEntry(gr3,"t = 3PI/2","AP");
leg->AddEntry(gr4,"t = PI","AP");

leg->SetFillColor(0);

leg->Draw();


//---------------------------------------------2d graph --------------------------------//

TCanvas *c2 = new TCanvas("c2","The FillRandom example",200,50,900,700);
TCanvas *c4 = new TCanvas("c4","The FillRandom example",200,50,900,700);

TGraph2D *gr11 = new TGraph2D();
TGraph2D *gr55 = new TGraph2D();

string titlefile = "Energy Explicit flow1; x; y; z";
c = titlefile.c_str();	
gr11->SetTitle(c);




titlefile = "Energy Implicit flow2; x; y; z";
c = titlefile.c_str();
gr55->SetTitle(c);	



	
			


int N = 0;
t = 0.0;
y = 0.0;
for (int j = 1; j < 1000; j++) 
{t = t + 0.01;
	
	for (int i = 0; y < h; ++i)
	{	
		
      y = i*(0.001)*sqrt(omega/(2.0*vis));               					
	  u = exp(-y*sqrt(omega/(2.0*vis)))*cos(omega*t-y*sqrt(omega/(2.0*vis)));

      gr11->SetPoint(N,y,t,u);
      gr55->SetPoint(N,y,t,u);
      ++N;
	 }     
}




/*
c2->cd();   
   gStyle->SetPalette(1);
   gr1->Draw("colz");
   
c4->cd();   
   gStyle->SetPalette(1);
   gr5->Draw("colz");   */
   
   c2->cd();   
   gStyle->SetPalette(1);
   gr11->Draw("colz");
   
   c4->cd();   
   gStyle->SetPalette(1);
   gr55->Draw("surf1z");   
   	
	

return 1;		
}


