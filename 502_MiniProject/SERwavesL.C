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
#include <TGraph.h>
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
	string titlefile2 = "6mm C-series Single PE Response px2; Time (ns); Amplitude (mV)"; 
	c = titlefile2.c_str();	
	
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
float omega = 100.0;
float vis = 1.0;
float y, u;

		
			
			//Data filling *********************************
	
			float thump = 0;
				
			for (int i = 0; i < 1000; ++i)
			{		y = i*0.01;
                    t = 0.1;
					u = exp(-y*sqrt(omega/(2.0*vis)))*cos(omega*t-y*sqrt(omega/(2.0*vis)));
					gr1->SetPoint(i,u,y);		
			}


		
}//end of for loop	
	


    gr1->Draw("AP");
	//gr1->GetYaxis()->SetRangeUser(0,yaxis);
    gr5->Draw("samep");

	//gr2->GetYaxis()->SetRangeUser(0,yaxis);
    gr2->Draw("samep");

	//gr3->GetYaxis()->SetRangeUser(0,yaxis);
    gr3->Draw("samep");
 
	//gr4->GetYaxis()->SetRangeUser(0,yaxis);
    gr4->Draw("samep");
	

	TLegend *leg = new TLegend(0.6,0.3,0.7,0.55);

leg->AddEntry(gr1,"-27V","AP");
leg->AddEntry(gr2,"-27.5V","AP");
leg->AddEntry(gr3,"-28V","AP");
leg->AddEntry(gr4,"-28.5V","AP");
leg->AddEntry(gr5,"-29V","AP");
leg->SetFillColor(0);

leg->Draw();

	
	

		
}
