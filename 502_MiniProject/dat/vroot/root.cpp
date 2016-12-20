#include "root.hpp"

using namespace std;



int run() 
{ 
TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);



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
	
			
			//Data filling *********************************
float h = 8.0*sqrt(vis/omega);
	

				
			for (int i = 0; y < h; ++i)
			{		omega = 1.0;
                   
                    y = i*(0.001)*sqrt(1.0/(2.0*vis));
                    t = 0.0;
					u = exp(-y*sqrt(omega/(2.0*vis)))*cos(omega*t-y*sqrt(omega/(2.0*vis)));
					gr1->SetPoint(i,u,y);		
		
				omega = 5.0;
                    t = (0.0);
					u = exp(-y*sqrt(omega/(2.0*vis)))*cos(omega*t-y*sqrt(omega/(2.0*vis)));
					gr2->SetPoint(i,u,y);		
		
		omega = 20.0;
                    t = (0.0);
					u = exp(-y*sqrt(omega/(2.0*vis)))*cos(omega*t-y*sqrt(omega/(2.0*vis)));
					gr3->SetPoint(i,u,y);		
	omega = 40.0;
	
                    t = (0);
					u = exp(-y*sqrt(omega/(2.0*vis)))*cos(omega*t-y*sqrt(omega/(2.0*vis)));
					gr4->SetPoint(i,u,y);		
	
omega = 80.0;

                    t = (0);
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

leg->AddEntry(gr1,"w = 1","AP");
leg->AddEntry(gr2,"w = 5","AP");
leg->AddEntry(gr3,"w = 20","AP");
leg->AddEntry(gr4,"w = 40","AP");
leg->AddEntry(gr5,"w = 80","AP");

leg->SetFillColor(0);

leg->Draw(); 
//---------------------------------------------2d graph --------------------------------//

TCanvas *c2 = new TCanvas("c2","The FillRandom example",200,50,900,700);
TCanvas *c4 = new TCanvas("c4","The FillRandom example",200,50,900,700);

TGraph2D *gr11 = new TGraph2D();
TGraph2D *gr55 = new TGraph2D();

string titlefile = "Energy Explicit flow1; t; y; u";
c = titlefile.c_str();	
gr11->SetTitle(c);




titlefile = "Energy Implicit flow2; t; y; u";
c = titlefile.c_str();
gr55->SetTitle(c);	


			


int N = 0;
t = 0.0;
y = 0.0;
for (int j = 0; j < 100; j++) 
{t = t + 0.1;
	
	for (int i = 0; i < 70; ++i)
	{	
		
      y = i*(0.1)*sqrt(omega/(2.0*vis));               					
	  u = exp(-y*sqrt(omega/(2.0*vis)))*cos(omega*t-y*sqrt(omega/(2.0*vis)));

      gr11->SetPoint(N,t,y,u);
      gr55->SetPoint(N,t,y,u);
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





//--------------------get calibration curve----------------------//

void get_cali(datset & myset, float & d0, float & d1, float & l0, float & l1)
{
TCanvas *c3 = new TCanvas("c3","The FillRandom example",200,10,900,700);
c3->cd();
TGraph *g = new TGraph();
TGraph *g1 = new TGraph();
TGraph *g2 = new TGraph();

for(int i = 0; i < myset.datvec.size(); ++i)
{ 
dat temp = myset.datvec.at(i);

g1->SetPoint(i,temp.d2, temp.d1);
g2->SetPoint(i,temp.d3, temp.d1);

}

g->SetPoint(0,-3, 100);
g->SetPoint(1,3, 0);



TF1 *tfit = new TF1("tfit", "pol1");
tfit->SetLineColor(2);

g1->Fit(tfit);

TF1 *tfit2 = new TF1("tfit2", "pol1");
tfit2->SetLineColor(4);
g2->Fit(tfit2);

d0 = tfit->GetParameter(0);
d1 = tfit->GetParameter(1);

l0 = tfit2->GetParameter(0);
l1 = tfit2->GetParameter(1);





printf("lift equation: Load = %f + %f   drag eqauttion: Load = %f + %f\n",l1,l0,d1,d0);


g1->SetMarkerColor(9);
g1->SetMarkerStyle(22);

g2->SetMarkerColor(9);
g2->SetMarkerStyle(21);
gStyle->SetOptFit();

g->SetTitle("Drag Coefficient and Lift Coefficient Calibration; Voltage (V); Load (N)");
g->Draw("AP");
g1->Draw("SAMEP");
g2->Draw("SAMEP");

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);

leg1->AddEntry(g1,"Drag Calibration","AP");
leg1->AddEntry(g2,"Lift Calibration","AP");
leg1->AddEntry(tfit,"Drag Linear fit","l");
leg1->AddEntry(tfit2,"Lift Linear fit","l");

//leg->AddEntry(fitb,"this one","l");
leg1->Draw();


}



void plot_pics(datset myset)
{
TCanvas *c15 = new TCanvas("c15","The FillRandom example",200,50,900,700);
string titlefile = "Poisson Analytical; AOA (Degrees); Channel 0-35; Gauge Pressure (KPa)";
const char* c; 
c = titlefile.c_str();	
TGraph2D *gr = new TGraph2D();





int N =0;

for(int j = 0; j < myset.datvec.size(); ++j)
{int p = -4; 
     
   

     
    dat temp = myset.datvec.at(j); 
    
   
      gr->SetPoint(N,p,j,((789.0)*-9.81*sin(PI/3.0)*(temp.d1-8.85)*0.0254)/1000);
      ++N;
      p = p+2;
      gr->SetPoint(N,p,j,((789.0)*-9.81*sin(PI/3.0)*(temp.d2-8.85)*0.0254)/1000);
      ++N;
      p = p+2;
      gr->SetPoint(N,p,j,((789.0)*-9.81*sin(PI/3.0)*(temp.d3-8.85)*0.0254)/1000);
      ++N;
      p = p+2;
      gr->SetPoint(N,p,j,((789.0)*-9.81*sin(PI/3.0)*(temp.d4-8.85)*0.0254)/1000);
      ++N;
      p = p+2;
      gr->SetPoint(N,p,j,((789.0)*-9.81*sin(PI/3.0)*(temp.d5-8.85)*0.0254)/1000);
      ++N;
      p = p+2;
      gr->SetPoint(N,p,j,((789.0)*-9.81*sin(PI/3.0)*(temp.d6-8.85)*0.0254)/1000);
      ++N;
      p = p+2;
      gr->SetPoint(N,p,j,((789.0)*-9.81*sin(PI/3.0)*(temp.d7-8.85)*0.0254)/1000);
      ++N;
      p = p+2;
      gr->SetPoint(N,p,j,((789.0)*-9.81*sin(PI/3.0)*(temp.d8-8.85)*0.0254)/1000);
      ++N;
      p = p+2;
      gr->SetPoint(N,p,j,((789.0)*-9.81*sin(PI/3.0)*(temp.d9-8.85)*0.0254)/1000);
      ++N;
      p = p+2;
      gr->SetPoint(N,p,j,((789.0)*-9.81*sin(PI/3.0)*(temp.d10-8.85)*0.0254)/1000);
      ++N;
      p = p+2;
      gr->SetPoint(N,p,j,((789.0)*-9.81*sin(PI/3.0)*(temp.d11-8.85)*0.0254)/1000);
      ++N;
      p = p+2;
      gr->SetPoint(N,p,j,((789.0)*-9.81*sin(PI/3.0)*(temp.d12-8.85)*0.0254)/1000);
      ++N;
      p = p+2;

      
}





gr->SetTitle(c);

titlefile = "Poisson Analytical; AOA (Degrees); Channel 0-35; z";
c = titlefile.c_str();
	

   gStyle->SetPalette(1);
   //gr->Draw("colz");
   gr->Draw("colz"); 

   

   
}



//--------------------------------------reynolds---------------------------------------------///
void plot_rey(datset myset)
{

TCanvas *c24 = new TCanvas("c24","The FillRandom example",200,10,900,700);
c24->cd();
TGraph *g = new TGraph();
TGraphErrors *g1 = new TGraphErrors();
TGraphErrors *g2 = new TGraphErrors();

int n= 0;
for(int i = 0; i < myset.datvec.size(); ++i)
{ 
dat temp = myset.datvec.at(i);

float drag = temp.d2*myset.d1+myset.d0;
float cd = (2.0*drag)/(1.2041*pow(4.03*sqrt(30.0),2)*0.3075*0.6858);

float reynolds = (4.03*sqrt(temp.d1)*0.3075)/0.00001511;


g1->SetPoint(i,reynolds,cd);


}
g->SetPoint(0,1000, 1);
g->SetPoint(1,1000000, 0);

/*TF1 *tfit = new TF1("myfit","([0]*(1/((x-[1]))))+[2]", 5, 17);
tfit->SetParameter(0, 1);
tfit->SetParameter(1, 35);
tfit->SetParLimits(1, 18.75, 40);
tfit->SetParameter(2, 0);

tfit->SetLineColor(2);
tfit->SetLineWidth(1);
g1->Fit(tfit, "", "", -7, 18);*/

TF1 *tfit2 = new TF1("tfit2", "pol3");
tfit2->SetLineColor(4);
tfit2->SetLineWidth(1);
g1->Fit(tfit2, "", "", 200000, 900000);

	g1->SetMarkerColor(46);
	g1->SetMarkerStyle(25);

	g2->SetMarkerColor(9);
	g2->SetMarkerStyle(26);

	
//gStyle->SetOptFit();

g->SetTitle("Drag Coefficient vs Reynolds; Reynolds Number;                    Cd");
g->Draw("ap");
g1->Draw("samep");


TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(g1,"Drag Coefficient","AP");

//leg2->AddEntry(tfit,"Drag Poly fit","l");
leg2->AddEntry(tfit2,"Polynomial fit","l");


//leg->AddEntry(fitb,"this one","l");
leg2->Draw();

}






//--------------------------------------pressure cp---------------------------------------------///
void plot_press(datset myset)
{

TCanvas *c23 = new TCanvas("c23","The FillRandom example",200,10,900,700);
c23->cd();
TGraph *g = new TGraph();
TGraphErrors *g1 = new TGraphErrors();
TGraphErrors *g2 = new TGraphErrors();

int n= 0;
for(int i = 0; i < myset.datvec.size(); ++i)
{ 
dat temp = myset.datvec.at(i);

//----------pressure-----------//

float pss = ((789.0)*-9.81*sin(PI/3.0)*(temp.d2*0.0254));
float pss1 = -2*(pss)/((1.2041*pow(4.03*sqrt(30.0),2)));


float x = (temp.d1/100.0)*0.3075;

if(i < 24){
g1->SetPoint(i,x,pss1);
//g1->SetPointError(i,0, dragerror); 
}

if(i>=24){
g2->SetPoint(n,x, pss1);
//g2->SetPointError(i,0, lifterror); 
++n;
}

}

g->SetPoint(0,-0.01, 20);
g->SetPoint(1,0.3075, -5);

/*TF1 *tfit = new TF1("myfit","([0]*(1/((x-[1]))))+[2]", 5, 17);
tfit->SetParameter(0, 1);
tfit->SetParameter(1, 35);
tfit->SetParLimits(1, 18.75, 40);
tfit->SetParameter(2, 0);

tfit->SetLineColor(2);
tfit->SetLineWidth(1);
g1->Fit(tfit, "", "", -7, 18);

TF1 *tfit2 = new TF1("tfit2", "pol6");
tfit2->SetLineColor(4);
tfit2->SetLineWidth(1);
g2->Fit(tfit2, "", "", -7, 18);*/

	g1->SetMarkerColor(46);
	g1->SetMarkerStyle(25);

	g2->SetMarkerColor(9);
	g2->SetMarkerStyle(26);

	
//gStyle->SetOptFit();

g->SetTitle("Drag  and Lift Coefficient vs AOA; chord x (m); -Cp");
g->Draw("ap");
g2->Draw("samePc");
//g1->Draw("SAMEP");
g1->Draw("SAMEPc");


TLegend *leg2 = new TLegend(0.7282851,0.1008902,0.8997773,0.2299703);

leg2->AddEntry(g1,"Suction Surface","APl");
leg2->AddEntry(g2,"Pressure Surface","APl");
//leg2->AddEntry(tfit,"Drag Poly fit","l");
//leg2->AddEntry(tfit2,"Lift Poly fit","l");


//leg->AddEntry(fitb,"this one","l");
leg2->Draw();


//integration
float tarea = 0;

for(int i = 0; i <23; ++i)
{ 
dat temp1 = myset.datvec.at(i);
dat temp2 = myset.datvec.at(i+1);

//----------pressure-----------//

float ps1 = ((789.0)*-9.81*sin(PI/3.0)*(temp1.d2*0.0254));
float ps2 = ((789.0)*-9.81*sin(PI/3.0)*(temp2.d2*0.0254));

float pss1 = -(ps1)/(0.5*(1.2041*pow(4.03,2)*30.0));
float pss2 = -(ps2)/(0.5*(1.2041*pow(4.03,2)*30.0));

float chx = (temp2.d1/100.0)*0.3075-(temp1.d1/100.0)*0.3075;
printf("pss2 %f pss1 %f chx %f  x2 %f x1 %f\n", pss2, pss1, chx, temp2.d1, temp1.d1);

float area = ((pss2+pss1)/2.0)*chx;
tarea = area+tarea;

}

cout << "tarea: " << tarea << "\n";

for(int i = 24; i < myset.datvec.size()-1; ++i)
{ 
dat temp1 = myset.datvec.at(i);
dat temp2 = myset.datvec.at(i+1);

//----------pressure-----------//
float ps1 = ((789.0)*-9.81*sin(PI/3.0)*(temp1.d2*0.0254));
float ps2 = ((789.0)*-9.81*sin(PI/3.0)*(temp2.d2*0.0254));

float pss1 = -(ps1-296.0)/(0.5*(1.2041*pow(4.03,2)*30.0));
float pss2 = -(ps2-296.0)/(0.5*(1.2041*pow(4.03,2)*30.0));
float chx = (temp2.d1/100.0)*0.3075-(temp1.d1/100.0)*0.3075;
printf("pss2 %f pss1 %f chx %f  x2 %f x1 %f\n", pss2, pss1, chx, temp2.d1, temp1.d1);

float area = ((pss2+pss1)/2.0)*chx;
tarea = area+tarea;

}
cout << "tarea: " << tarea << "\n";
}
//---------------------------cl vs cd vs alpha-----------------------//


void plotcd_cl_a(datset myset)
{

TCanvas *c22 = new TCanvas("c22","The FillRandom example",250,10,900,700);
c22->cd();
TGraph *g = new TGraph();
TGraphErrors *g1 = new TGraphErrors();
TGraphErrors *g2 = new TGraphErrors();
TGraphErrors *g3 = new TGraphErrors();
TGraphErrors *g4 = new TGraphErrors();

for(int i = 0; i < myset.datvec.size(); ++i)
{ 
dat temp = myset.datvec.at(i);

//-----------drag------------//

float drag = temp.d2*myset.d1+myset.d0;
float cd = (2.0*drag)/(1.2041*pow(4.03*sqrt(30.0),2)*0.3075*0.6858);

float dragerror = 10*cd*sqrt(pow((temp.d3*myset.d1+myset.d0)/(temp.d2*myset.d1+myset.d0),2)+0.000024124); 
cout << "\nerror: " << dragerror;
g1->SetPoint(i,temp.d1,cd);
g1->SetPointError(i,0, dragerror); 

//------------lift----------//

float lift = temp.d4*myset.d1+myset.d0;
float cl = (2.0*lift)/(1.2041*pow(4.03*sqrt(30.0),2)*0.3075*0.6858);

float lifterror = 10*cl*sqrt(pow((temp.d5*myset.l1+myset.l0)/(temp.d4*myset.d1+myset.d0),2)+0.000024124);

g2->SetPoint(i,temp.d1, cl);
g2->SetPointError(i,0, lifterror); 



//------------theory--------------//
float tcl = 2*PI*(1.0+(0.77*0.111))*(((temp.d1/360)*2*PI)+(2.0*0.024))*(0.6858);

g3->SetPoint(i,temp.d1, tcl);


}

g->SetPoint(0,-3, 100);
g->SetPoint(1,3, 0);
g4->SetPoint(0,6, 0.388642);

TF1 *tfit = new TF1("myfit","([0]*(1/((x-[1]))))+[2]", 5, 17);
tfit->SetParameter(0, 1);
tfit->SetParameter(1, 35);
tfit->SetParLimits(1, 18.75, 40);
tfit->SetParameter(2, 0);

tfit->SetLineColor(2);
tfit->SetLineWidth(2);
tfit->SetLineStyle(9);
g1->Fit(tfit, "", "", -7, 18);

TF1 *tfit2 = new TF1("tfit2", "pol6");
tfit2->SetLineColor(4);
tfit2->SetLineWidth(2);
tfit2->SetLineStyle(10);
g2->Fit(tfit2, "", "", -7, 18);

	g1->SetMarkerColor(1);
	g1->SetMarkerStyle(21);
	g1->SetMarkerSize(1.0);

	g2->SetMarkerColor(1);
	g2->SetMarkerStyle(22);
	g2->SetMarkerSize(1.25);
	
	g4->SetMarkerColor(1);
	g4->SetMarkerStyle(33);
	g4->SetMarkerSize(2);
	
	g3->SetLineColor(8);
	g3->SetLineWidth(2);
	
//gStyle->SetOptFit();

g2->SetTitle("Drag  and Lift Coefficient vs AOA; AOA (Degrees);");
g2->Draw("AP");
//g1->Draw("SAMEP");
g1->Draw("SAMEP");
g3->Draw("SAMEPc");
g4->Draw("SAMEP");

TLegend *leg2 = new TLegend(0.7282851,0.1008902,0.8997773,0.2299703);

leg2->AddEntry(g1,"Drag Coefficient","APl");
leg2->AddEntry(g2,"Lift Coefficient","APl");
leg2->AddEntry(tfit,"Drag fit","l");
leg2->AddEntry(tfit2,"Lift fit","l");
leg2->AddEntry(g3,"Theoretical Lift Coefficient","l");
leg2->AddEntry(g4,"Integrated Pressure","P");


//leg->AddEntry(fitb,"this one","l");
leg2->Draw();


}




//---------------------getdata----------------------------------//
void drawt_data(datset & myset)
{
 for(int i = 0; i < myset.datvec.size(); ++i)
    {
    dat temp;
    temp = myset.datvec.at(i);
    printf("writing d1: %f d2: %f d3: %f d4: %f d5: %f d6: %f d7: %f d8: %f d9: %f d10: %f d11: %f d12: %f \n", temp.d1, temp.d2, temp.d3, temp.d4, temp.d5, temp.d6 , temp.d7, temp.d8, temp.d9, temp.d10, temp.d11, temp.d12); 
    }


}


void get_data(datset & myset, string cadfile, int dsize)
{  
cout << "loading file: " << cadfile;
const char * path;

char norm;

        //enter .txt and .wfi files to be used
	string wfifile= cadfile;
	path = wfifile.c_str();

dat temp;
//find vertex

std::cout << "\nRegular .txt Import: ";
	
		FILE * file = fopen(path, "r");
		if( file == NULL ){
			printf("File not found with path given\n");
			getchar();
			return;
		}

		int ind = 0;
		int ggg = 0;

		while( 1 ){
			
			char lineHeader[500];
			// read the first word of the line
			int res = fscanf(file, "%s", lineHeader);
			if (res == EOF)
				break; // EOF = End Of File. Quit the loop.

			
			if ( strcmp( lineHeader, "data" ) == 0 ){
				
				int matches = fscanf(file, "%f %f %f %f %f %f %f %f %f %f %f %f \n", &temp.d1, &temp.d2, &temp.d3, &temp.d4, &temp.d5, &temp.d6 , &temp.d7, &temp.d8, &temp.d9, &temp.d10, &temp.d11, &temp.d12);
			
				if (matches != dsize){
					printf("File can't be read by our simple parser :-( problem with vertices\n");
					return;
				}				
				myset.datvec.push_back(temp);
								
			}
		

		}
      
cout << "  size:" << myset.datvec.size() << "\nFinshed \n\n";		
}//end of getData()















//--------------------------3D graphs----------------------------//

void draw_3DgraphP(carray myarray)
{
float DIM2 = myarray.DIM1;

TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);



string titlefile = "Poisson Numerical; x; y; z";
const char* c; 
c = titlefile.c_str();	

TGraph2D *gr = new TGraph2D();


gr->SetTitle(c);

titlefile = "Poisson Analytical; x; y; z";
c = titlefile.c_str();
	


int N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIM2*(i-0.5);
	  float dy = DIM2*(j-0.5);
	  float T = myarray.mcell[i][j];

	      gr->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}


N = 0;

c1->cd();

   gStyle->SetPalette(1);
   //gr->Draw("colz");
   gr->Draw("surf1z"); 

   

   
}

//--------------------------3D graphs----------------------------//

void draw_3Dgraph(carray myarray, carray myarray2)
{
float DIM2 = myarray.DIM1;

TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);
TCanvas *c2 = new TCanvas("c2","The FillRandom example",200,50,900,700);
TCanvas *c4 = new TCanvas("c4","The FillRandom example",200,50,900,700);

c1->cd();
string titlefile = "Laplace Numerical; x; y; z";
const char* c; 
c = titlefile.c_str();	

TGraph2D *gr = new TGraph2D();
TGraph2D *gr1 = new TGraph2D();
TGraph2D *gr5 = new TGraph2D();
gr->SetTitle(c);

titlefile = "Laplace Analytical; x; y; z";
c = titlefile.c_str();
gr1->SetTitle(c);	

int N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIM2*(i-1);
	  float dy = DIM2*(j-1);
	      float T = (cos(PI*dx)*sinh(PI*dy))/sinh(PI);

	      gr1->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}

N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 0; j < myarray.sizey; j++) 
	{ float dx = DIM2*(i-1);
	  float dy = DIM2*(j-1);
	  float T = myarray.mcell[i][j];

	      gr->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}


N = 0;
for (int i = 1; i < myarray2.sizex-1; i++) 
{
	
	for (int j = 0; j < myarray2.sizey; j++) 
	{ float dx = DIM2*(i-1);
	  float dy = DIM2*(j-1);
	  float T = myarray2.mcell[i][j];

	      gr5->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}
c1->cd();

   gStyle->SetPalette(1);
   gr->Draw("colz");
c2->cd();   
   gStyle->SetPalette(1);
   gr1->Draw("colz");
   
c4->cd();   
   gStyle->SetPalette(1);
   gr5->Draw("colz");   
   
}


//--------------------------diff graphs----------------------------//

void draw_graph(carray myarray, carray myarray2)
{TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);  
 TCanvas *c2 = new TCanvas("c2","The FillRandom example",200,50,900,700); 
string titlefile;
const char* c;  
float DIM2 = myarray.DIM1;

c1->cd();
titlefile = "Convergence Behavior for a 10 x 10 mesh; Iterations (N); T(K)-T(K+1)";
c = titlefile.c_str();
TGraph *gr1 = new TGraph();	
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(7);
	gr1->SetTitle(c);  
TGraph *gr2 = new TGraph();	
	gr2->SetMarkerColor(3);
	gr2->SetMarkerStyle(7);
	gr2->SetTitle(c); 	
	
for (int i = 0; i < myarray.diff.size(); ++i)
	{	
	float tomp = myarray.diff.at(i);
	gr1->SetPoint(i,i,tomp);		
	}
for (int i = 0; i < myarray2.diff.size(); ++i)
	{	
	float temp = myarray2.diff.at(i);
	gr2->SetPoint(i,i,temp);		
	}

 c1->SetLogy();	
gr1->Draw("ALP");
gr2->Draw("sameLP");
	 

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);

leg1->AddEntry(gr1,"Gauss-Seidel","AP");
leg1->AddEntry(gr2,"Gauss-Seidel SOR 1.5","AP");

//leg->AddEntry(fitb,"this one","l");
leg1->Draw();

c2->cd();
titlefile = "L2 norm Converged Solution for 10 x 10 mesh; Iterations (N); L2";
c = titlefile.c_str();
TGraph *gr3 = new TGraph();	
	gr3->SetMarkerColor(4);
	gr3->SetMarkerStyle(7);
	gr3->SetTitle(c);  
TGraph *gr4 = new TGraph();	
	gr4->SetMarkerColor(3);
	gr4->SetMarkerStyle(7);
	gr4->SetTitle(c); 	
	
for (int i = 0; i < myarray.l2norm.size(); ++i)
	{	
	float tomp = myarray.l2norm.at(i);
	gr3->SetPoint(i,i,tomp);		
	}
for (int i = 0; i < myarray2.l2norm.size(); ++i)
	{	
	float temp = myarray2.l2norm.at(i);
	gr4->SetPoint(i,i,temp);		
	}


gr3->Draw("ALP");
gr4->Draw("sameLP");
	 

TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(gr3,"Gauss-Seidel","AP");
leg2->AddEntry(gr4,"Gauss-Seidel SOR 1.5","AP");

//leg->AddEntry(fitb,"this one","l");
leg2->Draw();


}


//--------------------------diff graphs----------------------------//

void draw_graph_diff3(carray myarray, carray myarray2, carray myarray3)
{TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);   
string titlefile;
const char* c;  
float DIM2 = myarray.DIM1;

c1->cd();
titlefile = "Convergence Behavior for a 20 x 20 mesh; Iterations (N); T(K)-T(K+1)";
c = titlefile.c_str();
TGraph *gr1 = new TGraph();	
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(7);
	gr1->SetLineColor(1);
	gr1->SetTitle(c);  
TGraph *gr2 = new TGraph();	
	gr2->SetMarkerColor(3);
	gr2->SetMarkerStyle(7);
	gr2->SetLineColor(1);
	gr2->SetTitle(c); 
TGraph *gr3 = new TGraph();	
	gr3->SetMarkerColor(2);
	gr3->SetMarkerStyle(7);
	gr3->SetLineColor(1);
	gr3->SetTitle(c); 		
	
for (int i = 0; i < myarray.diff.size(); ++i)
	{	
	float tomp = myarray.diff.at(i);
	gr1->SetPoint(i,i,tomp);		
	}
for (int i = 0; i < myarray2.diff.size(); ++i)
	{	
	float temp = myarray2.diff.at(i);
	gr2->SetPoint(i,i,temp);		
	}
for (int i = 0; i < myarray3.diff.size(); ++i)
	{	
	float temp = myarray3.diff.at(i);
	gr3->SetPoint(i,i,temp);		
	}
 
 c1->SetLogy();	

gr1->Draw("AP");
gr2->Draw("sameP");
gr3->Draw("sameP");
	 

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);

leg1->AddEntry(gr1,"SOR w=1","AP");
leg1->AddEntry(gr2,"SOR w=1.3","AP");
leg1->AddEntry(gr3,"SOR w=1.5","AP");

//leg->AddEntry(fitb,"this one","l");
leg1->Draw();



}

//--------------------------l2 graphs----------------------------//

void draw_graph_l2norm3(carray myarray, carray myarray2, carray myarray3)
{TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);   
string titlefile;
const char* c;  
float DIM2 = myarray.DIM1;

c1->cd();
titlefile = "Accuracy Behavior for a 40 x 40 mesh; Iterations (N); L2 norm                 ";
c = titlefile.c_str();
TGraph *gr1 = new TGraph();	
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(7);
	gr1->SetLineColor(1);
	gr1->SetTitle(c);  
TGraph *gr2 = new TGraph();	
	gr2->SetMarkerColor(3);
	gr2->SetMarkerStyle(7);
	gr2->SetLineColor(1);
	gr2->SetTitle(c); 
TGraph *gr3 = new TGraph();	
	gr3->SetMarkerColor(2);
	gr3->SetMarkerStyle(7);
	gr3->SetLineColor(1);
	gr3->SetTitle(c); 		
	
for (int i = 0; i < myarray.l2norm.size(); ++i)
	{	
	float tomp = myarray.l2norm.at(i);
	gr1->SetPoint(i,i,tomp);		
	}
for (int i = 0; i < myarray2.l2norm.size(); ++i)
	{	
	float temp = myarray2.l2norm.at(i);
	gr2->SetPoint(i,i,temp);		
	}
for (int i = 0; i < myarray3.l2norm.size(); ++i)
	{	
	float temp = myarray3.l2norm.at(i);
	gr3->SetPoint(i,i,temp);		
	}
 
 //c1->SetLogy();	

gr1->Draw("ACP");
gr2->Draw("sameCP");
gr3->Draw("sameCP");
	 

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);

leg1->AddEntry(gr1,"SOR w=1","AP");
leg1->AddEntry(gr2,"SOR w=1.3","AP");
leg1->AddEntry(gr3,"SOR w=1.5","AP");

//leg->AddEntry(fitb,"this one","l");
leg1->Draw();



}






