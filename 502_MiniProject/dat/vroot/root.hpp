#ifndef ROOT_INCLUDED
#define ROOT_INCLUDED

#include "TF1.h"
#include <TAxis.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TF1.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TApplication.h"
#include <TLatex.h>
#include <TImage.h>
#include <TRandom3.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TPad.h>
#include <TFrame.h>
#include "TH3.h"
#include "TNtuple.h"

#include <TRandom3.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <math.h> 
#include <thread>
#include <stdio.h>
#include <iterator>
#include <sys/stat.h>
#include <unistd.h>
#include "/home/jerin/cfd510/Poisson/numerical/numerical.hpp"

using namespace std;


#define PI 3.141592654


struct dat{
float d1 = 0;
float d2 = 0;
float d3 = 0;
float d4 = 0;
float d5 = 0;
float d6 = 0;
float d7 = 0;
float d8 = 0;
float d9 = 0;
float d10 = 0;
float d11 = 0;
float d12 = 0;

};

struct datset{
vector<dat> datvec;

float d0 = 0;
float d1 = 0;

float l0 = 0;
float l1 = 0;
};


void plot_pics(datset myset);

void get_data(datset & myset, string cadfile, int dsize);

void drawt_data(datset & myset);

void plot_rey(datset myset);

void plotcd_cl_a(datset myset);

void plot_press(datset myset);

void get_cali(datset & myset, float & d0, float & d1, float & l0, float & l1);

int run();







void draw_graph(carray myarray, carray myarray2);

void draw_3DgraphP(carray myarray);

void draw_3Dgraph(carray myarray, carray myarray2);

void draw_graph_diff3(carray myarray, carray myarray2, carray myarray3);

void draw_graph_l2norm3(carray myarray, carray myarray2, carray myarray3);













#endif
