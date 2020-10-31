#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include "TROOT.h"
#include "TClass.h"
#include "THashList.h"
#include "TH2.h"
#include "TVirtualPad.h"
#include "TF2.h"
#include "TProfile.h"
#include "TMath.h"

using namespace std;

//Zach Warner
//Proj. IV homework
//This program solves a 2-Dimensional Poisson Problem with Mixed Boundry conditions 


double pi = atan(1)*4;

// forward Fourier Transform
double fftsin(int M, int k, double fin[M])
{
  double tempfin;
  double fout; 
  tempfin = 0;
  for(int m = 0 ; m < M; m++)
    {
      tempfin = tempfin + fin[m] * sin( (m*k*pi) / M );
      fout = tempfin;
    }
  return fout; 
}

//generates random number between 0 and 1
double randd() {
  return (double)rand() / ((double)RAND_MAX + 1);
}

//Gaussian Random Number Generator
double gaussrand(double sig, double y0)
{
  double x1,x2,tempy;

  x1 = randd() ;
  x2 = randd() ;
  tempy = sqrt( -2 * log(x1) ) * cos(2*pi*x2);
  return (tempy * sig) + y0;
}



// ------  Main Program ------------

int main()
{

  int Np = 10000; // number of particles

  double x[Np],y[Np],px[Np],py[Np];
  double x0,y0,px0,py0, a1, b1,a2,b2,wl,wd,wr,wu;
  double sigmesh = 2;
  double sigp = 1;
  double dx = 0.2;
  double dy = 0.2;
  int sizex = 12/dx;
  int sizey = 12/dy;
  double meshx[sizex],meshy[sizey];
  int count = 0;
  int count2 = 0;
  int tempx, tempy;
  x0  = 0;
  y0  = 0;
  px0 = 0;
  py0 = 0;

  for(int i = 0; i < Np; i++)
    {
      x[i]  = gaussrand(sigmesh,x0);
      y[i]  = gaussrand(sigmesh,y0);
      //calc particles outside mesh
      if ((x[i] < -6 || x[i] > 6) || (y[i] < -6 || y[i] > 6))
	count ++;
      px[i] = gaussrand(sigp,px0);
      py[i] = gaussrand(sigp,py0);
      //cout << x[i] << "  " << y[i] << "  " << px[i] << "   " << py[i] << endl;
    }
 
   for(int i = 0; i < sizex +1; i++)
    {
      meshx[i] = i*dx - 6;
      meshy[i] = i*dy - 6;
    }
  double w[sizex][sizey], numfunc[sizex][sizey],rho[sizex][sizey];
    for(int i = 0; i < sizex; i++)
      {
	for(int j = 0; j < sizey; j++)
	  {
	    w[i][j] = 0;
	    numfunc[i][j] = 0;
	    rho[i][j] = 0;
	  }
      }
//Construct Particle Density @ grid points of mesh
  for(int i = 0; i < Np; i++)
    {
      //makes sure particle is in mesh
      if ((x[i] >= -6 && x[i] <= 6) && (y[i] >= -6 && y[i] <= 6))
	 {
	   tempx = x[i]/dx+30;
	   tempy = y[i]/dy+30;
	   // cout << tempx << "  " << tempy << endl; 
	   a1 = x[i] - tempx * dx;
	   a2 = 1.0 - a1; 
	   b1 = y[i] - tempy * dy;
	   b2 = 1.0 - b1;
	   wl = a2/dx;
	   wd = b2/dy;
	    w[tempx][tempy] = wl*wd;
	   numfunc[tempx][tempy] = numfunc[tempx][tempy] + w[tempx][tempy];
	   rho[tempx][tempy] = numfunc[tempx][tempy] / (dx*dy);
	   }
    }
 


  
// This section is for graphing 

  /*TGraph *gr1 = new TGraph(N-1,x,u_graph);
   TAxis *axis = gr1->GetXaxis();

   gr1->Draw("AC");
   axis->SetLimits(0,1); //x-axis
   gr1->GetHistogram()->SetMaximum(1.2);//y-axis
   gr1->GetHistogram()->SetMinimum(0);
   gr1->GetXaxis()->SetTitle("x"); 
   gr1->GetYaxis()->SetTitle("u(x,y");
   gr1->GetXaxis()->CenterTitle();
   gr1->GetYaxis()->CenterTitle();
   gr1->SetTitle("2-D Electrostatics Poisson Soln. for x @ y = 0.5");
   gr1->SetLineColor(4);
   gr1->SetMarkerColor(4);
   gr1->Draw("AC*");*/

  TH2D* h2D = new TH2D("h2", "rho(x,y) for PIC (sigma = 2)" ,sizex,-6,6,sizey,-6,6);
  for(int i = 1; i < sizex; i++)
      for(int j = 1; j < sizey; j++)	
	  h2D->SetBinContent(i,j,rho[i][j]);
    
  h2D->GetXaxis()->SetTitle("x"); 
  h2D->GetYaxis()->SetTitle("y");
  h2D->GetXaxis()->CenterTitle();
  h2D->GetYaxis()->CenterTitle();
  h2D->Draw("colz");
    
  return 0;
  
}
