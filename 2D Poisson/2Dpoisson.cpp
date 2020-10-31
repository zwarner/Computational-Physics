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
#include "TMatrixFBase.h"
#include "TMatrixDBase.h"
#include "THLimitsFinder.h"
#include "TMath.h"
#include "TObjString.h"
#include "TVirtualHistPainter.h"

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

//Main Program//

int main()
{

  double test_u[1000];
  double temp_x;
  double temp_y;

  /*  //  -----  PART A -----
  int N = 500;
  int M = 500;
  double h;

  double alphal = 1.0;
  double alphah = 1.0;
  double betal  = 0.0;
  double betah  = 0.0;
  double x0 = 0.0;
  double gammal[M], gammalhat[M], tempgammalhat;
  double gammah[M], gammahhat[M], tempgammahhat;
  double rho[N][M], rhohat[N][M], temprho[M];
  double w[N];
  double u[N][M] , uhat[N][M], tempuhat,tempu[M], u_graph[N];
  
  double yf = 1.0;
  double dx = 1.0/N;
  double dy = 1.0/M;
  double x[N], y[M];
  double a[N],c[N],b[N],xi[N],eta[N];
  double arg;

// Calc. constansts for tridiagonal linear eqn. that don't requiere w 

  a[1]   = 0.0;
  for(int i = 2; i < N; i++)
    a[i] = 1.0;

  c[N-1] = 0.0;
  for(int i = 1; i < N-1; i++)
    c[i] = 1.0;

// Mesh 

  for (int i = 0; i < N + 1; i++)
    x[i] = x0 + i * dx;
  for (int i = 0; i < M + 1; i++)
    y[i] = i * dy;
    
  for(int i = 0; i < M; i++)
    {
      gammal[i] = 0.0;
      gammah[i] = y[i] * (1.0 - y[i]);
    }
  
// loop for gammalhat FFT  
  for(int k = 1 ; k < M ; k++)
    gammalhat[k] = fftsin(M, k, gammal);
  
// loop for gammahhat FFt 
  for(int k = 1 ; k < M ; k++)
    gammahhat[k] = fftsin(M, k, gammah);
      
// loop for rhohat FFT 
  for(int n = 1; n < N ; n++)
   { 
     for(int m = 1; m < M; m++)
	temprho[m]=6.0*x[n]*y[m]*(1.0 - y[m]) - 2.0*pow(x[n],3);
       
     for(int k = 1; k < M ; k++)
        rhohat[n][k] = fftsin(M, k, temprho);
    }

// loop for 1D possoin solver for uhat from w 
   for(int k = 1; k < N; k++)
     {
       w[1]   = pow(dx,2)*rhohat[1][k] - (dx*gammalhat[k]) / (alphal*dx - betal);
       w[N-1] = pow(dx,2)*rhohat[N-1][k] - (dx*gammahhat[k]) / (alphah*dx -betah);
       for(int n = 2; n < N-1; n++)
	 w[n] = pow(dx,2)*rhohat[n][k];

       arg= (k*pi*dx);
       b[1]   = -( 2 + pow(arg,2) ) - betal / (alphal*dx - betal);
       b[N-1] = -( 2 + pow(arg,2) ) + betah / (alphah*dx + betah);
       for(int i = 2; i < N-1; i++)
	 b[i] = -( 2 + pow(arg,2) );

       xi[N-1]  = -a[N-1]/b[N-1];
       xi[1]    = 0;  
     // calc. remaiaing constants for tridiagonal linear eqn. 
       eta[N-1] = w[N-1]/b[N-1];
       for (int i = N-2; i > 1; i--)
	 {
	   xi[i]  = -a[i] / (b[i] + c[i] * xi[i+1]);
	   eta[i] = (w[i] - c[i] * eta[i+1]) / (b[i] + c[i] * xi[i+1]);
	 }
       eta[1]   = (w[1] - c[1] * eta[2]) / (b[1] + c[1] * xi[2]);

       uhat[0][k] = 0.0;
       for(int n = 1; n < N ; n++)
	 uhat[n][k] = xi[n] * uhat[n-1][k] + eta[n];  
       }
  
// loop for uhat to u 
 for(int n = 1; n < N ; n++)
   { 
     for(int i = 1; i < M; i++)
	tempu[i]=uhat[n][i];
	
     for(int m = 1; m < M ; m++)
       u[n][m] = (2.0/M) * fftsin(M, m, tempu);
   }

// loop for graphing
  for(int i = 1; i < N; i++)
      u_graph[i]=u[i][250];*/

  /*
// ------  PART B -----
  int N = 500;
  int M = 500;
  double h;

  double alphal = 1.0;
  double alphah = 1.0;
  double betal  = 0.0;
  double betah  = 0.0;
  double x0 = 0.0, y0 = 0.0, xf = 1.0, yf = 1.0;
  double gammal[M], gammalhat[M], tempgammalhat;
  double gammah[M], gammahhat[M], tempgammahhat;
  double rho[N][M], rhohat[N][M], temprho[M];
  double w[N];
  double u[N][M] , uhat[N][M], tempuhat,tempu[M], u_graph[N];
  
  double dx = abs(xf-x0)/N;
  double dy = abs(yf-y0)/M;
  double x[N], y[M];
  double a[N],c[N],b[N],xi[N],eta[N];
  double arg;
  double deltafunc = - 1.0 / (dx*dy);
  double E_field;

// Calc. constansts for tridiagonal linear eqn. that don't requiere w 

  a[1]   = 0.0;
  for(int i = 2; i < N; i++)
    a[i] = 1.0;

  c[N-1] = 0.0;
  for(int i = 1; i < N-1; i++)
    c[i] = 1.0;

// Mesh 

  for (int i = 0; i < N + 1; i++)
    x[i] = x0 + i * dx;
  for (int i = 0; i < M + 1; i++)
    y[i] = i * dy;
    
  for(int i = 0; i < M; i++)
    {
      gammal[i] = 0.0;
      gammah[i] = 0.0;
    }
  
// loop for gammalhat FFT  
  for(int k = 0 ; k < M ; k++)
    gammalhat[k] = fftsin(M, k, gammal);
  
// loop for gammahhat FFt 
  for(int k = 0 ; k < M ; k++)
    gammahhat[k] = fftsin(M, k, gammah);
      
// loop for rhohat FFT 
  for(int n = 1; n < N ; n++)
   { 
     for(int k = 1; k < M ; k++)
       {
	 if(n == N/2 && k == M/2) // location of wire is at x = 0.5 and y = 0.5
	   temprho[k] = deltafunc;
	 else
	   temprho[k] = 0;
       }
     for(int k = 1; k < M; k++)
        rhohat[n][k] = fftsin(M, k, temprho);
    }

// loop for 1D possoin solver for uhat from w 
   for(int k = 1; k < N; k++)
     {
       w[1]   = pow(dx,2)*rhohat[1][k] - (dx*gammalhat[k]) / (alphal*dx - betal);
       w[N-1] = pow(dx,2)*rhohat[N-1][k] - (dx*gammahhat[k]) / (alphah*dx -betah);
       for(int n = 2; n < N-1; n++)
	 w[n] = pow(dx,2)*rhohat[n][k];

       arg= (k*pi*dx);
       b[1]   = -( 2 + pow(arg,2) ) - betal / (alphal*dx - betal);
       b[N-1] = -( 2 + pow(arg,2) ) + betah / (alphah*dx + betah);
       for(int i = 2; i < N-1; i++)
	 b[i] = -( 2 + pow(arg,2) );

       xi[N-1]  = -a[N-1]/b[N-1];
       xi[1]    = 0;  
     // calc. remaiaing constants for tridiagonal linear eqn. 
       eta[N-1] = w[N-1]/b[N-1];
       for (int i = N-2; i > 1; i--)
	 {
	   xi[i]  = -a[i] / (b[i] + c[i] * xi[i+1]);
	   eta[i] = (w[i] - c[i] * eta[i+1]) / (b[i] + c[i] * xi[i+1]);
	 }
       eta[1]   = (w[1] - c[1] * eta[2]) / (b[1] + c[1] * xi[2]);

       uhat[0][k] = 0.0;
       for(int n = 1; n < N ; n++)
	 uhat[n][k] = xi[n] * uhat[n-1][k] + eta[n];  
       }
  
// loop for uhat to u 
 for(int n = 1; n < N ; n++)
   { 
     for(int i = 1; i < M; i++)
	tempu[i]=uhat[n][i];
	
     for(int m = 1; m < M ; m++)
       u[n][m] = (2.0/M) * fftsin(M, m, tempu);
   }

 //calc E-field
 E_field = -(u[251][250] - u[250][250])/dx;
 cout <<"E_field @ x = " << 251.0/N << " and y = " << 250.0/N << " : "<<  E_field << endl;
 E_field = -(u[253][250] - u[252][250])/dx;
 cout <<"E_field @ x = " << 253.0/N << " and y = " << 250.0/N << " : "<<  E_field << endl;
 E_field = -(u[255][250] - u[254][250])/dx;
 cout <<"E_field @ x = " << 255.0/N << " and y = " << 250.0/N << " : "<<  E_field << endl;*/
 
 
// ------  PART C -----
  int N = 500;
  int M = 500;
  double h;

  double alphal = 1.0;
  double alphah = 1.0;
  double betal  = 0.0;
  double betah  = 0.0;
  double x0 = 0.0, y0 = 0.0, xf = 1.0, yf = 1.0;
  double gammal[M], gammalhat[M], tempgammalhat;
  double gammah[M], gammahhat[M], tempgammahhat;
  double rho[N][M], rhohat[N][M], temprho[M];
  double w[N];
  double u[N][M] , uhat[N][M], tempuhat,tempu[M], u_graph[N];
  
  double dx = abs(xf-x0)/N;
  double dy = abs(yf-y0)/M;
  double x[N], y[M];
  double a[N],c[N],b[N],xi[N],eta[N];
  double arg;
  double deltafunc = - 1.0 / (dx*dy);

// Calc. constansts for tridiagonal linear eqn. that don't requiere w 

  a[1]   = 0.0;
  for(int i = 2; i < N; i++)
    a[i] = 1.0;

  c[N-1] = 0.0;
  for(int i = 1; i < N-1; i++)
    c[i] = 1.0;

// Mesh 

  for (int i = 0; i < N + 1; i++)
    x[i] = x0 + i * dx;
  for (int i = 0; i < M + 1; i++)
    y[i] = i * dy;
    
  for(int i = 0; i < M; i++)
    {
      gammal[i] = 0.0;
      gammah[i] = 0.0;
    }
  
// loop for gammalhat FFT  
  for(int k = 0 ; k < M ; k++)
    gammalhat[k] = fftsin(M, k, gammal);
  
// loop for gammahhat FFt 
  for(int k = 0 ; k < M ; k++)
    gammahhat[k] = fftsin(M, k, gammah);
      
// loop for rhohat FFT 
  for(int n = 1; n < N ; n++)
   { 
     for(int k = 1; k < M ; k++)
       {
	 if((n == 100 && k == 100)||(n == N/2 && k == M/2)) // location of wires are x = y = 0.5 & y = x = 0
	   temprho[k] = deltafunc;
	 else
	   temprho[k] = 0;
       }
     for(int k = 1; k < M; k++)
        rhohat[n][k] = fftsin(M, k, temprho);
    }

// loop for 1D possoin solver for uhat from w 
   for(int k = 1; k < N; k++)
     {
       w[1]   = pow(dx,2)*rhohat[1][k] - (dx*gammalhat[k]) / (alphal*dx - betal);
       w[N-1] = pow(dx,2)*rhohat[N-1][k] - (dx*gammahhat[k]) / (alphah*dx -betah);
       for(int n = 2; n < N-1; n++)
	 w[n] = pow(dx,2)*rhohat[n][k];

       arg= (k*pi*dx);
       b[1]   = -( 2 + pow(arg,2) ) - betal / (alphal*dx - betal);
       b[N-1] = -( 2 + pow(arg,2) ) + betah / (alphah*dx + betah);
       for(int i = 2; i < N-1; i++)
	 b[i] = -( 2 + pow(arg,2) );

       xi[N-1]  = -a[N-1]/b[N-1];
       xi[1]    = 0;  
     // calc. remaiaing constants for tridiagonal linear eqn. 
       eta[N-1] = w[N-1]/b[N-1];
       for (int i = N-2; i > 1; i--)
	 {
	   xi[i]  = -a[i] / (b[i] + c[i] * xi[i+1]);
	   eta[i] = (w[i] - c[i] * eta[i+1]) / (b[i] + c[i] * xi[i+1]);
	 }
       eta[1]   = (w[1] - c[1] * eta[2]) / (b[1] + c[1] * xi[2]);

       uhat[0][k] = 0.0;
       for(int n = 1; n < N ; n++)
	 uhat[n][k] = xi[n] * uhat[n-1][k] + eta[n];  
       }
  
// loop for uhat to u 
 for(int n = 1; n < N ; n++)
   { 
     for(int i = 1; i < M; i++)
	tempu[i]=uhat[n][i];
	
     for(int m = 1; m < M ; m++)
       u[n][m] = (2.0/M) * fftsin(M, m, tempu);
   }

 
// loop of annalytical solution for compariosn with solution generated from poisson solver of Part A
  temp_x=0.0;
  temp_y=0.5;
  h = 1.0/N;
  for(int i =0; i < N; i++)
    {
      test_u[i] = temp_y * (1-temp_y) * pow(temp_x,3);
      temp_x = temp_x+h;
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

  TH2D* h2D = new TH2D("h2", "u(x,y) for 2 Wires (dx = dy = 1/N where N = 200)" ,N,0,xf,M,0,yf);
  for(int i = 1; i < N; i++)
      for(int j = 1; j < M; j++)	
	  h2D->SetBinContent(i,j,u[i][j]);
    
  h2D->GetXaxis()->SetTitle("x"); 
  h2D->GetYaxis()->SetTitle("y");
  h2D->GetXaxis()->CenterTitle();
  h2D->GetYaxis()->CenterTitle();
  h2D->Draw("colz");
    
  return 0;
  
}
