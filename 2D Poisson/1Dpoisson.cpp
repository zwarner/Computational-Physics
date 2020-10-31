#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>

using namespace std;

//Zach Warner
//Proj. IV homework
/*This program solves a 1-Dimensional Poisson Problem with Mixed Boundry conditions for rho(x)=1-2x^2
*/
	 
int main()
{

  double test_u[1000];
  double test_x[1000];
  double temp_x;
  int N = 1000;
  double temph;
  //constants for boundry conditions
  double alpha_l=1.0;
  double alpha_h=1.0;
  double beta_l=-1.0;
  double beta_h=1.0;
  double gamma_l=1.0;
  double gamma_h=1.0;
  double h=(1.0-0.0)/N; //delta x value (i.e. grid intervals)

  //coefficiens of the tridiagonal linear equation
  double a[N]; 
  double b[N];
  double c[N];
  double w[N];
  double rho[N];

  //coefficeints of recusive relationships
  double xi[N];
  double eta[N];

  //Poisson solution
  double u[N];
  double u_graph[N];
  double x_graph[N];

  //calculating coefficeints of tridiagonal linear equation
  temph=h;
  for(int i=1; i < N; i++)
    {
      rho[i]=1.0-2.0*pow(temph,2);
      temph = temph+h;
    }
      
  a[1]=0.0;
  for(int i=2; i < N; i++)
    {
      a[i]=1.0;
    }

  b[1]=-2.0-beta_l/(alpha_l*h-beta_l);
  cout << b[1];
  for(int i=2; i < N-1; i++)
    {
      b[i]=-2.0;
    }
  b[N-1]=-2+beta_h/(alpha_h*h+beta_h);
  
  c[N-1]=0.0;
  for(int i=1; i < N-1; i++)
    {
      c[i]=1.0;
    }

  w[1]=rho[1]*pow(h,2)-(gamma_l*h)/(alpha_l*h-beta_l);
  for(int i=2;i<N-1;i++)
    {
      w[i]= rho[i]*pow(h,2);
    }
  w[N-1]=rho[N-1]*pow(h,2)-(gamma_h*h)/(alpha_h*h+beta_h);

  //calculating coefficients for recursive relationship
  xi[N-1]=-a[N-1]/b[N-1];
  eta[N-1]=w[N-1]/b[N-1];
  
  for(int i=N-2; i>0; i--)
    {
      xi[i]=-a[i]/(b[i]+c[i]*xi[i+1]);
      eta[i]=(w[i]-c[i]*eta[i+1])/(b[i]+c[i]*xi[i+1]);
    }
  
  //caluclating solution to poisson equation
  u[0]=7.0/9.0;
  for(int i=1;i<N;i++)
    {
      u[i]=xi[i]*u[i-1]+eta[i];
      cout << u[i] << endl;
    }


  for(int i=0;i<N-1;i++)
    {
      u_graph[i]=u[i+1];
    }

  
  //loop of annalytical solution for compariosn with solution generated from poison solver
  temp_x=0.0;
  for(int i =0; i<N;i++)
    {
      test_x[i]=temp_x;
      test_u[i]=-(1.0/6)*pow(test_x[i],4)+(1.0/2)*pow(test_x[i],2)-(2.0/9)*test_x[i]+(7.0/9);
	temp_x=temp_x+h;
    }


  
 //This section is for graphing

   TGraph *gr1 = new TGraph(N-1,test_x,test_u);
   TAxis *axis = gr1->GetXaxis();

   gr1->Draw("AC");//Draws forward historesis loop
   axis->SetLimits(0,1); //x-axis
   gr1->GetHistogram()->SetMaximum(1);//y-axis
   gr1->GetHistogram()->SetMinimum(0.5);
   gr1->GetXaxis()->SetTitle("x"); 
   gr1->GetYaxis()->SetTitle("u(x)");
   gr1->GetXaxis()->CenterTitle();
   gr1->GetYaxis()->CenterTitle();
   gr1->SetTitle("Solution to Poisson Equation (Analytical)");
   gr1->SetLineColor(4);
   gr1->SetMarkerColor(4);
   gr1->Draw("AC");
    

  return 0;
  
}
