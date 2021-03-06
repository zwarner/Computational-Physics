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
//Project III
/*This program solves the damped driven pendulum using RK4 method and does the following:
(1) Calculates the first three perieod-doubling bifurcation. Plot the bifurcation in parameter space and in phase space. 
(2) calculate the frequency power spectrum for the obrit 
(3) Calculate and plot a chaoitic orbit and its frequency spectrum
*/


double pi=atan(1.0)*4.0;
double pi2=atan(1.0)*8.0;


double f1(double v, double theta,double t,double Q)
{
  double mu = 1/Q;
  double a = 1.5;
  double nu = 2.0/3.0;
 return -mu*v-sin(theta)+a*cos(nu*t);
}
	 
int main()
{

  //fourier transform
  double N=17000; //number of steps
  double t_ft_graph[40000];
  double y_graph[40000];
  double P_graph[39999], w_graph[39999];
  double t = 0;
  double A[39999],B[39999],P[39999],w[39999],f[39999];
  double temp,tempA,tempB,tempf, tempw;

  //rk4
  int i,count;
  double x0,t0,xt,a,mu,nu,v0,E0,b0,Q,v_bi,theta;
  double h;// step size RK4
  double x_RK4[5],v_RK4[5];
  double h2,t1,t2;
  double n; // number of steps for RK4
  double x_graph[30001],v_graph[30001],t_graph[30001],steps[30001],E_graph[30001];
  double v_update,x_update;
  double x[30001];
  double v[30001],vft[30001];
  double b[30001];
  double x_bi[27],Q_bi[27];
  double temp_t,temp_x,temp_v;
  double k1,k2,k3,k4,l1,l2,l3,l4; //k-values to be used in the RK4 method
 
  Q=1.38;
    x0 = 0;    // Initial angle
    t0 = 0;    // Initial time
    v0 = 0;  // Initial angular velocity
    a=1.5;
    mu=1/Q;
    nu=2.0/3.0;
    t1=0;
    t2=300;
    h=.01;
    n=(t2-t1)/h;
    b0=a*sqrt(pow(nu,2)*pow(mu,2)+pow((1-pow(nu,2)),2));
    
    h2 = 0.5*h;  // used in RK4 method
  
    
    x[0] = x0;//inital conditions for angle
    v[0] = v0;//i.c. angular velocity
    temp_t=0;

      
      for(int j=0;j<n;j++)   // RK 4th order method
      {
	temp_v=v[j];
	temp_x=x[j];
	steps[j]=j;
	
	if(temp_t>130) //eliminating transient motion
	{
	  i=j-13000;
	  vft[i]=temp_v;
	  x_graph[i]=x[j]/pi;
	  v_graph[i]=v[j]/1.5;
	 
	}
	
        k1=h*v[j];
	l1=h*f1(temp_v,temp_x,temp_t,Q);
        k2=h*(v[j]+.5*l1);
	l2=h*f1(temp_v+.5*l1,temp_x+.5*k1,temp_t+h2,Q);
        k3=h*(v[j]+.5*l2);
	l3=h*f1(temp_v+.5*l2,temp_x+.5*k2,temp_t+h2,Q);
        k4=h*(v[j]+l3);
	l4=h*f1(temp_v+l3,temp_x+k3,temp_t+h,Q);
        x[j+1]=(x[j]+(k1+2*k2+2*k3+k4)/6);
        v[j+1]=(v[j]+(l1+2*l2+2*l3+l4)/6);
        
	t_graph[j]=temp_t;
	temp_t=temp_t+h;
	
       }

	
     for(int n = 0; n<N-1; n++)//forward Fourier transform
    {
      tempA=0;
      tempB=0;
      tempw=(pi2*n)/N;

      for(int k=0; k<N-1; k++)
	{
     
	  temp = tempw*k;
	  tempA = tempA + v_graph[k]*cos(temp);
	  tempB = tempB - v_graph[k]*sin(temp);
	}
     
      A[n] = tempA;
      B[n] = tempB;
      P[n] = pow(tempA,2) + pow(tempB,2);
      P_graph[n]=log10(P[n]);
      w_graph[n] = tempw/(h*pi2);
      }
     
 //This section is for graphing

   TGraph *gr1 = new TGraph(N-1,w_graph,P_graph);
   TAxis *axis = gr1->GetXaxis();

   // gr1->Draw("AC");//Draws forward historesis loop
   axis->SetLimits(0,100); //x-axis
   gr1->GetHistogram()->SetMaximum(10);//y-axis
   gr1->GetHistogram()->SetMinimum(0);
   gr1->GetXaxis()->SetTitle("w/2pi"); 
   gr1->GetYaxis()->SetTitle("log(P)");
   gr1->GetXaxis()->CenterTitle();
   gr1->GetYaxis()->CenterTitle();
   gr1->SetTitle("Frequency Power Spectrum (Q=1.38)");
   gr1->SetLineColor(0);
   gr1->SetMarkerColor(4);
   gr1->Draw("AC*");
  
  return 0;
  
}
