#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <TGraph.h>
#include <TCanvas.h>

using namespace std;


//Zach Warner
//This program calculates the numerical solution for harmonic oscillator using both Euler's and the 4th-order Runge-Kutta method
//Plotting done using root
//angular frequence = mass = 1 thus momentum = velocity.


double pi=atan(1)*4;

int main()
{
  double xt0,t0,pt0,pt,xt,t,h,xtold,ptold;
  double x_euler[5], p_euler[5], x_RK4[5],p_RK4[5],err_euler[5],err_RK4[5];
  int i=0;
  int j=0;
  int k=0;
  double h2;
  double n; // number of steps
  double n_graph[5]; //number of steps for graphing
  double p_update,x_update;
  double x_actual=-cos(1+pi/2);
  double p_actual=sin(1+pi/2);
 
  double x[10001];
  double p[10001];
  double kx1,kx2,kx3,kx4,kp1,kp2,kp3,kp4; //k-values to be used in the RK4 method.

  cout << endl << "x expected: " << x_actual << endl;
  cout << "p expected: " << p_actual << endl << endl;
  

  for (h = 0.0001; h <= 1; h = h*10) // loop for step size
    {
       xt0 = 0;    // Initial position
       t0 = 0;     // Initial time
       pt0 = 1;    // Initial velocity (momenutm)
       h2 = h*0.5;  // used in RK4 method
  
       xtold = xt0;   // Sets initial conditions
       ptold = pt0;
       n = 1/h;
       n_graph[k]=log10(n); //log(1/h) for use in graph
       x[0] = 0;
       p[0] = 1;

       cout << "h = " << h << endl;
       cout << "n = " << n << endl;
       
        for (t=t0; t<1.0; t+=h)//Eulers method
	{
	  pt = ptold - h*xtold;
	  xt = xtold + h*ptold;
	  xtold = xt;
	  ptold = pt;
      	}

       for(j=0;j<n;j++)   // RK 4th order method
	 {
	   kx1=h*p[j];
	   kp1=-h*x[j];
	   kx2=h*(p[j]+h2*kp1);
	   kp2=-h*(x[j]+h2*kx1);
	   kx3=h*(p[j]+h2*kp2);
	   kp3=-h*(x[j]+h2*kx2);
	   kx4=h*(p[j]+h*kp3);
	   kp4=-h*(x[j]+h*kx3);
	   x[j+1]=(x[j]+(kx1+2*kx2+2*kx3+kx4)/6);
	   p[j+1]=(p[j]+(kp1+2*kp2+2*kp3+kp4)/6);
	   p_update=p[j+1];
	   x_update=x[j+1];
	 }

       x_euler[k]=xtold;
       p_euler[k]=ptold;
	   
       cout << "   Euler(x): " << x_euler[k] << endl;
       cout << "   Euler(p): " << p_euler[k] << endl;
       
       x_RK4[k] = x_update; // stores RK 4th order values calculated for each step size (h)
       p_RK4[k] = p_update; 
       cout << "   RK4(x): " << x_RK4[k] << endl;
       cout << "   RK4(p): " << p_RK4[k] << endl;
       k++;  
    }

  for(int i=0;i<5;i++) //calculates the error for 
    {
      err_euler[i]=sqrt(pow(x_actual-x_euler[i],2)+pow(p_actual-p_euler[i],2))/sqrt(pow(x_actual,2)+pow(p_actual,2));
      err_RK4[i]=sqrt(pow(x_actual-x_RK4[i],2)+pow(p_actual-p_RK4[i],2))/sqrt(pow(x_actual,2)+pow(p_actual,2));
      err_euler[i]=log10(err_euler[i]);
      err_RK4[i]=log10(err_RK4[i]);
      cout << "euler error: " << err_euler[i] << endl;
      cout << "RK4 error: " << err_RK4[i] << endl << endl;
    }

  //This section graphs log(|error|) vs log(n)

  TGraph *gr1 = new TGraph(5,n_graph,err_euler);
  TGraph *gr2 = new TGraph(5,n_graph,err_RK4);

   gr1->Draw("AC*");//Draws euler error
   gr1->GetXaxis()->SetTitle("log(1/h)");
   gr1->GetYaxis()->SetTitle("log(|error|)");
   gr1->GetXaxis()->CenterTitle();
   gr1->GetYaxis()->CenterTitle();
   gr1->Draw("AC*");
   gr2->SetLineWidth(1);
   gr2->SetLineColor(2);
   gr2->Draw("C*");//Draws RK4 Error
   
   //Euler method gave largest error since it is a first order approximation and RK is 4th order. Also the error grows smaller with increasing step size. This also makes sense since step size determins the amount of information you are considering between the ranges of t=0 to t=1
  
  return 0;
  
}
