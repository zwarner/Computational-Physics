#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <time.h>

using namespace std;

//Zach Warner
//Proj. V homework
/*This program generates Gaussian Random numbers and simlates random walks in one dimension
*/

//generates random number between 0 and 1
double randd() {
  return (double)rand() / ((double)RAND_MAX + 1);
}

double pi = atan(1)*4;

int main()
{

// ------ Homework 1 -----
  
  int N = 1000;
  double  x1[N], x2[N], y[N];
  double y0 = 0;
  double sig = 1;
  double tempy;
  double moment, tempmoment;
  double avgy, tempavgy;
  int n = 2;
  
  //loop for Gaussian Random number generator
  for(int i = 0; i < N; i++)
    {
      x1[i] = randd();
      x2[i] = randd();
      tempy = sqrt( -2 * log(x1[i]) ) * cos( 2*pi*x2[i] );
      y[i] = (tempy * sig) + y0;
    }

  tempavgy=0;
  for(int i = 0; i < N; i++)
      tempavgy = tempavgy + y[i];
  
  avgy = tempavgy / N; 
  //cout << "avgerage = " <<  avgy << endl;

  tempmoment = 0;
  for(int i = 0; i < N; i++)
      tempmoment = tempmoment + pow((y[i]-avgy),n);
    
  moment = tempmoment / N;
  // cout << "momemnt(" << n << ") = " << moment;

// ----- Homework 2 -----

  int Np = 5000;
  double h = 1.0;
  double deltat = 1.0;

  double end[Np]; //final location of particle n       
  double temppos[Np], t[N], tempt, avgx[N], moment2[N];
  int rho[2*N + 1];
  double prob[2*N + 1],S[N];
  double tempsum, sum, temps;

  int m;
  double x;
  double D = pow(h,2) / (2 * deltat);

 for(int i = 0; i < Np; i++)      
   {
     end[i] = 0;
     temppos[i]=0;
     t[i]=0;
   }

 //stepping
 tempt = 0;
 for(int i = 0; i < N; i++)
   {

   for(int i = 0; i < 2*N ; i++)
     {
       prob[i] = 0;
       rho[i] = 0;
     }
   for(int j = 0; j < Np; j++)
     {
       m = 0;
       x = randd();
	 
       if(x<=0.5) m = 1;
       if(x>0.5)  m = -1;

       temppos[j] = temppos[j] + m ; //updates each particles position
       end[j] = temppos[j];
     }
     
    tempsum=0;
    for(int k = 0; k < Np; k++) // sums particle positions for avg.
      tempsum = tempsum + end[k];
     
    avgx[i] = tempsum / Np; // avg. position of particles
    t[i]    = i + 1;        // updates time
    moment2[i] = 2 * D * t[i]; // calculates 2nd order moment

    //calculates number of particles at each position
    for(int j = -N; j <= N; j++)
      {
        for(int k = 0; k < Np; k++)
	  if(j==end[k]) rho[j+N]++; // thus rho[0] = -N 
      }

     //calculates probability of finding particle at position
     for(int j = 0; j < (2*N+1); j++)
         prob[j] = rho[j]*1.0/Np;
        
     //calcualtes entropy at time i
     temps=0;
     for(int j = 0; j < (2*N+1); j++)
       {
	 if(prob[j]!=0) temps = temps + prob[j]*log10(prob[j]);
	 else temps = temps;
       }

     S[i] = - temps;
   }
 
 //This section is for graphing

 TGraph *gr1 = new TGraph(N,t,S);
 TAxis *axis = gr1->GetXaxis();

   gr1->Draw("AC");//Draws forward historesis loop
   axis->SetLimits(0,N); //x-axis
   gr1->GetHistogram()->SetMaximum(2);//y-axis
   gr1->GetHistogram()->SetMinimum(0);
   gr1->GetXaxis()->SetTitle("t"); 
   gr1->GetYaxis()->SetTitle("S");
   gr1->GetXaxis()->CenterTitle();
   gr1->GetYaxis()->CenterTitle();
   gr1->SetTitle("Entropy as a Finction of Time (5000 particles)");
   gr1->SetLineColor(4);
   gr1->SetMarkerColor(4);
   gr1->Draw("AC");
    
  return 0;
  
}
