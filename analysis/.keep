
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>
#include <functional>
#include <fstream>

using namespace std;

#define SNR 20       // SNR [dB]
#define M 1    // M+1 collision wihtin a slot
#define sigma_h2 1 // variance of Rayleigh fading

using namespace std;

class Slotted_ALOHA
{
public:
  vector<double> epsilon_n;
  double p;
  double q;
  double G;
  double snr;
  double theta;
  double b;
  Slotted_ALOHA(){
    for(int i=0; i<M+1; i++)
      {
        epsilon_n.push_back(0);
      }
    
      theta=1.0;
      snr = 1/(pow(10,-1*((double)SNR/10))); // SNR (linear value)
  }
  ~Slotted_ALOHA(){}
  
};

class Edge
{
public:
  int degree;
  double coef;

  Edge(int a, double b)
  {
    degree = a;
    coef = b;
  }
  ~Edge(){}
  
};

class Dist
{
public:
  vector<Edge> lambda;
  vector<Edge> rho;
  void add_Ledge(int a, double b)
  {
    lambda.push_back(Edge(a,b));
  }
  void add_Redge(int a, double b)
  {
    rho.push_back(Edge(a,b));
  }
  
  Dist()
  {}
  ~Dist(){}
};


double factorial(int n)
{
 double x=1;
    for (int i = 1; i <= n; i++) {
        x = x*i;
    }
    return x;
}


double Decryption(int t,int r,Slotted_ALOHA& a){
    double dec=1.0;

        dec=pow((1+a.theta),-t*(r-(t+1)*0.5))*exp(-1*(pow((a.theta+1),t)-1)/(sigma_h2*a.snr));
        return dec;
}

double n_Ncollision_fading_eps(double snr,Slotted_ALOHA& a,int m){
    double eps=0.0;
    double sum=0.0;
    double R=m+1;
    double D=0;
    double dec=0;
    double K=0;

    for(int t=1;t<=R;t++){
        D=Decryption(t,R,a);
        dec+=factorial(R-1)/factorial(R-t)*D;
    }
    eps=1-dec;
//    eps=0;
    return eps;
}


void set_epsilon_n(Slotted_ALOHA& a)
{
 
  for(int i=0; i<=M;i++)
   {
       a.epsilon_n[i] =n_Ncollision_fading_eps(a.snr, a, i);
   }
}

double func_lambda(Dist& d, double x)
{
  double lam = 0.0;
  for(vector<Edge>::iterator e = d.lambda.begin(); e != d.lambda.end(); e++)
    {
      lam += e->coef * pow(x, e->degree);
    }
  return lam; 
}

double func_rho(Dist& d, double G, double x)
{
  double dl = 0;
  for(vector<Edge>::iterator e = d.lambda.begin(); e != d.lambda.end();e++)
    {
      dl+=e->coef*(e->degree);
    }
    
//    cout<<dl<<endl;
  return exp(-1*G*dl*(1-x));
}

double sum_epsilon(Dist& d, double x, int m, Slotted_ALOHA& a)
{
  double sum = 0;
  double dl = 0;
  for(vector<Edge>::iterator e = d.lambda.begin(); e != d.lambda.end(); e++)
    {
      dl += e->coef * e->degree;
    }
//    cout<<dl <<endl;
  for(int i=0; i<=m;i++)
    {
      sum += ( (1-a.epsilon_n[i])/factorial(i))*pow(a.G*dl*x, i);    
    }
  return sum;
}

// node-perspective PLR
double plr(double x, Dist& Nd)
{
  double sum = 0;
  double lsum = 0;
 
  for(vector<Edge>::iterator p = Nd.lambda.begin(); p!=Nd.lambda.end(); p++)
    {
      sum += (p->coef)*pow(x, p->degree);
    }
  return sum;
}


void density_evolution_fading(Dist& Ed, Dist& Nd, int m, Slotted_ALOHA& a)
{
  double empty=0;
  double Gh = (double)(m+1);
  double Gl = 0.0;
  double q = 0;
  set_epsilon_n(a);
  cout << "Parameter: ";
    double dl = 0;
    for(vector<Edge>::iterator e = Nd.lambda.begin(); e != Nd.lambda.end();e++)
      {
        dl+=e->coef*(e->degree);
      }
  cout << "SNR " << a.snr << ", kappa= " << m+1 << ", d " << dl << endl;
  for(int i=0; i<=m;i++) cout << "epsilon " << i << " - " << a.epsilon_n[i] << endl;
  for(int i=0; ;i++)
    {
      a.G = (Gh+Gl)/2;
      a.q=1;
      a.p=1;
      q=1;
      empty = 0;
      for(int j=0; ;j++)
        {
          empty = plr(a.q, Nd);
          a.q = 1-func_rho(Nd, a.G, 1-func_lambda(Ed,a.q))*sum_epsilon(Nd, func_lambda(Ed,a.q), m, a );
          q = plr(a.q, Nd);
      
          if(fabs(q - empty) <=  pow(10, -7)) break;
        }

      if(fabs(Gh - Gl) <= pow(10, -7)) break;
       
      if(q <= pow(10, -2)) {
        Gl = a.G;
      }
      else {
        Gh = a.G;
      }
    }
  cout << "G " << a.G << endl;
}

void file_input(Dist& E_dist, Dist& N_dist, char c[])
{

  ifstream inputfile(c);
  if(!inputfile)
    {
      cout << "Failed to read the file" << endl;
    }
  cout << "Input file is " << c << endl;
  int DL;
  int DR;
  double lsum=0;
  double rsum=0;
  inputfile >> DL;
  vector<int> l_degree(DL);
  vector<double> l_coef(DL);
    
  cout << "<Edge distribution>" << endl;
  for(int i=0; i<DL; i++)
    {
      inputfile >> l_degree[i] >> l_coef[i];
      E_dist.add_Ledge(l_degree[i], l_coef[i]);
      cout << l_degree[i] << " " << l_coef[i] << endl;
      lsum+=(l_coef[i]/(l_degree[i]+1));
    }
    
  cout << "<Node distribution>" << endl;
  for(int i=0; i<DL; i++)
    {
      N_dist.add_Ledge(l_degree[i]+1, (l_coef[i]/(l_degree[i]+1))/lsum);
      cout << l_degree[i]+1 << " " << (l_coef[i]/(l_degree[i]+1))/lsum << endl;
    }
  inputfile.close();
  cout << endl;

}


int main(int argc, char* argv[])
{ 
  Dist Ed;
  Dist Nd;
  Slotted_ALOHA a;

  file_input(Ed, Nd, argv[1]);

  cout << "++++++++++AWGN cannel++++++++++" << endl;
  density_evolution_fading(Ed, Nd, 0, a);//awgn
  cout << endl;
    
  cout << "++++++++++Rayleigh Fading cannel++++++++++" << endl;
  density_evolution_fading(Ed, Nd, M, a);//fading

  return 0;
  
}
