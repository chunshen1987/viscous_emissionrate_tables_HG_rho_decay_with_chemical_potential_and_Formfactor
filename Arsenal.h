#ifndef ARSENAL_H
#define ARSENAL_H

#include<fstream>
#include<string>
#include<vector>
#include<cstdlib>

using namespace std;

void interpolation1D_linear(double* Tb_x, double* Tb_y, double* xval, double* yval, int nTblength, int nvalpoint);

double Simpson_sum(double* , int, double);
//double interpolation2D_bilinear(Table2D* , double , double , int);

vector< vector<double>* >* readBlockData(istream &stream_in);
void releaseBlockData(vector< vector<double>* >* data);
vector<double> stringToDoubles(string str);

void
gauss (int npts, int job, double a, double b, double xpts[], double weights[]); //generate points and weights for gaussian quadrature, code is in gauss.cpp

double interpCubicDirect(vector<double>* x, vector<double>* y, double xx, bool allowExtrapolation=false);
double interpCubicMono(vector<double>* x, vector<double>* y, double xx, bool allowExtrapolation=false);
double interpLinearDirect(vector<double>* x, vector<double>* y, double xx, bool allowExtrapolation=false);
double interpLinearMono(vector<double>* x, vector<double>* y, double xx, bool allowExtrapolation=false);
double interpNearestDirect(vector<double>* x, vector<double>* y, double xx, bool allowExtrapolation=false);
double interpNearestMono(vector<double>* x, vector<double>* y, double xx, bool allowExtrapolation=false);

inline double interpCubic4Points(double A0, double A1, double A2, double A3, double dx, double xx)
// Assuming that A0, A1, A2, A3 are the values located at 0, dx, 2*dx, and 3*dx, interpolate the value
// at xx using cubic polynomials.
{
    double deltaX = xx - dx;
    return (-A0+3.0*A1-3.0*A2+A3)/(6.0*dx*dx*dx)*deltaX*deltaX*deltaX
            + (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX
            - (2.0*A0+3.0*A1-6.0*A2+A3)/(6.0*dx)*deltaX
            + A1;
}

double invertFunc(double (*func)(double), double y, double xL, double xR, double dx, double x0, double relative_accuracy=1e-10);

double invertTableDirect_hook(double xx);
double invertTableDirect(vector<double>* x, vector<double>* y, double y0, double x0, double relative_accuracy=1e-10);

long binarySearch(vector<double>* A, double value, bool skip_out_of_range=false);

void outputFunctionerror(string function_name, string massage, double value, int level);

string toLower(string str);
string trim(string str);
double stringToDouble(string str);


#endif
