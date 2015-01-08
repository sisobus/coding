#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <cstdlib>
using namespace std;

#define MAX_DATA_SIZE   100000
//#define MAX_DIM         10
#ifndef MAX_DIM
    int MAX_DIM;
#endif
#define SPACING_VARIANCE_TH         0.0001
#define CORRELATION_COEFFICIENT_TH  0.6
double **data;
double **vps;

double dabs(double a){return a>0?a:-a;}

void initializing();
void readData(const string filename);
void readVpData(const string filename);
double squareDouble(double a);
double calculateDistance(double *p1,double *p2);
double calculateDistance(vector<double> p1,double *p2);
double calculateSpacingVariance(double *vp);
double calculateSpacingVariance(vector<double> vp);
double calculateCorrelationCoefficient(double *vp1,double *vp2);
double calculateCorrelationCoefficient(vector<double> vp1,vector<double> vp2);
vector<vector<double> > makeVeltkampVantagePoints();
void printVeltkampData(const string filename,vector<vector<double> > vp);
