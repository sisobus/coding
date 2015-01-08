#include "main.h"

int main(int argc,char *argv[]) {
#ifndef MAX_DIM
    int argument;
    sscanf(argv[1],"%d",&argument);
    MAX_DIM = argument;
    
    char version[111];
    strcpy(version,argv[2]);
#endif

    string data_filename = "";
    char t[1111];
    sprintf(t,"data/data_%d.txt",MAX_DIM);
    data_filename = string (t);

    string vp_filename = "";
    sprintf(t,"vp/d%d.vp",MAX_DIM);
    vp_filename = string (t);

    string veltkamp_filename = "";
    sprintf(t,"veltkamp/d%d.veltkamp",MAX_DIM);
    veltkamp_filename = string (t);

    initializing();
    readData(data_filename);

//    for ( int i = 50 ; i <= 300 ; i += 50  ){

        char s[1111];
        sprintf(s,"%s/%d.dat",version,MAX_DIM);
        string nowInputFileName = s;
        vector<vector<double> > nowPoints;
        FILE *fp = fopen(nowInputFileName.c_str(),"r");
        for ( int j = 0 ; j < MAX_DIM ; j++ ) {
            vector<double> nv;
            for ( int k = 0 ; k < MAX_DIM ; k++ ) {
                double t;
                fscanf(fp,"%lf",&t);
                nv.push_back(t);
            }
            nowPoints.push_back(nv);
        }
        puts("read ok!");
        fclose(fp);
        sprintf(s,"%s/%d.out",version,MAX_DIM);
        string nowOutputFileName = s;
        printVeltkampData(nowOutputFileName,nowPoints);
//    }
    /*
       string vp_filename = "";
       sprintf(t,"vp/d%d.vp",MAX_DIM);
       vp_filename = string (t);

       string veltkamp_filename = "";
       sprintf(t,"veltkamp/d%d.veltkamp",MAX_DIM);
       veltkamp_filename = string (t);

       initializing();
       readData(data_filename);
       readVpData(vp_filename);

#ifdef DEBUG
for ( int i = 0 ; i < MAX_DIM ; i++ ) 
for ( int j = i+1 ; j < MAX_DIM ; j++ ) 
printf("%d %d cc : %lf\n",i,j,calculateCorrelationCoefficient(vps[i],vps[j]));
#endif
vector<vector<double > > veltkamp = makeVeltkampVantagePoints();
#ifdef DEBUG
for ( int i = 0 ; i < veltkamp.size() ; i++ ) {
for ( int j = 0 ; j < veltkamp[i].size() ; j++ ) 
printf("%lf ",veltkamp[i][j]);
puts("");
}
#endif
printVeltkampData(veltkamp_filename,veltkamp);
     */

return 0;
}
void initializing() {
    data = (double **)malloc(sizeof(double*)*MAX_DATA_SIZE);
    for ( int i = 0 ; i < MAX_DATA_SIZE ; i++ ) 
        data[i] = (double *)malloc(sizeof(double)*MAX_DIM);
    vps = (double **)malloc(sizeof(double*)*MAX_DIM);
    for ( int i = 0 ; i < MAX_DIM ; i++ ) 
        vps[i] = (double *)malloc(sizeof(double)*MAX_DIM);
    puts("initializing ok!");
}
void readData(const string filename) {
    FILE *fp = fopen(filename.c_str(),"r");
    for ( int i = 0 ; i < MAX_DATA_SIZE ; i++ ) 
        for ( int j = 0 ; j < MAX_DIM ; j++ ) 
            fscanf(fp,"%lf",&data[i][j]);
    fclose(fp);
    puts("readData ok!");
}
void readVpData(const string filename) {
    FILE *fp = fopen(filename.c_str(),"r");
    for ( int i = 0 ; i < MAX_DIM ; i++ ) 
        for ( int j = 0 ; j < MAX_DIM ; j++ ) 
            fscanf(fp,"%lf",&vps[i][j]);
    fclose(fp);
    puts("readVpData ok!");
}
double squareDouble(double a) {
    return a*a;
}
double calculateDistance(double *p1,double *p2) {
    double ret = 0;
    for ( int i = 0 ; i < MAX_DIM ; i++ ) 
        ret += squareDouble(p1[i] - p2[i]);
    return sqrt(ret);
}
double calculateDistance(vector<double> p1,double *p2) {
    double ret = 0;
    for ( int i = 0 ; i < MAX_DIM ; i++ )
        ret += squareDouble(p1[i] - p2[i]);
    return sqrt(ret);
}
double calculateSpacingVariance(double *vp) {
    double *d = (double *)malloc(sizeof(double)*MAX_DATA_SIZE);
    double average = 0;
    for ( int j = 0 ; j < MAX_DATA_SIZE ; j++ ) {
        d[j] = calculateDistance(vp,data[j]);
        average += d[j];
    }
    average /= MAX_DATA_SIZE;
    sort(d,d+MAX_DATA_SIZE);
    double space = 0;
    for ( int j = 0 ; j < MAX_DATA_SIZE-1 ; j++ ) 
        space += squareDouble((d[j+1]-d[j])-average);
    space *= (1/(double)(MAX_DATA_SIZE-1));
    free(d);
    return space/MAX_DATA_SIZE;
}
double calculateSpacingVariance(vector<double> vp) {
    double *d = (double *)malloc(sizeof(double)*MAX_DATA_SIZE);
    double average = 0;
    for ( int j = 0 ; j < MAX_DATA_SIZE ; j++ ) {
        d[j] = calculateDistance(vp,data[j]);
        average += d[j];
    }
    average /= MAX_DATA_SIZE;
    sort(d,d+MAX_DATA_SIZE);
    double space = 0;
    for ( int j = 0 ; j < MAX_DATA_SIZE-1 ; j++ )
        space += squareDouble((d[j+1]-d[j])-average);
    space *= (1/(double)(MAX_DATA_SIZE-1));
    free(d);
    return space/MAX_DATA_SIZE; 
}
double calculateCorrelationCoefficient(double *vp1,double *vp2) {
    double sumOfSubD1iAndD2i = 0;
    double sumOfD1i = 0;
    double sumOfD2i = 0;
    double sumOfD1iSquare = 0;
    double sumOfD2iSquare = 0;
    for ( int i = 0 ; i < MAX_DATA_SIZE ; i++ ) {
        double d1i = calculateDistance(vp1,data[i]);
        double d2i = calculateDistance(vp2,data[i]);

        sumOfSubD1iAndD2i += (d1i*d2i);
        sumOfD1i += d1i;
        sumOfD2i += d2i;
        sumOfD1iSquare += squareDouble(d1i);
        sumOfD2iSquare += squareDouble(d2i);
    }
    double ret = 
        ( MAX_DATA_SIZE*sumOfSubD1iAndD2i - sumOfD1i * sumOfD2i )
        /   ( sqrt(MAX_DATA_SIZE*sumOfD1iSquare - squareDouble(sumOfD1i)) 
                * sqrt(MAX_DATA_SIZE*sumOfD2iSquare - squareDouble(sumOfD2i)) );
    return ret;
}
double calculateCorrelationCoefficient(vector<double> vp1,vector<double> vp2) {
    double sumOfSubD1iAndD2i = 0;
    double sumOfD1i = 0;
    double sumOfD2i = 0;
    double sumOfD1iSquare = 0;
    double sumOfD2iSquare = 0;
    for ( int i = 0 ; i < MAX_DATA_SIZE ; i++ ) {
        double d1i = calculateDistance(vp1,data[i]);
        double d2i = calculateDistance(vp2,data[i]);

        sumOfSubD1iAndD2i += (d1i*d2i);
        sumOfD1i += d1i;
        sumOfD2i += d2i;
        sumOfD1iSquare += squareDouble(d1i);
        sumOfD2iSquare += squareDouble(d2i);
    }
    double ret =
        ( MAX_DATA_SIZE*sumOfSubD1iAndD2i - sumOfD1i * sumOfD2i )
        /   ( sqrt(MAX_DATA_SIZE*sumOfD1iSquare - squareDouble(sumOfD1i))
                * sqrt(MAX_DATA_SIZE*sumOfD2iSquare - squareDouble(sumOfD2i)) );
    return ret; 
}
vector<vector<double> > makeVeltkampVantagePoints() {
    vector<vector<double> > nvps;

    srand((unsigned)time(NULL));
    bool *vis = (bool *)malloc(sizeof(bool)*MAX_DATA_SIZE);
    for ( int i = 0 ; i < MAX_DATA_SIZE ; i++ ) 
        vis[i] = false;
    for ( int i = 0 ; i < MAX_DIM ; i++ ) {
        int now;
        while ( vis[now = rand()%MAX_DATA_SIZE] );
        vis[now] = true;
        vector<double> vp;
        for ( int j = 0 ; j < MAX_DIM ; j++ ) 
            vp.push_back(data[now][j]);
        nvps.push_back(vp);
    }
    bool change = true;
    while ( change ) {
        change = false;
        vector<vector<double> >::iterator it;
        for ( it = nvps.begin() ; it != nvps.end() ; it++ ) {
            double currentSpacingVariance = calculateSpacingVariance(*it);
            if ( currentSpacingVariance > SPACING_VARIANCE_TH ) {
                nvps.erase(it);
                int now;
                while ( vis[now = rand()%MAX_DATA_SIZE] );
                vis[now] = true;
                vector<double> vp;
                for ( int j = 0 ; j < MAX_DIM ; j++ ) 
                    vp.push_back(data[now][j]);
                nvps.push_back(vp);
                it = nvps.begin();
                change = true;
            }
        }
        vector<vector<double> >::iterator it1,it2;
        for ( it1 = nvps.begin() ; it1 != nvps.end() ; it1++ ) {
            for ( it2 = nvps.begin() ; it2 != nvps.end() ; it2++ ) {
                if ( it1 == it2 ) continue;
                double currentCorrelationCoefficient = dabs(calculateCorrelationCoefficient(*it1,*it2));
                if ( currentCorrelationCoefficient > CORRELATION_COEFFICIENT_TH ) {
                    double it1CorrelationCoefficient = calculateSpacingVariance(*it1);
                    double it2CorrelationCoefficient = calculateSpacingVariance(*it2);
                    if ( it1CorrelationCoefficient > it2CorrelationCoefficient ) {
                        nvps.erase(it1);
                        int now;
                        while ( vis[now = rand()%MAX_DATA_SIZE] );
                        vis[now] = true;
                        vector<double> vp;
                        for ( int j = 0 ; j < MAX_DIM ; j++ ) 
                            vp.push_back(data[now][j]);
                        nvps.push_back(vp);
                        it1 = nvps.begin();
                        it2 = nvps.begin();
                    } else {
                        nvps.erase(it2);
                        int now;
                        while ( vis[now = rand()%MAX_DATA_SIZE] );
                        vis[now] = true;
                        vector<double> vp;
                        for ( int j = 0 ; j < MAX_DIM ; j++ ) 
                            vp.push_back(data[now][j]);
                        nvps.push_back(vp);
                        it1 = nvps.begin();
                        it2 = nvps.begin();
                    }
                    change = true;
                }
            }
        }
    }
    return nvps;
}
void printVeltkampData(const string filename,vector<vector<double> > vp) {
    FILE *fp = fopen(filename.c_str(),"w");
    double totalSpacingVariance = 0;
    double totalCorrelationCoefficient = 0;
    int cnt = 0;
    puts("printVeltkampData in");
    vector<double> correlationCoefficient;
    for ( int i = 0 ; i < (int)vp.size() ; i++ ) {
        printf("%d end\n",i);
        double currentSpacingVariance = calculateSpacingVariance(vp[i]);
        totalSpacingVariance += currentSpacingVariance;
        for ( int j = i+1 ; j < (int)vp.size() ; j++ ) {
            double currentCorrelationCoefficient = calculateCorrelationCoefficient(vp[i],vp[j]);
            correlationCoefficient.push_back(dabs(currentCorrelationCoefficient));
            printf("%lf\n",currentCorrelationCoefficient);
            totalCorrelationCoefficient += dabs(currentCorrelationCoefficient);
            cnt++;
        }
    }
    puts("printVeltkampData out");
    double averageSpacingVariance = totalSpacingVariance/(int)vp.size();
    double averageCorrelationCoefficient = totalCorrelationCoefficient/cnt;

    sort(correlationCoefficient.begin(),correlationCoefficient.end());
    fprintf(fp,"top 5 Correlation Coefficient\n");
    for ( int i = 0 ; i < 5 ; i++ ) 
        fprintf(fp,"%lf ",correlationCoefficient[(int)correlationCoefficient.size()-i-1]);
    fprintf(fp,"\n");
    fprintf(fp,"bottom 5 Correlation Coefficient\n");
    for ( int i = 0 ; i < 5 ; i++ ) 
        fprintf(fp,"%lf ",correlationCoefficient[i]);
    fprintf(fp,"\n");

    fprintf(fp,"Spacing Variance : %lf\n",averageSpacingVariance);
    fprintf(fp,"Correlation Coefficient : %lf\n",averageCorrelationCoefficient);
    fclose(fp);
}
