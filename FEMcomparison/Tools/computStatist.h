#ifndef computStatist_h
#define computStatist_h
#include<iomanip>
#include<iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <pzvec.h>

using namespace std;

double mean(vector<unsigned long long> vect);
double standDeviation(vector<unsigned long long>&);
double CoefVariation(vector<unsigned long long>&);
//printTable(int dim,int nthreads, int refLevel, bool MKL_FEMcomparison);

inline double mean(vector<unsigned long long> vect){
    vector<unsigned long long>::iterator it;
    double sum=0;
    for(it=vect.begin(); it!=vect.end(); it++)
       sum+=*it;
    return sum/vect.size();
}

inline double standDeviation(vector<unsigned long long>& vect){
    double sum=0;
    vector<unsigned long long>::iterator it;
    double media=mean(vect);
    for(it=vect.begin(); it!=vect.end(); it++)
        sum+=pow(*it-media,2.);
    return pow(sum/vect.size(),0.5);
}

inline double CoefVariation(vector<unsigned long long>& vect){
    double desvStand=standDeviation(vect);
    double media=mean(vect);

    return desvStand/media;
}


inline void printTableAssemble(int dim,bool MKL_FEMcomparison, int refLevel,int nthreads, int test,vector<unsigned long long> vect)
{
    auto time=std::time(nullptr);
    std::stringstream flujo;
    flujo<<std::put_time(std::localtime(&time), "%F_%T");
    flujo<<".csv";
    auto filename=flujo.str();
    std::replace(filename.begin(),filename.end(),':','_');
    std::ofstream ofs(filename, std::ios_base::app);
    ofs<<"H1Hybrid - dim,";
    ofs<<"MKL contribute,";
    ofs<<"Refinement level,";
    ofs<<"Number of threads,";
    ofs<<"Number of tests,";
    ofs<<"Assemble Average time,";
    ofs<<"Coef. Variation %,\n";
    
    ofs<<dim<<",";
    ofs<<MKL_FEMcomparison<<",";
    ofs<<refLevel<<",  ";
    ofs<<nthreads<<",";
    ofs<<test<<",";
    ofs<<mean(vect)<<",";
    ofs<<100*CoefVariation(vect)<<",\n";
    ofs.close();
}
#endif
