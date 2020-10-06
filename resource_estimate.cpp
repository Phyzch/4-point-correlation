//
// Created by phyzch on 4/15/20.
// This part of code estimate the cost of resources in our program
//
#include"system.h"
#include"util.h"
using namespace std;
typedef struct{
    unsigned long size, resident,share,text,lib,data,dt;
}statm_t;
void estimate_memory_cost(ofstream & resource_output)
{ // will use /proc/self/statm file to estimate the resource cost for current process.
    struct timeval ru_utime; // user CPU time used
    struct timeval ru_stime; // system time used.
    int who= RUSAGE_SELF;
    struct rusage usage;
    int ret;
    ret= getrusage(who,&usage);
    if(ret==0){
         ;// output memory usage and CPU time usage here
         ru_utime=usage.ru_utime;
         ru_stime=usage.ru_stime;
         resource_output<<"System time useage: "<<ru_stime.tv_sec<<"."<<ru_stime.tv_usec<<endl;
         // total time is :tv_sec + (1.0/1000000) * tv_usec
         resource_output<<"User CPU time usage:  "<<ru_utime.tv_sec<<"."<<ru_utime.tv_usec<<endl;
         resource_output<<endl;
     }
    else{
         cout<<"Fail to output time and resource usage with getrusage function"<<endl;
     }
    statm_t result;
    // access data in "/proc/self/statm" which store information for memory usage of current process.
    const char * statmpath= "/proc/self/statm";
    FILE * Statm_file = fopen (statmpath,"r");
    if (! Statm_file){
        perror(statmpath);
    }
    fscanf(Statm_file,"%ld %ld  %ld %ld %ld %ld %ld", & result.size,& result.resident,& result.share,& result.text,& result.lib, & result.data,& result.dt);
    // The unit for memory here is pagesize. From command line using getconf PAGESIZE you can see the page size you get. Mine here is 4kB, so I will use this one.
    resource_output<<"The total program size at this time point:  "<<double(result.size)/250<<"MB"<<endl;
    resource_output<<"The resident memory size in RAM:  "<<double(result.resident)/250<<"MB"<<endl;

    resource_output<<"The total data usage  (data + stack):  "<<double(result.data)/250<<"MB"<<endl;
    resource_output<<endl;
    fclose(Statm_file);
}
