// CVodeSimTestSimple.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>

#include <cuda_runtime.h>

#include "enums/enum_mar_cell_MKII.cuh"

#include "modules/globals.hpp"
#include "modules/commons.hpp"


clock_t START_TIMER;

char buffer[255];
double concs [4] = {0.0, 33.0, 66.0, 99.0};
// variables for I/O
FILE* fp_vm;
FILE* fp_gate;
    
// timer section
clock_t tic()
/* start timer*/
{
    return START_TIMER = clock();
}

void toc(clock_t start)
{
  /*stop timer*/
    std::cout
        << "Elapsed time: "
        << (clock() - start) / (double)CLOCKS_PER_SEC << "s"
        << std::endl;
}


void get_IC50_data_from_file(const char* file_name);

__global__ void initConsts(double* CONSTANTS, double* RATES, double *STATES){
    STATES[0] = -86.2;
CONSTANTS[0] = 8.314;
CONSTANTS[1] = 310;
CONSTANTS[2] = 96.485;
CONSTANTS[3] = 185;
CONSTANTS[4] = 16404;
CONSTANTS[5] = 10;
CONSTANTS[6] = 1000;
CONSTANTS[7] = 1;
CONSTANTS[8] = -52;
CONSTANTS[9] = 0.03;
CONSTANTS[10] = 5.4;
CONSTANTS[11] = 140;
STATES[1] = 138.3;
STATES[2] = 11.6;
CONSTANTS[12] = 2;
STATES[3] = 0.0002;
CONSTANTS[13] = 5.405;
CONSTANTS[14] = 0.096;
STATES[4] = 0;
STATES[5] = 1;
CONSTANTS[15] = 0.062;
STATES[6] = 0;
CONSTANTS[16] = 14.838;
STATES[7] = 0;
STATES[8] = 0.75;
STATES[9] = 0.75;
CONSTANTS[17] = 0.00029;
CONSTANTS[18] = 0.175;
STATES[10] = 0;
STATES[11] = 1;
STATES[12] = 1;
CONSTANTS[19] = 0.000592;
CONSTANTS[20] = 0.294;
STATES[13] = 1;
STATES[14] = 0;
CONSTANTS[21] = 1.362;
CONSTANTS[22] = 1;
CONSTANTS[23] = 40;
CONSTANTS[24] = 1000;
CONSTANTS[25] = 0.1;
CONSTANTS[26] = 2.5;
CONSTANTS[27] = 0.35;
CONSTANTS[28] = 1.38;
CONSTANTS[29] = 87.5;
CONSTANTS[30] = 0.825;
CONSTANTS[31] = 0.0005;
CONSTANTS[32] = 0.0146;
STATES[15] = 0.2;
STATES[16] = 1;
CONSTANTS[33] = 2;
CONSTANTS[34] = 0.016464;
CONSTANTS[35] = 0.25;
CONSTANTS[36] = 0.008232;
CONSTANTS[37] = 0.00025;
CONSTANTS[38] = 8e-5;
CONSTANTS[39] = 0.000425;
CONSTANTS[40] = 0.15;
CONSTANTS[41] = 0.001;
CONSTANTS[42] = 10;
CONSTANTS[43] = 0.3;
CONSTANTS[44] = 1094;
CONSTANTS[45] = 2.00000;
}

__global__ void computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[8] = 1.00000/(1.00000+exp((STATES[0]+20.0000)/7.00000));
ALGEBRAIC[21] =  1125.00*exp(- pow(STATES[0]+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+exp((25.0000 - STATES[0])/10.0000));
RATES[11] = (ALGEBRAIC[8] - STATES[11])/ALGEBRAIC[21];
ALGEBRAIC[10] = 1.00000/(1.00000+exp((STATES[0]+20.0000)/5.00000));
ALGEBRAIC[23] =  85.0000*exp(- pow(STATES[0]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((STATES[0] - 20.0000)/5.00000))+3.00000;
RATES[13] = (ALGEBRAIC[10] - STATES[13])/ALGEBRAIC[23];
ALGEBRAIC[11] = 1.00000/(1.00000+exp((20.0000 - STATES[0])/6.00000));
ALGEBRAIC[24] =  9.50000*exp(- pow(STATES[0]+40.0000, 2.00000)/1800.00)+0.800000;
RATES[14] = (ALGEBRAIC[11] - STATES[14])/ALGEBRAIC[24];
ALGEBRAIC[12] = (STATES[3]<0.000350000 ? 1.00000/(1.00000+pow(STATES[3]/0.000350000, 6.00000)) : 1.00000/(1.00000+pow(STATES[3]/0.000350000, 16.0000)));
ALGEBRAIC[25] = (ALGEBRAIC[12] - STATES[16])/CONSTANTS[33];
RATES[16] = (ALGEBRAIC[12]>STATES[16]&&STATES[0]>- 60.0000 ? 0.00000 : ALGEBRAIC[25]);
ALGEBRAIC[1] = 1.00000/(1.00000+exp((- 26.0000 - STATES[0])/7.00000));
ALGEBRAIC[14] = 450.000/(1.00000+exp((- 45.0000 - STATES[0])/10.0000));
ALGEBRAIC[27] = 6.00000/(1.00000+exp((STATES[0]+30.0000)/11.5000));
ALGEBRAIC[36] =  1.00000*ALGEBRAIC[14]*ALGEBRAIC[27];
RATES[4] = (ALGEBRAIC[1] - STATES[4])/ALGEBRAIC[36];
ALGEBRAIC[2] = 1.00000/(1.00000+exp((STATES[0]+88.0000)/24.0000));
ALGEBRAIC[15] = 3.00000/(1.00000+exp((- 60.0000 - STATES[0])/20.0000));
ALGEBRAIC[28] = 1.12000/(1.00000+exp((STATES[0] - 60.0000)/20.0000));
ALGEBRAIC[37] =  1.00000*ALGEBRAIC[15]*ALGEBRAIC[28];
RATES[5] = (ALGEBRAIC[2] - STATES[5])/ALGEBRAIC[37];
ALGEBRAIC[3] = 1.00000/(1.00000+exp((- 5.00000 - STATES[0])/14.0000));
ALGEBRAIC[16] = 1100.00/ pow((1.00000+exp((- 10.0000 - STATES[0])/6.00000)), 1.0 / 2);
ALGEBRAIC[29] = 1.00000/(1.00000+exp((STATES[0] - 60.0000)/20.0000));
ALGEBRAIC[38] =  1.00000*ALGEBRAIC[16]*ALGEBRAIC[29];
RATES[6] = (ALGEBRAIC[3] - STATES[6])/ALGEBRAIC[38];
ALGEBRAIC[4] = 1.00000/pow(1.00000+exp((- 56.8600 - STATES[0])/9.03000), 2.00000);
ALGEBRAIC[17] = 1.00000/(1.00000+exp((- 60.0000 - STATES[0])/5.00000));
ALGEBRAIC[30] = 0.100000/(1.00000+exp((STATES[0]+35.0000)/5.00000))+0.100000/(1.00000+exp((STATES[0] - 50.0000)/200.000));
ALGEBRAIC[39] =  1.00000*ALGEBRAIC[17]*ALGEBRAIC[30];
RATES[7] = (ALGEBRAIC[4] - STATES[7])/ALGEBRAIC[39];
ALGEBRAIC[5] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[18] = (STATES[0]<- 40.0000 ?  0.0570000*exp(- (STATES[0]+80.0000)/6.80000) : 0.00000);
ALGEBRAIC[31] = (STATES[0]<- 40.0000 ?  2.70000*exp( 0.0790000*STATES[0])+ 310000.*exp( 0.348500*STATES[0]) : 0.770000/( 0.130000*(1.00000+exp((STATES[0]+10.6600)/- 11.1000))));
ALGEBRAIC[40] = 1.00000/(ALGEBRAIC[18]+ALGEBRAIC[31]);
RATES[8] = (ALGEBRAIC[5] - STATES[8])/ALGEBRAIC[40];
ALGEBRAIC[6] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[19] = (STATES[0]<- 40.0000 ? (( ( - 25428.0*exp( 0.244400*STATES[0]) -  6.94800e-06*exp( - 0.0439100*STATES[0]))*(STATES[0]+37.7800))/1.00000)/(1.00000+exp( 0.311000*(STATES[0]+79.2300))) : 0.00000);
ALGEBRAIC[32] = (STATES[0]<- 40.0000 ? ( 0.0242400*exp( - 0.0105200*STATES[0]))/(1.00000+exp( - 0.137800*(STATES[0]+40.1400))) : ( 0.600000*exp( 0.0570000*STATES[0]))/(1.00000+exp( - 0.100000*(STATES[0]+32.0000))));
ALGEBRAIC[41] = 1.00000/(ALGEBRAIC[19]+ALGEBRAIC[32]);
RATES[9] = (ALGEBRAIC[6] - STATES[9])/ALGEBRAIC[41];
ALGEBRAIC[7] = 1.00000/(1.00000+exp((- 5.00000 - STATES[0])/7.50000));
ALGEBRAIC[20] = 1.40000/(1.00000+exp((- 35.0000 - STATES[0])/13.0000))+0.250000;
ALGEBRAIC[33] = 1.40000/(1.00000+exp((STATES[0]+5.00000)/5.00000));
ALGEBRAIC[42] = 1.00000/(1.00000+exp((50.0000 - STATES[0])/20.0000));
ALGEBRAIC[45] =  1.00000*ALGEBRAIC[20]*ALGEBRAIC[33]+ALGEBRAIC[42];
RATES[10] = (ALGEBRAIC[7] - STATES[10])/ALGEBRAIC[45];
ALGEBRAIC[9] = 1.00000/(1.00000+pow(STATES[3]/0.000325000, 8.00000));
ALGEBRAIC[22] = 0.100000/(1.00000+exp((STATES[3] - 0.000500000)/0.000100000));
ALGEBRAIC[34] = 0.200000/(1.00000+exp((STATES[3] - 0.000750000)/0.000800000));
ALGEBRAIC[43] = (ALGEBRAIC[9]+ALGEBRAIC[22]+ALGEBRAIC[34]+0.230000)/1.46000;
ALGEBRAIC[46] = (ALGEBRAIC[43] - STATES[12])/CONSTANTS[45];
RATES[12] = (ALGEBRAIC[43]>STATES[12]&&STATES[0]>- 60.0000 ? 0.00000 : ALGEBRAIC[46]);
ALGEBRAIC[58] = (( (( CONSTANTS[21]*CONSTANTS[10])/(CONSTANTS[10]+CONSTANTS[22]))*STATES[2])/(STATES[2]+CONSTANTS[23]))/(1.00000+ 0.124500*exp(( - 0.100000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))+ 0.0353000*exp(( - STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])));
ALGEBRAIC[13] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[11]/STATES[2]);
ALGEBRAIC[53] =  CONSTANTS[16]*pow(STATES[7], 3.00000)*STATES[8]*STATES[9]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[54] =  CONSTANTS[17]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[59] = ( CONSTANTS[24]*( exp(( CONSTANTS[27]*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*pow(STATES[2], 3.00000)*CONSTANTS[12] -  exp(( (CONSTANTS[27] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*pow(CONSTANTS[11], 3.00000)*STATES[3]*CONSTANTS[26]))/( (pow(CONSTANTS[29], 3.00000)+pow(CONSTANTS[11], 3.00000))*(CONSTANTS[28]+CONSTANTS[12])*(1.00000+ CONSTANTS[25]*exp(( (CONSTANTS[27] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))));
RATES[2] = ( - (ALGEBRAIC[53]+ALGEBRAIC[54]+ 3.00000*ALGEBRAIC[58]+ 3.00000*ALGEBRAIC[59])*CONSTANTS[3])/( CONSTANTS[4]*CONSTANTS[2]);
ALGEBRAIC[26] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[10]/STATES[1]);
ALGEBRAIC[47] = 0.100000/(1.00000+exp( 0.0600000*((STATES[0] - ALGEBRAIC[26]) - 200.000)));
ALGEBRAIC[48] = ( 3.00000*exp( 0.000200000*((STATES[0] - ALGEBRAIC[26])+100.000))+ 1.00000*exp( 0.100000*((STATES[0] - ALGEBRAIC[26]) - 10.0000)))/(1.00000+exp( - 0.500000*(STATES[0] - ALGEBRAIC[26])));
ALGEBRAIC[49] = ALGEBRAIC[47]/(ALGEBRAIC[47]+ALGEBRAIC[48]);
ALGEBRAIC[50] =  CONSTANTS[13]*ALGEBRAIC[49]* pow((CONSTANTS[10]/5.40000), 1.0 / 2)*(STATES[0] - ALGEBRAIC[26]);
ALGEBRAIC[57] =  CONSTANTS[20]*STATES[14]*STATES[13]*(STATES[0] - ALGEBRAIC[26]);
ALGEBRAIC[51] =  CONSTANTS[14]* pow((CONSTANTS[10]/5.40000), 1.0 / 2)*STATES[4]*STATES[5]*(STATES[0] - ALGEBRAIC[26]);
ALGEBRAIC[35] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log((CONSTANTS[10]+ CONSTANTS[9]*CONSTANTS[11])/(STATES[1]+ CONSTANTS[9]*STATES[2]));
ALGEBRAIC[52] =  CONSTANTS[15]*pow(STATES[6], 2.00000)*(STATES[0] - ALGEBRAIC[35]);
ALGEBRAIC[55] = ( (( CONSTANTS[18]*STATES[10]*STATES[11]*STATES[12]*4.00000*STATES[0]*pow(CONSTANTS[2], 2.00000))/( CONSTANTS[0]*CONSTANTS[1]))*( STATES[3]*exp(( 2.00000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) -  0.341000*CONSTANTS[12]))/(exp(( 2.00000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) - 1.00000);
ALGEBRAIC[44] =  (( 0.500000*CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[12]/STATES[3]);
ALGEBRAIC[56] =  CONSTANTS[19]*(STATES[0] - ALGEBRAIC[44]);
ALGEBRAIC[61] = ( CONSTANTS[32]*(STATES[0] - ALGEBRAIC[26]))/(1.00000+exp((25.0000 - STATES[0])/5.98000));
ALGEBRAIC[60] = ( CONSTANTS[30]*STATES[3])/(STATES[3]+CONSTANTS[31]);
ALGEBRAIC[0] = (VOI -  floor(VOI/CONSTANTS[6])*CONSTANTS[6]>=CONSTANTS[5]&&VOI -  floor(VOI/CONSTANTS[6])*CONSTANTS[6]<=CONSTANTS[5]+CONSTANTS[7] ? CONSTANTS[8] : 0.00000);
RATES[0] = - (ALGEBRAIC[50]+ALGEBRAIC[57]+ALGEBRAIC[51]+ALGEBRAIC[52]+ALGEBRAIC[55]+ALGEBRAIC[58]+ALGEBRAIC[53]+ALGEBRAIC[54]+ALGEBRAIC[59]+ALGEBRAIC[56]+ALGEBRAIC[61]+ALGEBRAIC[60]+ALGEBRAIC[0]);
RATES[1] = ( - ((ALGEBRAIC[50]+ALGEBRAIC[57]+ALGEBRAIC[51]+ALGEBRAIC[52]+ALGEBRAIC[61]+ALGEBRAIC[0]) -  2.00000*ALGEBRAIC[58])*CONSTANTS[3])/( CONSTANTS[4]*CONSTANTS[2]);
ALGEBRAIC[62] =  (( CONSTANTS[34]*pow(STATES[15], 2.00000))/(pow(CONSTANTS[35], 2.00000)+pow(STATES[15], 2.00000))+CONSTANTS[36])*STATES[10]*STATES[16];
ALGEBRAIC[63] = CONSTANTS[39]/(1.00000+pow(CONSTANTS[37], 2.00000)/pow(STATES[3], 2.00000));
ALGEBRAIC[64] =  CONSTANTS[38]*(STATES[15] - STATES[3]);
ALGEBRAIC[65] = (( (- ((ALGEBRAIC[55]+ALGEBRAIC[56]+ALGEBRAIC[60]) -  2.00000*ALGEBRAIC[59])/( 2.00000*CONSTANTS[4]*CONSTANTS[2]))*CONSTANTS[3]+ALGEBRAIC[64]) - ALGEBRAIC[63])+ALGEBRAIC[62];
ALGEBRAIC[67] = 1.00000/(1.00000+( CONSTANTS[40]*CONSTANTS[41])/pow(STATES[3]+CONSTANTS[41], 2.00000));
RATES[3] =  ALGEBRAIC[65]*ALGEBRAIC[67];
ALGEBRAIC[66] =  (CONSTANTS[4]/CONSTANTS[44])*(ALGEBRAIC[63] - (ALGEBRAIC[62]+ALGEBRAIC[64]));
ALGEBRAIC[68] = 1.00000/(1.00000+( CONSTANTS[42]*CONSTANTS[43])/pow(STATES[15]+CONSTANTS[43], 2.00000));
RATES[15] =  ALGEBRAIC[66]*ALGEBRAIC[68];
}
__global__ void computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[8] = 1.00000/(1.00000+exp((STATES[0]+20.0000)/7.00000));
ALGEBRAIC[21] =  1125.00*exp(- pow(STATES[0]+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+exp((25.0000 - STATES[0])/10.0000));
ALGEBRAIC[10] = 1.00000/(1.00000+exp((STATES[0]+20.0000)/5.00000));
ALGEBRAIC[23] =  85.0000*exp(- pow(STATES[0]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((STATES[0] - 20.0000)/5.00000))+3.00000;
ALGEBRAIC[11] = 1.00000/(1.00000+exp((20.0000 - STATES[0])/6.00000));
ALGEBRAIC[24] =  9.50000*exp(- pow(STATES[0]+40.0000, 2.00000)/1800.00)+0.800000;
ALGEBRAIC[12] = (STATES[3]<0.000350000 ? 1.00000/(1.00000+pow(STATES[3]/0.000350000, 6.00000)) : 1.00000/(1.00000+pow(STATES[3]/0.000350000, 16.0000)));
ALGEBRAIC[25] = (ALGEBRAIC[12] - STATES[16])/CONSTANTS[33];
ALGEBRAIC[1] = 1.00000/(1.00000+exp((- 26.0000 - STATES[0])/7.00000));
ALGEBRAIC[14] = 450.000/(1.00000+exp((- 45.0000 - STATES[0])/10.0000));
ALGEBRAIC[27] = 6.00000/(1.00000+exp((STATES[0]+30.0000)/11.5000));
ALGEBRAIC[36] =  1.00000*ALGEBRAIC[14]*ALGEBRAIC[27];
ALGEBRAIC[2] = 1.00000/(1.00000+exp((STATES[0]+88.0000)/24.0000));
ALGEBRAIC[15] = 3.00000/(1.00000+exp((- 60.0000 - STATES[0])/20.0000));
ALGEBRAIC[28] = 1.12000/(1.00000+exp((STATES[0] - 60.0000)/20.0000));
ALGEBRAIC[37] =  1.00000*ALGEBRAIC[15]*ALGEBRAIC[28];
ALGEBRAIC[3] = 1.00000/(1.00000+exp((- 5.00000 - STATES[0])/14.0000));
ALGEBRAIC[16] = 1100.00/ pow((1.00000+exp((- 10.0000 - STATES[0])/6.00000)), 1.0 / 2);
ALGEBRAIC[29] = 1.00000/(1.00000+exp((STATES[0] - 60.0000)/20.0000));
ALGEBRAIC[38] =  1.00000*ALGEBRAIC[16]*ALGEBRAIC[29];
ALGEBRAIC[4] = 1.00000/pow(1.00000+exp((- 56.8600 - STATES[0])/9.03000), 2.00000);
ALGEBRAIC[17] = 1.00000/(1.00000+exp((- 60.0000 - STATES[0])/5.00000));
ALGEBRAIC[30] = 0.100000/(1.00000+exp((STATES[0]+35.0000)/5.00000))+0.100000/(1.00000+exp((STATES[0] - 50.0000)/200.000));
ALGEBRAIC[39] =  1.00000*ALGEBRAIC[17]*ALGEBRAIC[30];
ALGEBRAIC[5] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[18] = (STATES[0]<- 40.0000 ?  0.0570000*exp(- (STATES[0]+80.0000)/6.80000) : 0.00000);
ALGEBRAIC[31] = (STATES[0]<- 40.0000 ?  2.70000*exp( 0.0790000*STATES[0])+ 310000.*exp( 0.348500*STATES[0]) : 0.770000/( 0.130000*(1.00000+exp((STATES[0]+10.6600)/- 11.1000))));
ALGEBRAIC[40] = 1.00000/(ALGEBRAIC[18]+ALGEBRAIC[31]);
ALGEBRAIC[6] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[19] = (STATES[0]<- 40.0000 ? (( ( - 25428.0*exp( 0.244400*STATES[0]) -  6.94800e-06*exp( - 0.0439100*STATES[0]))*(STATES[0]+37.7800))/1.00000)/(1.00000+exp( 0.311000*(STATES[0]+79.2300))) : 0.00000);
ALGEBRAIC[32] = (STATES[0]<- 40.0000 ? ( 0.0242400*exp( - 0.0105200*STATES[0]))/(1.00000+exp( - 0.137800*(STATES[0]+40.1400))) : ( 0.600000*exp( 0.0570000*STATES[0]))/(1.00000+exp( - 0.100000*(STATES[0]+32.0000))));
ALGEBRAIC[41] = 1.00000/(ALGEBRAIC[19]+ALGEBRAIC[32]);
ALGEBRAIC[7] = 1.00000/(1.00000+exp((- 5.00000 - STATES[0])/7.50000));
ALGEBRAIC[20] = 1.40000/(1.00000+exp((- 35.0000 - STATES[0])/13.0000))+0.250000;
ALGEBRAIC[33] = 1.40000/(1.00000+exp((STATES[0]+5.00000)/5.00000));
ALGEBRAIC[42] = 1.00000/(1.00000+exp((50.0000 - STATES[0])/20.0000));
ALGEBRAIC[45] =  1.00000*ALGEBRAIC[20]*ALGEBRAIC[33]+ALGEBRAIC[42];
ALGEBRAIC[9] = 1.00000/(1.00000+pow(STATES[3]/0.000325000, 8.00000));
ALGEBRAIC[22] = 0.100000/(1.00000+exp((STATES[3] - 0.000500000)/0.000100000));
ALGEBRAIC[34] = 0.200000/(1.00000+exp((STATES[3] - 0.000750000)/0.000800000));
ALGEBRAIC[43] = (ALGEBRAIC[9]+ALGEBRAIC[22]+ALGEBRAIC[34]+0.230000)/1.46000;
ALGEBRAIC[46] = (ALGEBRAIC[43] - STATES[12])/CONSTANTS[45];
ALGEBRAIC[58] = (( (( CONSTANTS[21]*CONSTANTS[10])/(CONSTANTS[10]+CONSTANTS[22]))*STATES[2])/(STATES[2]+CONSTANTS[23]))/(1.00000+ 0.124500*exp(( - 0.100000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))+ 0.0353000*exp(( - STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])));
ALGEBRAIC[13] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[11]/STATES[2]);
ALGEBRAIC[53] =  CONSTANTS[16]*pow(STATES[7], 3.00000)*STATES[8]*STATES[9]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[54] =  CONSTANTS[17]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[59] = ( CONSTANTS[24]*( exp(( CONSTANTS[27]*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*pow(STATES[2], 3.00000)*CONSTANTS[12] -  exp(( (CONSTANTS[27] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*pow(CONSTANTS[11], 3.00000)*STATES[3]*CONSTANTS[26]))/( (pow(CONSTANTS[29], 3.00000)+pow(CONSTANTS[11], 3.00000))*(CONSTANTS[28]+CONSTANTS[12])*(1.00000+ CONSTANTS[25]*exp(( (CONSTANTS[27] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))));
ALGEBRAIC[26] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[10]/STATES[1]);
ALGEBRAIC[47] = 0.100000/(1.00000+exp( 0.0600000*((STATES[0] - ALGEBRAIC[26]) - 200.000)));
ALGEBRAIC[48] = ( 3.00000*exp( 0.000200000*((STATES[0] - ALGEBRAIC[26])+100.000))+ 1.00000*exp( 0.100000*((STATES[0] - ALGEBRAIC[26]) - 10.0000)))/(1.00000+exp( - 0.500000*(STATES[0] - ALGEBRAIC[26])));
ALGEBRAIC[49] = ALGEBRAIC[47]/(ALGEBRAIC[47]+ALGEBRAIC[48]);
ALGEBRAIC[50] =  CONSTANTS[13]*ALGEBRAIC[49]* pow((CONSTANTS[10]/5.40000), 1.0 / 2)*(STATES[0] - ALGEBRAIC[26]);
ALGEBRAIC[57] =  CONSTANTS[20]*STATES[14]*STATES[13]*(STATES[0] - ALGEBRAIC[26]);
ALGEBRAIC[51] =  CONSTANTS[14]* pow((CONSTANTS[10]/5.40000), 1.0 / 2)*STATES[4]*STATES[5]*(STATES[0] - ALGEBRAIC[26]);
ALGEBRAIC[35] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log((CONSTANTS[10]+ CONSTANTS[9]*CONSTANTS[11])/(STATES[1]+ CONSTANTS[9]*STATES[2]));
ALGEBRAIC[52] =  CONSTANTS[15]*pow(STATES[6], 2.00000)*(STATES[0] - ALGEBRAIC[35]);
ALGEBRAIC[55] = ( (( CONSTANTS[18]*STATES[10]*STATES[11]*STATES[12]*4.00000*STATES[0]*pow(CONSTANTS[2], 2.00000))/( CONSTANTS[0]*CONSTANTS[1]))*( STATES[3]*exp(( 2.00000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) -  0.341000*CONSTANTS[12]))/(exp(( 2.00000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) - 1.00000);
ALGEBRAIC[44] =  (( 0.500000*CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[12]/STATES[3]);
ALGEBRAIC[56] =  CONSTANTS[19]*(STATES[0] - ALGEBRAIC[44]);
ALGEBRAIC[61] = ( CONSTANTS[32]*(STATES[0] - ALGEBRAIC[26]))/(1.00000+exp((25.0000 - STATES[0])/5.98000));
ALGEBRAIC[60] = ( CONSTANTS[30]*STATES[3])/(STATES[3]+CONSTANTS[31]);
ALGEBRAIC[0] = (VOI -  floor(VOI/CONSTANTS[6])*CONSTANTS[6]>=CONSTANTS[5]&&VOI -  floor(VOI/CONSTANTS[6])*CONSTANTS[6]<=CONSTANTS[5]+CONSTANTS[7] ? CONSTANTS[8] : 0.00000);
ALGEBRAIC[62] =  (( CONSTANTS[34]*pow(STATES[15], 2.00000))/(pow(CONSTANTS[35], 2.00000)+pow(STATES[15], 2.00000))+CONSTANTS[36])*STATES[10]*STATES[16];
ALGEBRAIC[63] = CONSTANTS[39]/(1.00000+pow(CONSTANTS[37], 2.00000)/pow(STATES[3], 2.00000));
ALGEBRAIC[64] =  CONSTANTS[38]*(STATES[15] - STATES[3]);
ALGEBRAIC[65] = (( (- ((ALGEBRAIC[55]+ALGEBRAIC[56]+ALGEBRAIC[60]) -  2.00000*ALGEBRAIC[59])/( 2.00000*CONSTANTS[4]*CONSTANTS[2]))*CONSTANTS[3]+ALGEBRAIC[64]) - ALGEBRAIC[63])+ALGEBRAIC[62];
ALGEBRAIC[67] = 1.00000/(1.00000+( CONSTANTS[40]*CONSTANTS[41])/pow(STATES[3]+CONSTANTS[41], 2.00000));
ALGEBRAIC[66] =  (CONSTANTS[4]/CONSTANTS[44])*(ALGEBRAIC[63] - (ALGEBRAIC[62]+ALGEBRAIC[64]));
ALGEBRAIC[68] = 1.00000/(1.00000+( CONSTANTS[42]*CONSTANTS[43])/pow(STATES[15]+CONSTANTS[43], 2.00000));
}

drug_t ic50;
__device__ drug_t *d_ic50;
// double ic50[2000][14];
// double *d_ic50[2000][14];

double *d_concs[4];
__device__ double *d_time_step;

// __global__ void toc(clock_t start = START_TIMER);

__global__ void check_data(){
  printf("check data: \n");
  int idx = 14;
  for(int sample_index=0; sample_index<idx; sample_index++){
        printf("%lf|", d_ic50[2][sample_index]);
        }
     //   printf("\n \n");
}

__global__ void set_time_step(
  /*
  as 'adaptive' solver, we need the time step to change in the middle of 
  the process
  since we need to change almost every function to void, I change the 
  return time_step to 
  cudaMemCopy the time_step, 
  */
    double TIME,
    double time_point,
    double max_time_step,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC) {
    double time_step = 0.005;

    if (TIME <= time_point || (TIME - floor(TIME / CONSTANTS[stim_period]) * CONSTANTS[stim_period]) <= time_point) {
        //printf("TIME <= time_point ms\n");
        //return time_step;
        memcpy(d_time_step, &time_step, sizeof(double));
        __syncthreads(); //equivalent to break
        //printf("dV = %lf, time_step = %lf\n",RATES[V] * time_step, time_step);
    }
    else {
        //printf("TIME > time_point ms\n");
        if (std::abs(RATES[V] * time_step) <= 0.2) {//Slow changes in V
            //printf("dV/dt <= 0.2\n");
            time_step = std::abs(0.8 / RATES[V]);
            //Make sure time_step is between 0.005 and max_time_step
            if (time_step < 0.005) {
                time_step = 0.005;
            }
            else if (time_step > max_time_step) {
                time_step = max_time_step;
            }
            //printf("dV = %lf, time_step = %lf\n",std::abs(RATES[V] * time_step), time_step);
        }
        else if (std::abs(RATES[V] * time_step) >= 0.8) {//Fast changes in V
            //printf("dV/dt >= 0.8\n");
            time_step = std::abs(0.2 / RATES[V]);
            while (std::abs(RATES[V] * time_step) >= 0.8 && 0.005 < time_step && time_step < max_time_step) {
                time_step = time_step / 10.0;
                //printf("dV = %lf, time_step = %lf\n",std::abs(RATES[V] * time_step), time_step);
            }
        }
        // return time_step;
        memcpy(d_time_step, &time_step, sizeof(double));
    }
}
__global__ void solveAnalytical(double* CONSTANTS, double* RATES, double *STATES, double *ALGEBRAIC, double dt)
{
  ////==============
  ////Exact solution
  ////==============
  ////INa
  STATES[m] = ALGEBRAIC[mss] - (ALGEBRAIC[mss] - STATES[m]) * exp(-dt / ALGEBRAIC[tm]);
  STATES[hf] = ALGEBRAIC[hss] - (ALGEBRAIC[hss] - STATES[hf]) * exp(-dt / ALGEBRAIC[thf]);
  STATES[hs] = ALGEBRAIC[hss] - (ALGEBRAIC[hss] - STATES[hs]) * exp(-dt / ALGEBRAIC[ths]);
  STATES[j] = ALGEBRAIC[jss] - (ALGEBRAIC[jss] - STATES[j]) * exp(-dt / ALGEBRAIC[tj]);
  STATES[hsp] = ALGEBRAIC[hssp] - (ALGEBRAIC[hssp] - STATES[hsp]) * exp(-dt / ALGEBRAIC[thsp]);
  STATES[jp] = ALGEBRAIC[jss] - (ALGEBRAIC[jss] - STATES[jp]) * exp(-dt / ALGEBRAIC[tjp]);
  STATES[mL] = ALGEBRAIC[mLss] - (ALGEBRAIC[mLss] - STATES[mL]) * exp(-dt / ALGEBRAIC[tmL]);
  STATES[hL] = ALGEBRAIC[hLss] - (ALGEBRAIC[hLss] - STATES[hL]) * exp(-dt / CONSTANTS[thL]);
  STATES[hLp] = ALGEBRAIC[hLssp] - (ALGEBRAIC[hLssp] - STATES[hLp]) * exp(-dt / CONSTANTS[thLp]);
  ////Ito
  STATES[a] = ALGEBRAIC[ass] - (ALGEBRAIC[ass] - STATES[a]) * exp(-dt / ALGEBRAIC[ta]);
  STATES[iF] = ALGEBRAIC[iss] - (ALGEBRAIC[iss] - STATES[iF]) * exp(-dt / ALGEBRAIC[tiF]);
  STATES[iS] = ALGEBRAIC[iss] - (ALGEBRAIC[iss] - STATES[iS]) * exp(-dt / ALGEBRAIC[tiS]);
  STATES[ap] = ALGEBRAIC[assp] - (ALGEBRAIC[assp] - STATES[ap]) * exp(-dt / ALGEBRAIC[ta]);
  STATES[iFp] = ALGEBRAIC[iss] - (ALGEBRAIC[iss] - STATES[iFp]) * exp(-dt / ALGEBRAIC[tiFp]);
  STATES[iSp] = ALGEBRAIC[iss] - (ALGEBRAIC[iss] - STATES[iSp]) * exp(-dt / ALGEBRAIC[tiSp]);
  ////ICaL
  STATES[d] = ALGEBRAIC[dss] - (ALGEBRAIC[dss] - STATES[d]) * exp(-dt / ALGEBRAIC[td]);
  STATES[ff] = ALGEBRAIC[fss] - (ALGEBRAIC[fss] - STATES[ff]) * exp(-dt / ALGEBRAIC[tff]);
  STATES[fs] = ALGEBRAIC[fss] - (ALGEBRAIC[fss] - STATES[fs]) * exp(-dt / ALGEBRAIC[tfs]);
  STATES[fcaf] = ALGEBRAIC[fcass] - (ALGEBRAIC[fcass] - STATES[fcaf]) * exp(-dt / ALGEBRAIC[tfcaf]);
  STATES[fcas] = ALGEBRAIC[fcass] - (ALGEBRAIC[fcass] - STATES[fcas]) * exp(-dt / ALGEBRAIC[tfcas]);
  STATES[jca] = ALGEBRAIC[fcass] - (ALGEBRAIC[fcass] - STATES[jca]) * exp(- dt / CONSTANTS[tjca]);
  STATES[ffp] = ALGEBRAIC[fss] - (ALGEBRAIC[fss] - STATES[ffp]) * exp(-dt / ALGEBRAIC[tffp]);
  STATES[fcafp] = ALGEBRAIC[fcass] - (ALGEBRAIC[fcass] - STATES[fcafp]) * exp(-d / ALGEBRAIC[tfcafp]);
  STATES[nca] = ALGEBRAIC[anca] * CONSTANTS[k2n] / ALGEBRAIC[km2n] -
      (ALGEBRAIC[anca] * CONSTANTS[k2n] / ALGEBRAIC[km2n] - STATES[nca]) * exp(-ALGEBRAIC[km2n] * dt);
  ////IKr
  STATES[xrf] = ALGEBRAIC[xrss] - (ALGEBRAIC[xrss] - STATES[xrf]) * exp(-dt / ALGEBRAIC[txrf]);
  STATES[xrs] = ALGEBRAIC[xrss] - (ALGEBRAIC[xrss] - STATES[xrs]) * exp(-dt / ALGEBRAIC[txrs]);
  ////IKs
  STATES[xs1] = ALGEBRAIC[xs1ss] - (ALGEBRAIC[xs1ss] - STATES[xs1]) * exp(-dt / ALGEBRAIC[txs1]);
  STATES[xs2] = ALGEBRAIC[xs2ss] - (ALGEBRAIC[xs2ss] - STATES[xs2]) * exp(-dt / ALGEBRAIC[txs2]);
  ////IK1
  STATES[xk1] = ALGEBRAIC[xk1ss] - (ALGEBRAIC[xk1ss] - STATES[xk1]) * exp(-dt / ALGEBRAIC[txk1]);
  ////INaCa
  ////INaK
  ////IKb
  ////INab
  ////ICab
  ///IpCa
  ////Diffusion fluxes
  ////RyR receptors
  STATES[Jrelnp] = ALGEBRAIC[Jrel_inf] - (ALGEBRAIC[Jrel_inf] - STATES[Jrelnp]) * exp(-dt / ALGEBRAIC[tau_rel]);
  STATES[Jrelp] = ALGEBRAIC[Jrel_infp] - (ALGEBRAIC[Jrel_infp] - STATES[Jrelp]) * exp(-dt / ALGEBRAIC[tau_relp]);
  ////SERCA Pump
  ////Calcium translocation
  //
  ////=============================
  ////Approximated solution (Euler)
  ////=============================
  ////ICaL
  //STATES[jca] = STATES[jca] + RATES[jca] * dt;
  ////CaMK
  STATES[CaMKt] = STATES[CaMKt] + RATES[CaMKt] * dt;
  ////Membrane potential
  STATES[V] = STATES[V] + RATES[V] * dt;
  ////Ion Concentrations and Buffers
  STATES[nai] = STATES[nai] + RATES[nai] * dt;
  STATES[nass] = STATES[nass] + RATES[nass] * dt;
  STATES[ki] = STATES[ki] + RATES[ki] * dt;
  STATES[kss] = STATES[kss] + RATES[kss] * dt;
  STATES[cai] = STATES[cai] + RATES[cai] * dt;
  STATES[cass] = STATES[cass] + RATES[cass] * dt;
  STATES[cansr] = STATES[cansr] + RATES[cansr] * dt;
  STATES[cajsr] = STATES[cajsr] + RATES[cajsr] * dt; 
  //========================
  //Full Euler Approximation
  //========================
  //STATES[V] = STATES[V] + RATES[V] * dt;
  //STATES[CaMKt] = STATES[CaMKt] + RATES[CaMKt] * dt;
  //STATES[cass] = STATES[cass] + RATES[cass] * dt;
  //STATES[nai] = STATES[nai] + RATES[nai] * dt;
  //STATES[nass] = STATES[nass] + RATES[nass] * dt;
  //STATES[ki] = STATES[ki] + RATES[ki] * dt;
  //STATES[kss] = STATES[kss] + RATES[kss] * dt;
  //STATES[cansr] = STATES[cansr] + RATES[cansr] * dt;
  //STATES[cajsr] = STATES[cajsr] + RATES[cajsr] * dt;
  //STATES[cai] = STATES[cai] + RATES[cai] * dt;
  //STATES[m] = STATES[m] + RATES[m] * dt;
  //STATES[hf] = STATES[hf] + RATES[hf] * dt;
  //STATES[hs] = STATES[hs] + RATES[hs] * dt;
  //STATES[j] = STATES[j] + RATES[j] * dt;
  //STATES[hsp] = STATES[hsp] + RATES[hsp] * dt;
  //STATES[jp] = STATES[jp] + RATES[jp] * dt;
  //STATES[mL] = STATES[mL] + RATES[mL] * dt;
  //STATES[hL] = STATES[hL] + RATES[hL] * dt;
  //STATES[hLp] = STATES[hLp] + RATES[hLp] * dt;
  //STATES[a] = STATES[a] + RATES[a] * dt;
  //STATES[iF] = STATES[iF] + RATES[iF] * dt;
  //STATES[iS] = STATES[iS] + RATES[iS] * dt;
  //STATES[ap] = STATES[ap] + RATES[ap] * dt;
  //STATES[iFp] = STATES[iFp] + RATES[iFp] * dt;
  //STATES[iSp] = STATES[iSp] + RATES[iSp] * dt;
  //STATES[d] = STATES[d] + RATES[d] * dt;
  //STATES[ff] = STATES[ff] + RATES[ff] * dt;
  //STATES[fs] = STATES[fs] + RATES[fs] * dt;
  //STATES[fcaf] = STATES[fcaf] + RATES[fcaf] * dt;
  //STATES[fcas] = STATES[fcas] + RATES[fcas] * dt;
  //STATES[jca] = STATES[jca] + RATES[jca] * dt;
  //STATES[ffp] = STATES[ffp] + RATES[ffp] * dt;
  //STATES[fcafp] = STATES[fcafp] + RATES[fcafp] * dt;
  //STATES[nca] = STATES[nca] + RATES[nca] * dt;
  //STATES[xrf] = STATES[xrf] + RATES[xrf] * dt;
  //STATES[xrs] = STATES[xrs] + RATES[xrs] * dt;
  //STATES[xs1] = STATES[xs1] + RATES[xs1] * dt;
  //STATES[xs2] = STATES[xs2] + RATES[xs2] * dt;
  //STATES[xk1] = STATES[xk1] + RATES[xk1] * dt;
  //STATES[Jrelnp] = STATES[Jrelnp] + RATES[Jrelnp] * dt;
  //STATES[Jrelp] = STATES[Jrelp] + RATES[Jrelp] * dt;
}



__global__ void do_drug_sim_analytical(double conc,const param_t* p_param, 
const unsigned short sample_id)
{

  /*
  do drug effect simulation, loop will be replaced with kernel loops
  */
  double tcurr = 0.0, dt = 0.005, dt_set, tmax;
  double max_time_step = 1.0, time_point = 25.0;

  double ic50[14] = {2704.000000, 0.695400, 0.000000, 0.000000,
          50490.000000, 0.627700, 2371.000000, 1.984000,
          1947.000000, 1.473000, 12460.000000, 2.885000,
          53.100000, 1.075000};
  
  // files for storing results
  // time-series result
  FILE *vfp_m, *fp_inet, *fp_gate;

  // features
  double inet, qnet;

  // looping counter
  unsigned short idx = 14;
  
  // simulation parameters
  double dtw = 2.0;
  const char *drug_name = "bepridil";
  const double bcl = 2000;
  const double inet_vm_threshold = -88.0;
  const unsigned short pace_max = 10;
  const unsigned short celltype = 0.;
  const unsigned short last_pace_print = 3;
  const unsigned short last_drug_check_pace = 250;
  const unsigned int print_freq = (1./dt) * dtw;
  unsigned short pace_count = 0;
  unsigned short pace_steepest = 0;
  double* RATES;
  double* STATES;
  double* CONSTANTS;
  double* ALGEBRAIC;
  int num_of_algebraic = 69;
  int num_of_constants = 46;
  int num_of_rates = 17;
  int num_of_states = 17;

  RATES = (double*)malloc((num_of_rates)*sizeof(double));
  STATES = (double*)malloc((num_of_states)*sizeof(double));
  CONSTANTS = (double*)malloc((num_of_constants)*sizeof(double));
  ALGEBRAIC = (double*)malloc((num_of_algebraic)*sizeof(double));

  // apply some cell initialization
  initConsts<<<1,1>>>(CONSTANTS, RATES, STATES);
  //p_cell->initConsts( celltype, conc, ic50.data());
  CONSTANTS[stim_period] = bcl;

  // generate file for time-series output
  // snprintf(buffer, sizeof(buffer), "result/%s_%.2lf_vmcheck_smp%d.plt", 
  //           drug_name, conc, sample_id );
  // fp_vm = fopen( buffer, "w" );
  // snprintf(buffer, sizeof(buffer), "result/%s_%.2lf_gates_smp%d.plt",
  //           drug_name, conc, sample_id);
  // fp_gate = fopen(buffer, "w");
  // printf("drug name: %s , concentration: %.2lf , sample id: %d \n", drug_name, conc, sample_id);
  printf("\n");

  // printf(fp_vm, "%s %s\n", "Time", "Vm");
  //printf("Time: %s Vm: %s\n", "Time", "Vm");
  // fprintf(fp_gate, "Time %s\n", GATES_HEADER); //this is to write headers in results

  tmax = pace_max * bcl;

  while (tcurr < tmax) {
    // dt_set = set_time_step<<<1,1>>>(tcurr,
    //     		   time_point,
		//            max_time_step,
  	// 	         CONSTANTS,
		//            RATES,
		// 	         STATES,
		//            ALGEBRAIC);
    set_time_step<<<1,1>>>(tcurr,
        		   time_point,
		           max_time_step,
  		         CONSTANTS,
		           RATES,
			         STATES,
		           ALGEBRAIC);
              // cudaDeviceSynchronize();
    // printf("set time step\n");
    //printf("timestep pointer: %x \n",d_time_step);
    //dt_set = *d_time_step;
    dt_set = 0.0001;

    // // //Compute all rates at tcurr
    computeRates<<<1,1>>>(tcurr,
		          CONSTANTS,
            	RATES,
		          STATES,
            	ALGEBRAIC);
              // cudaDeviceSynchronize();
    // printf("compute rates at tcurr\n");

    //Compute the correct/accepted time step
    if (floor((tcurr + dt_set) / bcl) == floor(tcurr / bcl)) {
      dt = dt_set;
    }
    else {
      dt = (floor(tcurr / bcl) + 1) * bcl - tcurr;
    }

    //Compute the analytical solution
    solveAnalytical<<<1,1>>>(CONSTANTS, RATES, STATES, ALGEBRAIC, dt);
    //printf("solve analytical done\n");
    
    //=============//
    //Print results//
    //=============//
    // fprintf(fp_vm, "%lf %lf\n", tcurr, STATES[V]);
    // fprintf(fp_gate, "%lf ",tcurr);
    printf("tcurr: %lf States[v]: %lf\n", tcurr, STATES[V]);
    // printf("%lf \n \n",tcurr);    
    // for(idx = 0; idx < p_cell->gates_size; idx++){
    //   fprintf(fp_gate, "%lf ", p_cell->STATES[p_cell->GATES_INDICES[idx]]);
    // }
    // fprintf(fp_gate, "\n");
    printf("\n");
    
    //Next time step
    tcurr = tcurr + dt;
  }

  // clean the memories
  //fclose(fp_vm);
  //fclose(fp_gate);
}


//__global__ void Calculate(double d_ic50[11][14], double concs[4], Cellmodel *p_cell);
__global__ void Concentration(drug_t *d_ic50, double *concs[4]){
  
  /*
  uses block and thread in CUDA to replace concentration loop
  */

  // Get the thread ID.
  int sample_id = threadIdx.x;
  int conc_idx = blockIdx.x;
  //printf("doing calculation loop....\n");
  

  //for now, we hard code the concs
  double h_concs[4] = {0.0, 33.0, 66.0, 99.0};

  //memset(h_concs, -1, sizeof(h_concs));
  //printf("%lf", h_concs[1]);
  // cudaMemcpy(d_p_cell, p_cell, sizeof(Cellmodel), cudaMemcpyHostToDevice);
  // cudaMemcpy(h_concs, concs, 4*sizeof(double), cudaMemcpyDeviceToHost);

  // printf("concentration: %d -> value: %lf\n",conc_idx, h_concs[conc_idx]);
  // printf("Sample_ID: %d\n",sample_id );
  
  
  //       printf("\n");
        // for( const auto &conc: concs )
        // { // begin concentration loop
        // printf("Current Concentration: %lf  ",concs[a]);
        // // execute main simulation function
        // //do_drug_sim(conc, ic50[sample_id],
        // //            NULL, sample_id,
        // //            p_cell, ode_solver, cvode_firsttime);
        // // TODO @IritaSee: paralelise this loop that takes each data 
        
        //WARNING: concs still hard coded
       //do_drug_sim_analytical<<<1,1>>>(h_concs[conc_idx], *d_ic50[sample_id], NULL, sample_id);
       do_drug_sim_analytical<<<1,1>>>(h_concs[conc_idx], NULL, sample_id);

        // } // end concentration loop

}


int main()
{

    // input variables for cell simulation
    double bcl, dt;
    unsigned short pace;

    //prepare memory slots for ic_50 
    cudaSetDevice(0);
    cudaMalloc((drug_t**)&d_ic50, sizeof(drug_t));
    //perpare memory slots for concentration and copy it to the just created mem slots
    cudaMalloc((void**)&d_concs, 4*sizeof(double)); 
    cudaMemcpy(d_concs, concs, 4*sizeof(double), cudaMemcpyHostToDevice);
    //prepare memory slots for p_cell and copy it
    // cudaMalloc((void**)d_p_cell, sizeof(Cellmodel));
    // cudaMemcpy(d_p_cell, p_cell, sizeof(Cellmodel), cudaMemcpyHostToDevice);
    unsigned short idx;
    tic();
    snprintf(buffer, sizeof(buffer),
      "./drugs/bepridil/IC50_samples10.csv");
    //drug_t ic50 = get_IC50_data_from_file(buffer);
    //int data_row = sizeof(ic50)/sizeof(ic50[0]);
    int data_row = 10;
    get_IC50_data_from_file(buffer);
    if(sizeof(ic50)/sizeof(ic50[0]) == 0)
        printf("Something problem with the IC50 file!\n");
    else if(sizeof(ic50)/sizeof(ic50[0]) > 2000)
        printf("Too much input! Maximum sample data is 2000!\n");
    printf("start calculation....\n");
    // dim3 block(32,32);
    //dim3 grid ((columns+block.x-1)/block.x,(rows+block.y-1)/block.y);
    Concentration<<<4,data_row>>>(d_ic50, d_concs );  
    // Calculate(d_ic50, d_concs, d_p_cell );
    //concentration loop fails so i loop it altogether
    cudaDeviceSynchronize();
    toc(START_TIMER);
    // loop to do calculation in each data is replaced by this func
    
    // memory cleaning and finalize the program
    

    return 0;
}

void get_IC50_data_from_file(const char* file_name)
{
  /*get IC50 data from a file*/
  /*caution: keep it host function!*/
  FILE *fp_drugs;
  printf("Reading the data....\n");
  
  char *token;
  //std::array<double,14> temp_array; //make the d_ version as well?
  double temp_array[1][14];
  //unsigned short idx;
  unsigned int idx;

  if( (fp_drugs = fopen(file_name, "r")) == NULL){
    printf("Cannot open file %s\n",
      file_name);
    //return ic50;
  }

  int count = 0;

  fgets(buffer, sizeof(buffer), fp_drugs); // skip header
  while( fgets(buffer, sizeof(buffer), fp_drugs) != NULL )
  { // begin line reading
    token = strtok( buffer, "," );
    idx = 0;
    while( token != NULL )
    { // begin data tokenizing
      temp_array[0][idx] = strtod(token, NULL);
      token = strtok(NULL, ",");
      ic50[count][idx] = temp_array[0][idx];
      idx=idx+1;
    } // end data tokenizing
    for(int sample_index=0; sample_index<idx; sample_index++){
        printf("%lf|", ic50[count][sample_index]);
        }
        printf("\n \n");
    //ic50.push_back(temp_array);
    count = count+1;
  } // end line reading

  fclose(fp_drugs);

  //copy the ic50 to GPU memory
  printf("rows found: %d\n",idx);
  
  cudaMemcpy(d_ic50, ic50, idx * sizeof(drug_t), cudaMemcpyHostToDevice);
  cudaDeviceSynchronize();
  check_data<<<1,1>>>();
  //  printf("device memory sample contents: ");
  //  for(int sample_index=0; sample_index<idx; sample_index++){
  //       printf("%lf|", *d_ic50[1][sample_index]);
  //       }
  //       printf("\n \n");

  //return ic50;
}
