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
CONSTANTS[celltype] = 0;
CONSTANTS[R] = 8314;
CONSTANTS[T] = 310;
CONSTANTS[F] = 96485;
CONSTANTS[cm] = 1;
CONSTANTS[rad] = 0.0011;
CONSTANTS[L] = 0.01;
CONSTANTS[vcell] =  1000.00*3.14000*CONSTANTS[rad]*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[amp] = -80;
CONSTANTS[duration] = 0.5;
CONSTANTS[zna] = 1;
CONSTANTS[zca] = 2;
CONSTANTS[zk] = 1;
CONSTANTS[stim_start] = 10.0;
CONSTANTS[stim_end] = 100000000000000000;
CONSTANTS[stim_period] = 1000.0;
CONSTANTS[step_low] = -150.;
CONSTANTS[step_high] = 0;
CONSTANTS[step_start] = 10;
CONSTANTS[step_end] = 5000;
CONSTANTS[GNa] = 75;
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[nao] = 140;
CONSTANTS[mssV1] = 39.57;
CONSTANTS[mssV2] = 9.871;
CONSTANTS[mtD1] = 6.765;
CONSTANTS[mtD2] = 8.552;
CONSTANTS[mtV1] = 11.64;
CONSTANTS[mtV2] = 34.77;
CONSTANTS[mtV3] = 77.42;
CONSTANTS[mtV4] = 5.955;
CONSTANTS[hssV1] = 82.9;
CONSTANTS[hssV2] = 6.086;
CONSTANTS[Ahf] = 0.99;
CONSTANTS[Ahs] = 1.00000 - CONSTANTS[Ahf];
CONSTANTS[GNaL_b] = 0.0075;
CONSTANTS[GNaL] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GNaL_b]*0.600000 : CONSTANTS[GNaL_b]);
CONSTANTS[thL] = 200;
CONSTANTS[thLp] =  3.00000*CONSTANTS[thL];
CONSTANTS[PNab] = 3.75e-10;
CONSTANTS[Gto_b] = 0.02;
CONSTANTS[Gto] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[Gto_b]*4.00000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[Gto_b]*4.00000 : CONSTANTS[Gto_b]);
CONSTANTS[ko] = 5.4;
CONSTANTS[GKr_b] = 0.046;
CONSTANTS[GKr] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GKr_b]*1.30000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[GKr_b]*0.800000 : CONSTANTS[GKr_b]);
CONSTANTS[GKs_b] = 0.0034;
CONSTANTS[GKs] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GKs_b]*1.40000 : CONSTANTS[GKs_b]);
CONSTANTS[PKNa] = 0.01833;
CONSTANTS[GK1_b] = 0.1908;
CONSTANTS[GK1] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GK1_b]*1.20000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[GK1_b]*1.30000 : CONSTANTS[GK1_b]);
CONSTANTS[GKb_b] = 0.003;
CONSTANTS[GKb] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GKb_b]*0.600000 : CONSTANTS[GKb_b]);
CONSTANTS[Kmn] = 0.002;
CONSTANTS[k2n] = 1000;
CONSTANTS[tjca] = 75.0000;
CONSTANTS[Aff] = 0.600000;
CONSTANTS[Afs] = 1.00000 - CONSTANTS[Aff];
CONSTANTS[PCa_b] = 0.0001;
CONSTANTS[PCa] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[PCa_b]*1.20000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[PCa_b]*2.50000 : CONSTANTS[PCa_b]);
CONSTANTS[PCaK] =  0.000357400*CONSTANTS[PCa];
CONSTANTS[PCaNa] =  0.00125000*CONSTANTS[PCa];
CONSTANTS[PCap] =  1.10000*CONSTANTS[PCa];
CONSTANTS[PCaKp] =  0.000357400*CONSTANTS[PCap];
CONSTANTS[PCaNap] =  0.00125000*CONSTANTS[PCap];
CONSTANTS[cao] = 1.8;
CONSTANTS[PCab] = 2.5e-8;
CONSTANTS[GpCa] = 0.0005;
CONSTANTS[KmCap] = 0.0005;
CONSTANTS[kasymm] = 12.5;
CONSTANTS[kcaon] = 1.5e6;
CONSTANTS[kcaoff] = 5e3;
CONSTANTS[kna1] = 15;
CONSTANTS[kna2] = 5;
CONSTANTS[kna3] = 88.12;
CONSTANTS[qna] = 0.5224;
CONSTANTS[qca] = 0.167;
CONSTANTS[wnaca] = 5e3;
CONSTANTS[wna] = 6e4;
CONSTANTS[wca] = 6e4;
CONSTANTS[KmCaAct] = 150e-6;
CONSTANTS[Gncx_b] = 0.0008;
CONSTANTS[Gncx] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[Gncx_b]*1.10000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[Gncx_b]*1.40000 : CONSTANTS[Gncx_b]);
CONSTANTS[h10_i] = CONSTANTS[kasymm]+1.00000+ (CONSTANTS[nao]/CONSTANTS[kna1])*(1.00000+CONSTANTS[nao]/CONSTANTS[kna2]);
CONSTANTS[h11_i] = ( CONSTANTS[nao]*CONSTANTS[nao])/( CONSTANTS[h10_i]*CONSTANTS[kna1]*CONSTANTS[kna2]);
CONSTANTS[h12_i] = 1.00000/CONSTANTS[h10_i];
CONSTANTS[k1_i] =  CONSTANTS[h12_i]*CONSTANTS[cao]*CONSTANTS[kcaon];
CONSTANTS[k2_i] = CONSTANTS[kcaoff];
CONSTANTS[k5_i] = CONSTANTS[kcaoff];
CONSTANTS[h10_ss] = CONSTANTS[kasymm]+1.00000+ (CONSTANTS[nao]/CONSTANTS[kna1])*(1.00000+CONSTANTS[nao]/CONSTANTS[kna2]);
CONSTANTS[h11_ss] = ( CONSTANTS[nao]*CONSTANTS[nao])/( CONSTANTS[h10_ss]*CONSTANTS[kna1]*CONSTANTS[kna2]);
CONSTANTS[h12_ss] = 1.00000/CONSTANTS[h10_ss];
CONSTANTS[k1_ss] =  CONSTANTS[h12_ss]*CONSTANTS[cao]*CONSTANTS[kcaon];
CONSTANTS[k2_ss] = CONSTANTS[kcaoff];
CONSTANTS[k5_ss] = CONSTANTS[kcaoff];
CONSTANTS[k1p] = 949.5;
CONSTANTS[k2p] = 687.2;
CONSTANTS[k3p] = 1899;
CONSTANTS[k4p] = 639;
CONSTANTS[k1m] = 182.4;
CONSTANTS[k2m] = 39.4;
CONSTANTS[k3m] = 79300;
CONSTANTS[k4m] = 40;
CONSTANTS[Knai0] = 9.073;
CONSTANTS[Knao0] = 27.78;
CONSTANTS[delta] = -0.155;
CONSTANTS[Kki] = 0.5;
CONSTANTS[Kko] = 0.3582;
CONSTANTS[MgADP] = 0.05;
CONSTANTS[MgATP] = 9.8;
CONSTANTS[H] = 1e-7;
CONSTANTS[Kmgatp] = 1.698e-7;
CONSTANTS[eP] = 4.2;
CONSTANTS[Khp] = 1.698e-7;
CONSTANTS[Knap] = 224;
CONSTANTS[Kxkur] = 292;
CONSTANTS[a2] = CONSTANTS[k2p];
CONSTANTS[a4] = (( CONSTANTS[k4p]*CONSTANTS[MgATP])/CONSTANTS[Kmgatp])/(1.00000+CONSTANTS[MgATP]/CONSTANTS[Kmgatp]);
CONSTANTS[b1] =  CONSTANTS[k1m]*CONSTANTS[MgADP];
CONSTANTS[Pnak_b] = 30;
CONSTANTS[Pnak] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[Pnak_b]*0.900000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[Pnak_b]*0.700000 : CONSTANTS[Pnak_b]);
CONSTANTS[upScale] = (CONSTANTS[celltype]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[bt] = 4.75;
CONSTANTS[btp] =  1.25000*CONSTANTS[bt];
CONSTANTS[a_relp] =  0.500000*CONSTANTS[btp];
CONSTANTS[a_rel] =  0.500000*CONSTANTS[bt];
CONSTANTS[aCaMK] = 0.05;
CONSTANTS[bCaMK] = 0.00068;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[vjsr] =  0.00480000*CONSTANTS[vcell];
CONSTANTS[vnsr] =  0.0552000*CONSTANTS[vcell];
CONSTANTS[cmdnmax_b] = 0.05;
CONSTANTS[cmdnmax] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[cmdnmax_b]*1.30000 : CONSTANTS[cmdnmax_b]);
CONSTANTS[kmcmdn] = 0.00238;
CONSTANTS[trpnmax] = 0.07;
CONSTANTS[kmtrpn] = 0.0005;
CONSTANTS[BSLmax] = 1.124;
CONSTANTS[BSRmax] = 0.047;
CONSTANTS[KmBSR] = 0.00087;
CONSTANTS[KmBSL] = 0.0087;
CONSTANTS[csqnmax] = 10;
CONSTANTS[kmcsqn] = 0.8;
STATES[m] = 0;
STATES[j] = 1;
STATES[jp] = 1;
STATES[hf] = 1;
STATES[hs] = 1;
STATES[hsp] = 1;
STATES[V] = -87;
STATES[CaMKt] = 0;
STATES[cass] = 1e-4;
STATES[nai] = 7;
STATES[mL] = 0;
STATES[hL] = 1;
STATES[hLp] = 1;
STATES[a] = 0;
STATES[ap] = 0;
STATES[ki] = 145;
STATES[iF] = 1;
STATES[iS] = 1;
STATES[iFp] = 1;
STATES[iSp] = 1;
STATES[xrf] = 0;
STATES[xrs] = 0;
STATES[xs1] = 0;
STATES[xs2] = 0;
STATES[cai] = 1e-4;
STATES[xk1] = 1;
STATES[d] = 0;
STATES[ff] = 1;
STATES[fs] = 1;
STATES[fcaf] = 1;
STATES[nca] = 0;
STATES[jca] = 1;
STATES[fcas] = 1;
STATES[ffp] = 1;
STATES[fcafp] = 1;
STATES[kss] = 145;
STATES[nass] = 7;
STATES[cansr] = 1.2;
STATES[Jrelnp] = 0;
STATES[Jrelp] = 0;
STATES[cajsr] = 1.2;
}

__global__ void computeRates(double TIME, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[vffrt] = ( STATES[V]*CONSTANTS[F]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]);
ALGEBRAIC[vfrt] = ( STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]);
ALGEBRAIC[Istim] = (TIME>=CONSTANTS[stim_start]&&TIME<=CONSTANTS[stim_end]&&(TIME - CONSTANTS[stim_start]) -  floor((TIME - CONSTANTS[stim_start])/CONSTANTS[stim_period])*CONSTANTS[stim_period]<=CONSTANTS[duration] ? CONSTANTS[amp] : 0.00000);
ALGEBRAIC[mss] = 1.00000/(1.00000+exp(- (STATES[V]+CONSTANTS[mssV1])/CONSTANTS[mssV2]));
ALGEBRAIC[tm] = 1.00000/( CONSTANTS[mtD1]*exp((STATES[V]+CONSTANTS[mtV1])/CONSTANTS[mtV2])+ CONSTANTS[mtD2]*exp(- (STATES[V]+CONSTANTS[mtV3])/CONSTANTS[mtV4]));
ALGEBRAIC[hss] = 1.00000/(1.00000+exp((STATES[V]+CONSTANTS[hssV1])/CONSTANTS[hssV2]));
ALGEBRAIC[ths] = 1.00000/( 0.00979400*exp(- (STATES[V]+17.9500)/28.0500)+ 0.334300*exp((STATES[V]+5.73000)/56.6600));
ALGEBRAIC[thf] = 1.00000/( 1.43200e-05*exp(- (STATES[V]+1.19600)/6.28500)+ 6.14900*exp((STATES[V]+0.509600)/20.2700));
ALGEBRAIC[h] =  CONSTANTS[Ahf]*STATES[hf]+ CONSTANTS[Ahs]*STATES[hs];
ALGEBRAIC[jss] = ALGEBRAIC[hss];
ALGEBRAIC[tj] = 2.03800+1.00000/( 0.0213600*exp(- (STATES[V]+100.600)/8.28100)+ 0.305200*exp((STATES[V]+0.994100)/38.4500));
ALGEBRAIC[hssp] = 1.00000/(1.00000+exp((STATES[V]+89.1000)/6.08600));
ALGEBRAIC[thsp] =  3.00000*ALGEBRAIC[ths];
ALGEBRAIC[hp] =  CONSTANTS[Ahf]*STATES[hf]+ CONSTANTS[Ahs]*STATES[hsp];
ALGEBRAIC[tjp] =  1.46000*ALGEBRAIC[tj];
ALGEBRAIC[ENa] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[nao]/STATES[nai]);
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[fINap] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[INa] =  CONSTANTS[GNa]*(STATES[V] - ALGEBRAIC[ENa])*pow(STATES[m], 3.00000)*( (1.00000 - ALGEBRAIC[fINap])*ALGEBRAIC[h]*STATES[j]+ ALGEBRAIC[fINap]*ALGEBRAIC[hp]*STATES[jp]);
ALGEBRAIC[mLss] = 1.00000/(1.00000+exp(- (STATES[V]+42.8500)/5.26400));
ALGEBRAIC[tmL] = ALGEBRAIC[tm];
ALGEBRAIC[hLss] = 1.00000/(1.00000+exp((STATES[V]+87.6100)/7.48800));
ALGEBRAIC[hLssp] = 1.00000/(1.00000+exp((STATES[V]+93.8100)/7.48800));
ALGEBRAIC[fINaLp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[INaL] =  CONSTANTS[GNaL]*(STATES[V] - ALGEBRAIC[ENa])*STATES[mL]*( (1.00000 - ALGEBRAIC[fINaLp])*STATES[hL]+ ALGEBRAIC[fINaLp]*STATES[hLp]);
ALGEBRAIC[INab] = ( CONSTANTS[PNab]*ALGEBRAIC[vffrt]*( STATES[nai]*exp(ALGEBRAIC[vfrt]) - CONSTANTS[nao]))/(exp(ALGEBRAIC[vfrt]) - 1.00000);
ALGEBRAIC[ass] = 1.00000/(1.00000+exp(- (STATES[V] - 14.3400)/14.8200));
ALGEBRAIC[ta] = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- (STATES[V] - 18.4099)/29.3814)))+3.50000/(1.00000+exp((STATES[V]+100.000)/29.3814)));
ALGEBRAIC[iss] = 1.00000/(1.00000+exp((STATES[V]+43.9400)/5.71100));
ALGEBRAIC[delta_epi] = (CONSTANTS[celltype]==1.00000 ? 1.00000 - 0.950000/(1.00000+exp((STATES[V]+70.0000)/5.00000)) : 1.00000);
ALGEBRAIC[tiF_b] = 4.56200+1.00000/( 0.393300*exp(- (STATES[V]+100.000)/100.000)+ 0.0800400*exp((STATES[V]+50.0000)/16.5900));
ALGEBRAIC[tiS_b] = 23.6200+1.00000/( 0.00141600*exp(- (STATES[V]+96.5200)/59.0500)+ 1.78000e-08*exp((STATES[V]+114.100)/8.07900));
ALGEBRAIC[tiF] =  ALGEBRAIC[tiF_b]*ALGEBRAIC[delta_epi];
ALGEBRAIC[tiS] =  ALGEBRAIC[tiS_b]*ALGEBRAIC[delta_epi];
ALGEBRAIC[AiF] = 1.00000/(1.00000+exp((STATES[V] - 213.600)/151.200));
ALGEBRAIC[AiS] = 1.00000 - ALGEBRAIC[AiF];
ALGEBRAIC[i] =  ALGEBRAIC[AiF]*STATES[iF]+ ALGEBRAIC[AiS]*STATES[iS];
ALGEBRAIC[assp] = 1.00000/(1.00000+exp(- (STATES[V] - 24.3400)/14.8200));
ALGEBRAIC[dti_develop] = 1.35400+0.000100000/(exp((STATES[V] - 167.400)/15.8900)+exp(- (STATES[V] - 12.2300)/0.215400));
ALGEBRAIC[dti_recover] = 1.00000 - 0.500000/(1.00000+exp((STATES[V]+70.0000)/20.0000));
ALGEBRAIC[tiFp] = ALGEBRAIC[dti_develop] * ALGEBRAIC[dti_recover] * ALGEBRAIC[tiF];
ALGEBRAIC[tiSp] = ALGEBRAIC[dti_develop] * ALGEBRAIC[dti_recover] * ALGEBRAIC[tiS];
ALGEBRAIC[ip] =  ALGEBRAIC[AiF]*STATES[iFp]+ ALGEBRAIC[AiS]*STATES[iSp];
ALGEBRAIC[EK] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[ko]/STATES[ki]);
ALGEBRAIC[fItop] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Ito] =  CONSTANTS[Gto]*(STATES[V] - ALGEBRAIC[EK])*( (1.00000 - ALGEBRAIC[fItop])*STATES[a]*ALGEBRAIC[i]+ ALGEBRAIC[fItop]*STATES[ap]*ALGEBRAIC[ip]);
ALGEBRAIC[xrss] = 1.00000/(1.00000+exp(- (STATES[V]+8.33700)/6.78900));
ALGEBRAIC[txrf] = 12.9800+1.00000/( 0.365200*exp((STATES[V] - 31.6600)/3.86900)+ 4.12300e-05*exp(- (STATES[V] - 47.7800)/20.3800));
ALGEBRAIC[txrs] = 1.86500+1.00000/( 0.0662900*exp((STATES[V] - 34.7000)/7.35500)+ 1.12800e-05*exp(- (STATES[V] - 29.7400)/25.9400));
ALGEBRAIC[Axrf] = 1.00000/(1.00000+exp((STATES[V]+54.8100)/38.2100));
ALGEBRAIC[Axrs] = 1.00000 - ALGEBRAIC[Axrf];
ALGEBRAIC[xr] =  ALGEBRAIC[Axrf]*STATES[xrf]+ ALGEBRAIC[Axrs]*STATES[xrs];
ALGEBRAIC[rkr] = ( (1.00000/(1.00000+exp((STATES[V]+55.0000)/75.0000)))*1.00000)/(1.00000+exp((STATES[V] - 10.0000)/30.0000));
ALGEBRAIC[IKr] =  CONSTANTS[GKr]* pow((CONSTANTS[ko]/5.40000), 1.0 / 2)*ALGEBRAIC[xr]*ALGEBRAIC[rkr]*(STATES[V] - ALGEBRAIC[EK]);
ALGEBRAIC[xs1ss] = 1.00000/(1.00000+exp(- (STATES[V]+11.6000)/8.93200));
ALGEBRAIC[txs1] = 817.300+1.00000/( 0.000232600*exp((STATES[V]+48.2800)/17.8000)+ 0.00129200*exp(- (STATES[V]+210.000)/230.000));
ALGEBRAIC[xs2ss] = ALGEBRAIC[xs1ss];
ALGEBRAIC[txs2] = 1.00000/( 0.0100000*exp((STATES[V] - 50.0000)/20.0000)+ 0.0193000*exp(- (STATES[V]+66.5400)/31.0000));
ALGEBRAIC[KsCa] = 1.00000+0.600000/(1.00000+pow(3.80000e-05/STATES[cai], 1.40000));
ALGEBRAIC[EKs] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log((CONSTANTS[ko]+ CONSTANTS[PKNa]*CONSTANTS[nao])/(STATES[ki]+ CONSTANTS[PKNa]*STATES[nai]));
ALGEBRAIC[IKs] =  CONSTANTS[GKs]*ALGEBRAIC[KsCa]*STATES[xs1]*STATES[xs2]*(STATES[V] - ALGEBRAIC[EKs]);
ALGEBRAIC[xk1ss] = 1.00000/(1.00000+exp(- (STATES[V]+ 2.55380*CONSTANTS[ko]+144.590)/( 1.56920*CONSTANTS[ko]+3.81150)));
ALGEBRAIC[txk1] = 122.200/(exp(- (STATES[V]+127.200)/20.3600)+exp((STATES[V]+236.800)/69.3300));
ALGEBRAIC[rk1] = 1.00000/(1.00000+exp(((STATES[V]+105.800) -  2.60000*CONSTANTS[ko])/9.49300));
ALGEBRAIC[IK1] =  CONSTANTS[GK1]* pow(CONSTANTS[ko], 1.0 / 2)*ALGEBRAIC[rk1]*STATES[xk1]*(STATES[V] - ALGEBRAIC[EK]);
ALGEBRAIC[xkb] = 1.00000/(1.00000+exp(- (STATES[V] - 14.4800)/18.3400));
ALGEBRAIC[IKb] =  CONSTANTS[GKb]*ALGEBRAIC[xkb]*(STATES[V] - ALGEBRAIC[EK]);
ALGEBRAIC[dss] = 1.00000/(1.00000+exp(- (STATES[V]+3.94000)/4.23000));
ALGEBRAIC[td] = 0.600000+1.00000/(exp( - 0.0500000*(STATES[V]+6.00000))+exp( 0.0900000*(STATES[V]+14.0000)));
ALGEBRAIC[fss] = 1.00000/(1.00000+exp((STATES[V]+19.5800)/3.69600));
ALGEBRAIC[tff] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[V]+20.0000)/10.0000)+ 0.00450000*exp((STATES[V]+20.0000)/10.0000));
ALGEBRAIC[tfs] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[V]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[V]+5.00000)/6.00000));
ALGEBRAIC[f] =  CONSTANTS[Aff]*STATES[ff]+ CONSTANTS[Afs]*STATES[fs];
ALGEBRAIC[fcass] = ALGEBRAIC[fss];
ALGEBRAIC[tfcaf] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[V] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[V] - 4.00000)/7.00000));
ALGEBRAIC[tfcas] = 100.000+1.00000/( 0.000120000*exp(- STATES[V]/3.00000)+ 0.000120000*exp(STATES[V]/7.00000));
ALGEBRAIC[Afcaf] = 0.300000+0.600000/(1.00000+exp((STATES[V] - 10.0000)/10.0000));
ALGEBRAIC[Afcas] = 1.00000 - ALGEBRAIC[Afcaf];
ALGEBRAIC[fca] =  ALGEBRAIC[Afcaf]*STATES[fcaf]+ ALGEBRAIC[Afcas]*STATES[fcas];
ALGEBRAIC[tffp] =  2.50000*ALGEBRAIC[tff];
ALGEBRAIC[fp] =  CONSTANTS[Aff]*STATES[ffp]+ CONSTANTS[Afs]*STATES[fs];
ALGEBRAIC[tfcafp] =  2.50000*ALGEBRAIC[tfcaf];
ALGEBRAIC[fcap] =  ALGEBRAIC[Afcaf]*STATES[fcafp]+ ALGEBRAIC[Afcas]*STATES[fcas];
ALGEBRAIC[km2n] =  STATES[jca]*1.00000;
ALGEBRAIC[anca] = 1.00000/(CONSTANTS[k2n]/ALGEBRAIC[km2n]+pow(1.00000+CONSTANTS[Kmn]/STATES[cass], 4.00000));
ALGEBRAIC[PhiCaL] = ( 4.00000*ALGEBRAIC[vffrt]*( STATES[cass]*exp( 2.00000*ALGEBRAIC[vfrt]) -  0.341000*CONSTANTS[cao]))/(exp( 2.00000*ALGEBRAIC[vfrt]) - 1.00000);
ALGEBRAIC[PhiCaNa] = ( 1.00000*ALGEBRAIC[vffrt]*( 0.750000*STATES[nass]*exp( 1.00000*ALGEBRAIC[vfrt]) -  0.750000*CONSTANTS[nao]))/(exp( 1.00000*ALGEBRAIC[vfrt]) - 1.00000);
ALGEBRAIC[PhiCaK] = ( 1.00000*ALGEBRAIC[vffrt]*( 0.750000*STATES[kss]*exp( 1.00000*ALGEBRAIC[vfrt]) -  0.750000*CONSTANTS[ko]))/(exp( 1.00000*ALGEBRAIC[vfrt]) - 1.00000);
ALGEBRAIC[fICaLp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[ICaL] =  (1.00000 - ALGEBRAIC[fICaLp])*CONSTANTS[PCa]*ALGEBRAIC[PhiCaL]*STATES[d]*( ALGEBRAIC[f]*(1.00000 - STATES[nca])+ STATES[jca]*ALGEBRAIC[fca]*STATES[nca])+ ALGEBRAIC[fICaLp]*CONSTANTS[PCap]*ALGEBRAIC[PhiCaL]*STATES[d]*( ALGEBRAIC[fp]*(1.00000 - STATES[nca])+ STATES[jca]*ALGEBRAIC[fcap]*STATES[nca]);
ALGEBRAIC[ICaNa] =  (1.00000 - ALGEBRAIC[fICaLp])*CONSTANTS[PCaNa]*ALGEBRAIC[PhiCaNa]*STATES[d]*( ALGEBRAIC[f]*(1.00000 - STATES[nca])+ STATES[jca]*ALGEBRAIC[fca]*STATES[nca])+ ALGEBRAIC[fICaLp]*CONSTANTS[PCaNap]*ALGEBRAIC[PhiCaNa]*STATES[d]*( ALGEBRAIC[fp]*(1.00000 - STATES[nca])+ STATES[jca]*ALGEBRAIC[fcap]*STATES[nca]);
ALGEBRAIC[ICaK] =  (1.00000 - ALGEBRAIC[fICaLp])*CONSTANTS[PCaK]*ALGEBRAIC[PhiCaK]*STATES[d]*( ALGEBRAIC[f]*(1.00000 - STATES[nca])+ STATES[jca]*ALGEBRAIC[fca]*STATES[nca])+ ALGEBRAIC[fICaLp]*CONSTANTS[PCaKp]*ALGEBRAIC[PhiCaK]*STATES[d]*( ALGEBRAIC[fp]*(1.00000 - STATES[nca])+ STATES[jca]*ALGEBRAIC[fcap]*STATES[nca]);
ALGEBRAIC[ICab] = ( CONSTANTS[PCab]*4.00000*ALGEBRAIC[vffrt]*( STATES[cai]*exp( 2.00000*ALGEBRAIC[vfrt]) -  0.341000*CONSTANTS[cao]))/(exp( 2.00000*ALGEBRAIC[vfrt]) - 1.00000);
ALGEBRAIC[IpCa] = ( CONSTANTS[GpCa]*STATES[cai])/(CONSTANTS[KmCap]+STATES[cai]);
ALGEBRAIC[hna] = exp(( CONSTANTS[qna]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]));
ALGEBRAIC[hca] = exp(( CONSTANTS[qca]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]));
ALGEBRAIC[h1_i] = 1.00000+ (STATES[nai]/CONSTANTS[kna3])*(1.00000+ALGEBRAIC[hna]);
ALGEBRAIC[h2_i] = ( STATES[nai]*ALGEBRAIC[hna])/( CONSTANTS[kna3]*ALGEBRAIC[h1_i]);
ALGEBRAIC[h3_i] = 1.00000/ALGEBRAIC[h1_i];
ALGEBRAIC[h4_i] = 1.00000+ (STATES[nai]/CONSTANTS[kna1])*(1.00000+STATES[nai]/CONSTANTS[kna2]);
ALGEBRAIC[h5_i] = ( STATES[nai]*STATES[nai])/( ALGEBRAIC[h4_i]*CONSTANTS[kna1]*CONSTANTS[kna2]);
ALGEBRAIC[h6_i] = 1.00000/ALGEBRAIC[h4_i];
ALGEBRAIC[h7_i] = 1.00000+ (CONSTANTS[nao]/CONSTANTS[kna3])*(1.00000+1.00000/ALGEBRAIC[hna]);
ALGEBRAIC[h8_i] = CONSTANTS[nao]/( CONSTANTS[kna3]*ALGEBRAIC[hna]*ALGEBRAIC[h7_i]);
ALGEBRAIC[h9_i] = 1.00000/ALGEBRAIC[h7_i];
ALGEBRAIC[k3p_i] =  ALGEBRAIC[h9_i]*CONSTANTS[wca];
ALGEBRAIC[k3pp_i] =  ALGEBRAIC[h8_i]*CONSTANTS[wnaca];
ALGEBRAIC[k3_i] = ALGEBRAIC[k3p_i]+ALGEBRAIC[k3pp_i];
ALGEBRAIC[k4p_i] = ( ALGEBRAIC[h3_i]*CONSTANTS[wca])/ALGEBRAIC[hca];
ALGEBRAIC[k4pp_i] =  ALGEBRAIC[h2_i]*CONSTANTS[wnaca];
ALGEBRAIC[k4_i] = ALGEBRAIC[k4p_i]+ALGEBRAIC[k4pp_i];
ALGEBRAIC[k6_i] =  ALGEBRAIC[h6_i]*STATES[cai]*CONSTANTS[kcaon];
ALGEBRAIC[k7_i] =  ALGEBRAIC[h5_i]*ALGEBRAIC[h2_i]*CONSTANTS[wna];
ALGEBRAIC[k8_i] =  ALGEBRAIC[h8_i]*CONSTANTS[h11_i]*CONSTANTS[wna];
ALGEBRAIC[x1_i] =  CONSTANTS[k2_i]*ALGEBRAIC[k4_i]*(ALGEBRAIC[k7_i]+ALGEBRAIC[k6_i])+ CONSTANTS[k5_i]*ALGEBRAIC[k7_i]*(CONSTANTS[k2_i]+ALGEBRAIC[k3_i]);
ALGEBRAIC[x2_i] =  CONSTANTS[k1_i]*ALGEBRAIC[k7_i]*(ALGEBRAIC[k4_i]+CONSTANTS[k5_i])+ ALGEBRAIC[k4_i]*ALGEBRAIC[k6_i]*(CONSTANTS[k1_i]+ALGEBRAIC[k8_i]);
ALGEBRAIC[x3_i] =  CONSTANTS[k1_i]*ALGEBRAIC[k3_i]*(ALGEBRAIC[k7_i]+ALGEBRAIC[k6_i])+ ALGEBRAIC[k8_i]*ALGEBRAIC[k6_i]*(CONSTANTS[k2_i]+ALGEBRAIC[k3_i]);
ALGEBRAIC[x4_i] =  CONSTANTS[k2_i]*ALGEBRAIC[k8_i]*(ALGEBRAIC[k4_i]+CONSTANTS[k5_i])+ ALGEBRAIC[k3_i]*CONSTANTS[k5_i]*(CONSTANTS[k1_i]+ALGEBRAIC[k8_i]);
ALGEBRAIC[E1_i] = ALGEBRAIC[x1_i]/(ALGEBRAIC[x1_i]+ALGEBRAIC[x2_i]+ALGEBRAIC[x3_i]+ALGEBRAIC[x4_i]);
ALGEBRAIC[E2_i] = ALGEBRAIC[x2_i]/(ALGEBRAIC[x1_i]+ALGEBRAIC[x2_i]+ALGEBRAIC[x3_i]+ALGEBRAIC[x4_i]);
ALGEBRAIC[E3_i] = ALGEBRAIC[x3_i]/(ALGEBRAIC[x1_i]+ALGEBRAIC[x2_i]+ALGEBRAIC[x3_i]+ALGEBRAIC[x4_i]);
ALGEBRAIC[E4_i] = ALGEBRAIC[x4_i]/(ALGEBRAIC[x1_i]+ALGEBRAIC[x2_i]+ALGEBRAIC[x3_i]+ALGEBRAIC[x4_i]);
ALGEBRAIC[allo_i] = 1.00000/(1.00000+pow(CONSTANTS[KmCaAct]/STATES[cai], 2.00000));
ALGEBRAIC[JncxCa_i] =  ALGEBRAIC[E2_i]*CONSTANTS[k2_i] -  ALGEBRAIC[E1_i]*CONSTANTS[k1_i];
ALGEBRAIC[JncxNa_i] = ( 3.00000*( ALGEBRAIC[E4_i]*ALGEBRAIC[k7_i] -  ALGEBRAIC[E1_i]*ALGEBRAIC[k8_i])+ ALGEBRAIC[E3_i]*ALGEBRAIC[k4pp_i]) -  ALGEBRAIC[E2_i]*ALGEBRAIC[k3pp_i];
ALGEBRAIC[INaCa_i] =  0.800000*CONSTANTS[Gncx]*ALGEBRAIC[allo_i]*( CONSTANTS[zna]*ALGEBRAIC[JncxNa_i]+ CONSTANTS[zca]*ALGEBRAIC[JncxCa_i]);
ALGEBRAIC[h1_ss] = 1.00000+ (STATES[nass]/CONSTANTS[kna3])*(1.00000+ALGEBRAIC[hna]);
ALGEBRAIC[h2_ss] = ( STATES[nass]*ALGEBRAIC[hna])/( CONSTANTS[kna3]*ALGEBRAIC[h1_ss]);
ALGEBRAIC[h3_ss] = 1.00000/ALGEBRAIC[h1_ss];
ALGEBRAIC[h4_ss] = 1.00000+ (STATES[nass]/CONSTANTS[kna1])*(1.00000+STATES[nass]/CONSTANTS[kna2]);
ALGEBRAIC[h5_ss] = ( STATES[nass]*STATES[nass])/( ALGEBRAIC[h4_ss]*CONSTANTS[kna1]*CONSTANTS[kna2]);
ALGEBRAIC[h6_ss] = 1.00000/ALGEBRAIC[h4_ss];
ALGEBRAIC[h7_ss] = 1.00000+ (CONSTANTS[nao]/CONSTANTS[kna3])*(1.00000+1.00000/ALGEBRAIC[hna]);
ALGEBRAIC[h8_ss] = CONSTANTS[nao]/( CONSTANTS[kna3]*ALGEBRAIC[hna]*ALGEBRAIC[h7_ss]);
ALGEBRAIC[h9_ss] = 1.00000/ALGEBRAIC[h7_ss];
ALGEBRAIC[k3p_ss] =  ALGEBRAIC[h9_ss]*CONSTANTS[wca];
ALGEBRAIC[k3pp_ss] =  ALGEBRAIC[h8_ss]*CONSTANTS[wnaca];
ALGEBRAIC[k3_ss] = ALGEBRAIC[k3p_ss]+ALGEBRAIC[k3pp_ss];
ALGEBRAIC[k4p_ss] = ( ALGEBRAIC[h3_ss]*CONSTANTS[wca])/ALGEBRAIC[hca];
ALGEBRAIC[k4pp_ss] =  ALGEBRAIC[h2_ss]*CONSTANTS[wnaca];
ALGEBRAIC[k4_ss] = ALGEBRAIC[k4p_ss]+ALGEBRAIC[k4pp_ss];
ALGEBRAIC[k6_ss] =  ALGEBRAIC[h6_ss]*STATES[cass]*CONSTANTS[kcaon];
ALGEBRAIC[k7_ss] =  ALGEBRAIC[h5_ss]*ALGEBRAIC[h2_ss]*CONSTANTS[wna];
ALGEBRAIC[k8_ss] =  ALGEBRAIC[h8_ss]*CONSTANTS[h11_ss]*CONSTANTS[wna];
ALGEBRAIC[x1_ss] =  CONSTANTS[k2_ss]*ALGEBRAIC[k4_ss]*(ALGEBRAIC[k7_ss]+ALGEBRAIC[k6_ss])+ CONSTANTS[k5_ss]*ALGEBRAIC[k7_ss]*(CONSTANTS[k2_ss]+ALGEBRAIC[k3_ss]);
ALGEBRAIC[x2_ss] =  CONSTANTS[k1_ss]*ALGEBRAIC[k7_ss]*(ALGEBRAIC[k4_ss]+CONSTANTS[k5_ss])+ ALGEBRAIC[k4_ss]*ALGEBRAIC[k6_ss]*(CONSTANTS[k1_ss]+ALGEBRAIC[k8_ss]);
ALGEBRAIC[x3_ss] =  CONSTANTS[k1_ss]*ALGEBRAIC[k3_ss]*(ALGEBRAIC[k7_ss]+ALGEBRAIC[k6_ss])+ ALGEBRAIC[k8_ss]*ALGEBRAIC[k6_ss]*(CONSTANTS[k2_ss]+ALGEBRAIC[k3_ss]);
ALGEBRAIC[x4_ss] =  CONSTANTS[k2_ss]*ALGEBRAIC[k8_ss]*(ALGEBRAIC[k4_ss]+CONSTANTS[k5_ss])+ ALGEBRAIC[k3_ss]*CONSTANTS[k5_ss]*(CONSTANTS[k1_ss]+ALGEBRAIC[k8_ss]);
ALGEBRAIC[E1_ss] = ALGEBRAIC[x1_ss]/(ALGEBRAIC[x1_ss]+ALGEBRAIC[x2_ss]+ALGEBRAIC[x3_ss]+ALGEBRAIC[x4_ss]);
ALGEBRAIC[E2_ss] = ALGEBRAIC[x2_ss]/(ALGEBRAIC[x1_ss]+ALGEBRAIC[x2_ss]+ALGEBRAIC[x3_ss]+ALGEBRAIC[x4_ss]);
ALGEBRAIC[E3_ss] = ALGEBRAIC[x3_ss]/(ALGEBRAIC[x1_ss]+ALGEBRAIC[x2_ss]+ALGEBRAIC[x3_ss]+ALGEBRAIC[x4_ss]);
ALGEBRAIC[E4_ss] = ALGEBRAIC[x4_ss]/(ALGEBRAIC[x1_ss]+ALGEBRAIC[x2_ss]+ALGEBRAIC[x3_ss]+ALGEBRAIC[x4_ss]);
ALGEBRAIC[allo_ss] = 1.00000/(1.00000+pow(CONSTANTS[KmCaAct]/STATES[cass], 2.00000));
ALGEBRAIC[JncxCa_ss] =  ALGEBRAIC[E2_ss]*CONSTANTS[k2_ss] -  ALGEBRAIC[E1_ss]*CONSTANTS[k1_ss];
ALGEBRAIC[JncxNa_ss] = ( 3.00000*( ALGEBRAIC[E4_ss]*ALGEBRAIC[k7_ss] -  ALGEBRAIC[E1_ss]*ALGEBRAIC[k8_ss])+ ALGEBRAIC[E3_ss]*ALGEBRAIC[k4pp_ss]) -  ALGEBRAIC[E2_ss]*ALGEBRAIC[k3pp_ss];
ALGEBRAIC[INaCa_ss] =  0.200000*CONSTANTS[Gncx]*ALGEBRAIC[allo_ss]*( CONSTANTS[zna]*ALGEBRAIC[JncxNa_ss]+ CONSTANTS[zca]*ALGEBRAIC[JncxCa_ss]);
ALGEBRAIC[Knai] =  CONSTANTS[Knai0]*exp(( CONSTANTS[delta]*STATES[V]*CONSTANTS[F])/( 3.00000*CONSTANTS[R]*CONSTANTS[T]));
ALGEBRAIC[Knao] =  CONSTANTS[Knao0]*exp(( (1.00000 - CONSTANTS[delta])*STATES[V]*CONSTANTS[F])/( 3.00000*CONSTANTS[R]*CONSTANTS[T]));
ALGEBRAIC[P] = CONSTANTS[eP]/(1.00000+CONSTANTS[H]/CONSTANTS[Khp]+STATES[nai]/CONSTANTS[Knap]+STATES[ki]/CONSTANTS[Kxkur]);
ALGEBRAIC[a1] = ( CONSTANTS[k1p]*pow(STATES[nai]/ALGEBRAIC[Knai], 3.00000))/((pow(1.00000+STATES[nai]/ALGEBRAIC[Knai], 3.00000)+pow(1.00000+STATES[ki]/CONSTANTS[Kki], 2.00000)) - 1.00000);
ALGEBRAIC[a3] = ( CONSTANTS[k3p]*pow(CONSTANTS[ko]/CONSTANTS[Kko], 2.00000))/((pow(1.00000+CONSTANTS[nao]/ALGEBRAIC[Knao], 3.00000)+pow(1.00000+CONSTANTS[ko]/CONSTANTS[Kko], 2.00000)) - 1.00000);
ALGEBRAIC[b2] = ( CONSTANTS[k2m]*pow(CONSTANTS[nao]/ALGEBRAIC[Knao], 3.00000))/((pow(1.00000+CONSTANTS[nao]/ALGEBRAIC[Knao], 3.00000)+pow(1.00000+CONSTANTS[ko]/CONSTANTS[Kko], 2.00000)) - 1.00000);
ALGEBRAIC[b3] = ( CONSTANTS[k3m]*ALGEBRAIC[P]*CONSTANTS[H])/(1.00000+CONSTANTS[MgATP]/CONSTANTS[Kmgatp]);
ALGEBRAIC[b4] = ( CONSTANTS[k4m]*pow(STATES[ki]/CONSTANTS[Kki], 2.00000))/((pow(1.00000+STATES[nai]/ALGEBRAIC[Knai], 3.00000)+pow(1.00000+STATES[ki]/CONSTANTS[Kki], 2.00000)) - 1.00000);
ALGEBRAIC[x1] =  CONSTANTS[a4]*ALGEBRAIC[a1]*CONSTANTS[a2]+ ALGEBRAIC[b2]*ALGEBRAIC[b4]*ALGEBRAIC[b3]+ CONSTANTS[a2]*ALGEBRAIC[b4]*ALGEBRAIC[b3]+ ALGEBRAIC[b3]*ALGEBRAIC[a1]*CONSTANTS[a2];
ALGEBRAIC[x2] =  ALGEBRAIC[b2]*CONSTANTS[b1]*ALGEBRAIC[b4]+ ALGEBRAIC[a1]*CONSTANTS[a2]*ALGEBRAIC[a3]+ ALGEBRAIC[a3]*CONSTANTS[b1]*ALGEBRAIC[b4]+ CONSTANTS[a2]*ALGEBRAIC[a3]*ALGEBRAIC[b4];
ALGEBRAIC[x3] =  CONSTANTS[a2]*ALGEBRAIC[a3]*CONSTANTS[a4]+ ALGEBRAIC[b3]*ALGEBRAIC[b2]*CONSTANTS[b1]+ ALGEBRAIC[b2]*CONSTANTS[b1]*CONSTANTS[a4]+ ALGEBRAIC[a3]*CONSTANTS[a4]*CONSTANTS[b1];
ALGEBRAIC[x4] =  ALGEBRAIC[b4]*ALGEBRAIC[b3]*ALGEBRAIC[b2]+ ALGEBRAIC[a3]*CONSTANTS[a4]*ALGEBRAIC[a1]+ ALGEBRAIC[b2]*CONSTANTS[a4]*ALGEBRAIC[a1]+ ALGEBRAIC[b3]*ALGEBRAIC[b2]*ALGEBRAIC[a1];
ALGEBRAIC[E1] = ALGEBRAIC[x1]/(ALGEBRAIC[x1]+ALGEBRAIC[x2]+ALGEBRAIC[x3]+ALGEBRAIC[x4]);
ALGEBRAIC[E2] = ALGEBRAIC[x2]/(ALGEBRAIC[x1]+ALGEBRAIC[x2]+ALGEBRAIC[x3]+ALGEBRAIC[x4]);
ALGEBRAIC[E4] = ALGEBRAIC[x4]/(ALGEBRAIC[x1]+ALGEBRAIC[x2]+ALGEBRAIC[x3]+ALGEBRAIC[x4]);
ALGEBRAIC[E3] = ALGEBRAIC[x3]/(ALGEBRAIC[x1]+ALGEBRAIC[x2]+ALGEBRAIC[x3]+ALGEBRAIC[x4]);
ALGEBRAIC[JnakNa] =  3.00000*( ALGEBRAIC[E1]*ALGEBRAIC[a3] -  ALGEBRAIC[E2]*ALGEBRAIC[b3]);
ALGEBRAIC[JnakK] =  2.00000*( ALGEBRAIC[E4]*CONSTANTS[b1] -  ALGEBRAIC[E3]*ALGEBRAIC[a1]);
ALGEBRAIC[INaK] =  CONSTANTS[Pnak]*( CONSTANTS[zna]*ALGEBRAIC[JnakNa]+ CONSTANTS[zk]*ALGEBRAIC[JnakK]);
ALGEBRAIC[Jdiff] = (STATES[cass] - STATES[cai])/0.200000;
ALGEBRAIC[JdiffK] = (STATES[kss] - STATES[ki])/2.00000;
ALGEBRAIC[JdiffNa] = (STATES[nass] - STATES[nai])/2.00000;
ALGEBRAIC[Jupnp] = ( CONSTANTS[upScale]*0.00437500*STATES[cai])/(STATES[cai]+0.000920000);
ALGEBRAIC[Jupp] = ( CONSTANTS[upScale]*2.75000*0.00437500*STATES[cai])/((STATES[cai]+0.000920000) - 0.000170000);
ALGEBRAIC[fJupp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jleak] = ( 0.00393750*STATES[cansr])/15.0000;
ALGEBRAIC[Jup] = ( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak];
ALGEBRAIC[fJrelp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jrel_inf_temp] = ( CONSTANTS[a_rel]*- ALGEBRAIC[ICaL])/(1.00000+ 1.00000*pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[Jrel_inf] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_inf_temp]*1.70000 : ALGEBRAIC[Jrel_inf_temp]);
ALGEBRAIC[tau_relp_temp] = CONSTANTS[bt]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[tau_rel] = (ALGEBRAIC[tau_rel_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_rel_temp]);
ALGEBRAIC[Jrel_temp] = ( CONSTANTS[a_relp]*- ALGEBRAIC[ICaL])/(1.00000+pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[Jrel_infp] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_temp]*1.70000 : ALGEBRAIC[Jrel_temp]);
ALGEBRAIC[tau_rel_temp] = CONSTANTS[bt]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[tau_relp] = (ALGEBRAIC[tau_relp_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_relp_temp]);
ALGEBRAIC[Jrel] =  (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp];
ALGEBRAIC[Jtr] = (STATES[cansr] - STATES[cajsr])/100.000;
ALGEBRAIC[Bcai] = 1.00000/(1.00000+( CONSTANTS[cmdnmax]*CONSTANTS[kmcmdn])/pow(CONSTANTS[kmcmdn]+STATES[cai], 2.00000)+( CONSTANTS[trpnmax]*CONSTANTS[kmtrpn])/pow(CONSTANTS[kmtrpn]+STATES[cai], 2.00000));
ALGEBRAIC[Bcass] = 1.00000/(1.00000+( CONSTANTS[BSRmax]*CONSTANTS[KmBSR])/pow(CONSTANTS[KmBSR]+STATES[cass], 2.00000)+( CONSTANTS[BSLmax]*CONSTANTS[KmBSL])/pow(CONSTANTS[KmBSL]+STATES[cass], 2.00000));
ALGEBRAIC[Bcajsr] = 1.00000/(1.00000+( CONSTANTS[csqnmax]*CONSTANTS[kmcsqn])/pow(CONSTANTS[kmcsqn]+STATES[cajsr], 2.00000));
RATES[m] = (ALGEBRAIC[mss] - STATES[m])/ALGEBRAIC[tm];
RATES[j] = (ALGEBRAIC[jss] - STATES[j])/ALGEBRAIC[tj];
RATES[jp] = (ALGEBRAIC[jss] - STATES[jp])/ALGEBRAIC[tjp];
RATES[hf] = (ALGEBRAIC[hss] - STATES[hf])/ALGEBRAIC[thf];
RATES[hs] = (ALGEBRAIC[hss] - STATES[hs])/ALGEBRAIC[ths];
RATES[hsp] = (ALGEBRAIC[hssp] - STATES[hsp])/ALGEBRAIC[thsp];
RATES[mL] = (ALGEBRAIC[mLss] - STATES[mL])/ALGEBRAIC[tmL];
RATES[hL] = (ALGEBRAIC[hLss] - STATES[hL])/CONSTANTS[thL];
RATES[hLp] = (ALGEBRAIC[hLssp] - STATES[hLp])/CONSTANTS[thLp];
RATES[a] = (ALGEBRAIC[ass] - STATES[a])/ALGEBRAIC[ta];
RATES[ap] = (ALGEBRAIC[assp] - STATES[ap])/ALGEBRAIC[ta];
RATES[iF] = (ALGEBRAIC[iss] - STATES[iF])/ALGEBRAIC[tiF];
RATES[iS] = (ALGEBRAIC[iss] - STATES[iS])/ALGEBRAIC[tiS];
RATES[iFp] = (ALGEBRAIC[iss] - STATES[iFp])/ALGEBRAIC[tiFp];
RATES[iSp] = (ALGEBRAIC[iss] - STATES[iSp])/ALGEBRAIC[tiSp];
RATES[xrf] = (ALGEBRAIC[xrss] - STATES[xrf])/ALGEBRAIC[txrf];
RATES[xrs] = (ALGEBRAIC[xrss] - STATES[xrs])/ALGEBRAIC[txrs];
RATES[xs1] = (ALGEBRAIC[xs1ss] - STATES[xs1])/ALGEBRAIC[txs1];
RATES[xs2] = (ALGEBRAIC[xs2ss] - STATES[xs2])/ALGEBRAIC[txs2];
RATES[xk1] = (ALGEBRAIC[xk1ss] - STATES[xk1])/ALGEBRAIC[txk1];
RATES[d] = (ALGEBRAIC[dss] - STATES[d])/ALGEBRAIC[td];
RATES[ff] = (ALGEBRAIC[fss] - STATES[ff])/ALGEBRAIC[tff];
RATES[fs] = (ALGEBRAIC[fss] - STATES[fs])/ALGEBRAIC[tfs];
RATES[fcaf] = (ALGEBRAIC[fcass] - STATES[fcaf])/ALGEBRAIC[tfcaf];
RATES[nca] =  ALGEBRAIC[anca]*CONSTANTS[k2n] -  STATES[nca]*ALGEBRAIC[km2n];
RATES[jca] = (ALGEBRAIC[fcass] - STATES[jca])/CONSTANTS[tjca];
RATES[fcas] = (ALGEBRAIC[fcass] - STATES[fcas])/ALGEBRAIC[tfcas];
RATES[ffp] = (ALGEBRAIC[fss] - STATES[ffp])/ALGEBRAIC[tffp];
RATES[fcafp] = (ALGEBRAIC[fcass] - STATES[fcafp])/ALGEBRAIC[tfcafp];
RATES[Jrelnp] = (ALGEBRAIC[Jrel_inf] - STATES[Jrelnp])/ALGEBRAIC[tau_rel];
RATES[Jrelp] = (ALGEBRAIC[Jrel_infp] - STATES[Jrelp])/ALGEBRAIC[tau_relp];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
RATES[nai] = ( - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ 3.00000*ALGEBRAIC[INaCa_i]+ 3.00000*ALGEBRAIC[INaK]+ALGEBRAIC[INab])*CONSTANTS[Acap]*CONSTANTS[cm])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffNa]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[nass] = ( - (ALGEBRAIC[ICaNa]+ 3.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffNa];
RATES[ki] = ( - ((ALGEBRAIC[Ito]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[IKb]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[INaK])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffK]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[kss] = ( - ALGEBRAIC[ICaK]*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffK];
RATES[cai] =  ALGEBRAIC[Bcai]*((( - ((ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]) -  2.00000*ALGEBRAIC[INaCa_i])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vmyo]) - ( ALGEBRAIC[Jup]*CONSTANTS[vnsr])/CONSTANTS[vmyo])+( ALGEBRAIC[Jdiff]*CONSTANTS[vss])/CONSTANTS[vmyo]);
RATES[cass] =  ALGEBRAIC[Bcass]*((( - (ALGEBRAIC[ICaL] -  2.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vss])+( ALGEBRAIC[Jrel]*CONSTANTS[vjsr])/CONSTANTS[vss]) - ALGEBRAIC[Jdiff]);
RATES[cansr] = ALGEBRAIC[Jup] - ( ALGEBRAIC[Jtr]*CONSTANTS[vjsr])/CONSTANTS[vnsr];
RATES[cajsr] =  ALGEBRAIC[Bcajsr]*(ALGEBRAIC[Jtr] - ALGEBRAIC[Jrel]);
RATES[V] = - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ALGEBRAIC[Ito]+ALGEBRAIC[ICaL]+ALGEBRAIC[ICaNa]+ALGEBRAIC[ICaK]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[INaCa_i]+ALGEBRAIC[INaCa_ss]+ALGEBRAIC[INaK]+ALGEBRAIC[INab]+ALGEBRAIC[IKb]+ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]+ALGEBRAIC[Istim]);
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
/*__global__ void solveAnalytical(double dt)
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
*/


__global__ void do_drug_sim_analytical(double conc, drug_t *d_ic50, const param_t* p_param, 
const unsigned short sample_id)
{

  /*
  do drug effect simulation, loop will be replaced with kernel loops
  */
  double tcurr = 0.0, dt = 0.005, dt_set, tmax;
  double max_time_step = 1.0, time_point = 25.0;
  
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
  //cudaMalloc((drug_t**)&d_ic50, sizeof(drug_t));
  RATES = (double*)malloc((num_of_rates)*sizeof(double));
  STATES = (double*)malloc((num_of_states)*sizeof(double));
  CONSTANTS = (double*)malloc((num_of_constants)*sizeof(double));
  ALGEBRAIC = (double*)malloc((num_of_algebraic)*sizeof(double));

  // apply some cell initialization
  initConsts<<<1,1>>>(CONSTANTS, RATES, STATES);
  printf("constants: %lf rates: %lf states: %lf \n",CONSTANTS,RATES,STATES);
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
    // set_time_step<<<1,1>>>(tcurr,
    //     		   time_point,
		//            max_time_step,
  	// 	         CONSTANTS,
		//            RATES,
		// 	         STATES,
		//            ALGEBRAIC);
              // cudaDeviceSynchronize();
    // printf("set time step\n");
    //printf("timestep pointer: %x \n",d_time_step);
    //dt_set = *d_time_step;
    dt_set = 0.0001;

    // // //Compute all rates at tcurr
    // computeRates<<<1,1>>>(tcurr,
		//           CONSTANTS,
    //         	RATES,
		//           STATES,
    //         	ALGEBRAIC);
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
    //solveAnalytical<<<1,1>>>(dt);
    //printf("solve analytical done\n");
    
    //=============//
    //Print results//
    //=============//
    // fprintf(fp_vm, "%lf %lf\n", tcurr, STATES[V]);
    // fprintf(fp_gate, "%lf ",tcurr);
    //printf("tcurr: %lf States[V]: %lf\n", tcurr, STATES[V]);
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
       do_drug_sim_analytical<<<1,1>>>(h_concs[conc_idx], d_ic50, NULL, sample_id);

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
