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

#include "mar_cell_MKII.cuh"

#include "modules/globals.hpp"
#include "modules/commons.hpp"


clock_t START_TIMER;

char buffer[255];
double concs [4] = {0.0, 33.0, 66.0, 99.0};
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

drug_t ic50;
drug_t *d_ic50;
// double ic50[2000][14];
// double *d_ic50[2000][14];

double *d_concs[4];

// __global__ void toc(clock_t start = START_TIMER);

Cellmodel *p_cell = new mar_cell_MKII();
Cellmodel *d_p_cell = new mar_cell_MKII();

double set_time_step(
  /*
  as 'adaptive' solver, we need the time step to change in the middle of 
  the process
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
        return time_step;
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
        return time_step;
    }
}

__global__ void do_drug_sim_analytical(double conc,double ic50[14],const param_t* p_param, const unsigned short sample_id, Cellmodel *p_cell)
{

  /*
  do drug effect simulation, loop will be replaced with kernel loops
  */
  double tcurr = 0.0, dt = 0.005, dt_set, tmax;
  double max_time_step = 1.0, time_point = 25.0;

  // files for storing results
  // time-series result
  FILE *fp_vm, *fp_inet, *fp_gate;

  // features
  double inet, qnet;

  // looping counter
  unsigned short idx;
  
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

  // apply some cell initialization
  p_cell->initConsts();
  //p_cell->initConsts( celltype, conc, ic50.data());
  p_cell->CONSTANTS[stim_period] = bcl;

  // generate file for time-series output
  snprintf(buffer, sizeof(buffer), "result/%s_%.2lf_vmcheck_smp%d.plt", 
            drug_name, conc, sample_id );
  fp_vm = fopen( buffer, "w" );
  snprintf(buffer, sizeof(buffer), "result/%s_%.2lf_gates_smp%d.plt",
            drug_name, conc, sample_id);
  fp_gate = fopen(buffer, "w");

  fprintf(fp_vm, "%s %s\n", "Time", "Vm");
  fprintf(fp_gate, "Time %s\n", p_cell->GATES_HEADER);

  tmax = pace_max * bcl;

  while (tcurr < tmax) {
    dt_set = set_time_step(tcurr,
        		   time_point,
		           max_time_step,
  		         p_cell->CONSTANTS,
		           p_cell->RATES,
			         p_cell->STATES,
		           p_cell->ALGEBRAIC);

    //Compute all rates at tcurr
    p_cell->computeRates(tcurr,
		          p_cell->CONSTANTS,
            	p_cell->RATES,
		          p_cell->STATES,
            	p_cell->ALGEBRAIC);

    //Compute the correct/accepted time step
    if (floor((tcurr + dt_set) / bcl) == floor(tcurr / bcl)) {
      dt = dt_set;
    }
    else {
      dt = (floor(tcurr / bcl) + 1) * bcl - tcurr;
    }

    //Compute the analytical solution
    p_cell->solveAnalytical(dt);
    
    //=============//
    //Print results//
    //=============//
    fprintf(fp_vm, "%lf %lf\n", tcurr, p_cell->STATES[V]);
    fprintf(fp_gate, "%lf ",tcurr);
    for(idx = 0; idx < p_cell->gates_size; idx++){
      fprintf(fp_gate, "%lf ", p_cell->STATES[p_cell->GATES_INDICES[idx]]);
    }
    fprintf(fp_gate, "\n");
  
    //Next time step
    tcurr = tcurr + dt;
  }

  // clean the memories
  fclose(fp_vm);
  fclose(fp_gate);
}


//__global__ void Calculate(double d_ic50[11][14], double concs[4], Cellmodel *p_cell);
__global__ void Concentration(drug_t *d_ic50, double *concs[4], Cellmodel *p_cell ){
  
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

  printf("concentration: %d -> value: %lf\n",conc_idx, h_concs[conc_idx]);
  printf("Sample_ID: %d\n",sample_id );

  // for( const auto &it1 : d_ic50[sample_id] ){
  //       printf("%lf|", it1);
  //       }
  
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
        do_drug_sim_analytical<<<1,1>>>(h_concs[conc_idx], *d_ic50[sample_id], NULL, sample_id, p_cell);

        // } // end concentration loop

   
   delete p_cell;
}

double set_time_step(double TIME,
    double time_point,
    double max_time_step,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC);

int main()
{

    // input variables for cell simulation
    double bcl, dt;
    unsigned short pace;

    // variables for I/O
    FILE* fp_vm;
    FILE* fp_gate;
    
    //prepare memory slots for ic_50 
    cudaMalloc(&d_ic50, sizeof(drug_t));
    //perpare memory slots for concentration and copy it to the just created mem slots
    cudaMalloc((void**)&d_concs, 4*sizeof(double)); 
    cudaMemcpy(d_concs, concs, 4*sizeof(double), cudaMemcpyHostToDevice);
    //prepare memory slots for p_cell and copy it
    cudaMalloc((void**)d_p_cell, sizeof(Cellmodel));
    cudaMemcpy(d_p_cell, p_cell, sizeof(Cellmodel), cudaMemcpyHostToDevice);
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
    Concentration<<<4,data_row>>>(d_ic50, d_concs, d_p_cell );  
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
    //ic50.push_back(temp_array);
    count = count+1;
  } // end line reading

  fclose(fp_drugs);

  //copy the ic50 to GPU memory
  printf("rows found: %d\n",idx);
  cudaMemcpy(d_ic50, ic50, idx * sizeof(drug_t), cudaMemcpyHostToDevice);

  //return ic50;
}
