#ifndef CELL_HPP
#define CELL_HPP

class Cellmodel
{
protected:
  Cellmodel(){}
public:
  unsigned short algebraic_size;
  unsigned short constants_size;
  unsigned short states_size;
  unsigned short gates_size;
  double *ALGEBRAIC;
  double *CONSTANTS;
  double *RATES;
  double *STATES;
  char GATES_HEADER[255];
  unsigned short *GATES_INDICES;
  virtual ~Cellmodel() {}
  virtual __global__ void initConsts() = 0;
  virtual void initConsts(double type){}
  virtual void initConsts(bool is_dutta){}
  virtual void initConsts(double type, bool is_dutta){}
  virtual void initConsts(double type, double conc, double *hill){}
  virtual void initConsts(double type, double conc, double *hill, bool is_dutta){}
  virtual __global__ void computeRates(double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC) = 0;
  virtual __global__ void solveAnalytical(double dt) {};
};


#endif