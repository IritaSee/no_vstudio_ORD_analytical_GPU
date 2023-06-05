#ifndef MAR_CELL_MKII_HPP
#define MAR_CELL_MKII_HPP

#include "cellmodel.cuh"
#include "enums/enum_mar_cell_MKII.cuh"

class mar_cell_MKII : public Cellmodel{
public:
	mar_cell_MKII();
	~mar_cell_MKII();
	void initConsts();
	void initConsts(double celltype);
	void initConsts(double celltype, double conc, double *ic50 );
	void computeRates(double TIME, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);
	void solveAnalytical( double dt );
private:
	void ___applyDrugEffect(double conc, double *ic50, double epsilon);
	void ___initConsts();
};

#endif
