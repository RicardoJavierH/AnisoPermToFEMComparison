//
// Created by victor on 02/09/2020.
//

#ifndef FEMCOMPARISON_OUTPUT_H
#define FEMCOMPARISON_OUTPUT_H

#include "DataStructure.h"

//// Flush csv file with L2 and semi-H1 errors and rates
void FlushTable(PreConfig &eData,char *argv[]);

//// Erase all errors but L2 and semi-H1
void CleanErrors(std::string file);

//// The errors doesn't match its legend on generated error file
//// This method places the right error with the right legend
void InvertError(std::string file);

//// Fill csv file with L2 and semi-H1 errors and rates
void FillErrors(std::ofstream &table,std::string f,int mode);

//// Fill legend of csv file
void FillLegend(std::ofstream &table,int hash_count, int it_count);

//// Print start and current clock difference
void FlushTime(PreConfig &eData, clock_t start);

#endif //FEMCOMPARISON_OUTPUT_H
