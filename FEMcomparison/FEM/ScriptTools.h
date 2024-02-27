//
// Created by victor on 02/09/2020.
//

#ifndef SCRIPTOOLS_H
#define SCRIPTOOLS_H

#include "DataStructure.h"

void Generate(PreConfig &pConfig);
void Summarize(PreConfig &pConfig);
void GenerateThreadSpan(PreConfig &pConfig);
void GenerateRefSpan(PreConfig &pConfig);
std::string GetDataFileName(std::string &path, PreConfig &pConfig);
void CountDataPerSetup(std::vector<std::vector<int>*> &dataPerSetup, PreConfig &pConfig);

#endif //FEMCOMPARISON_OUTPUT_H
