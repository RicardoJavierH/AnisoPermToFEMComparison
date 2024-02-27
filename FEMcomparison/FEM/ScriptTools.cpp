#include "ScriptTools.h"
#include "InputTreatment.h"

void Generate(PreConfig &pConfig){
    
    pConfig.problem = "EAnisoPerm";
    std::ofstream fileStream;
    pConfig.speedUpOfstream = &fileStream;
    //pConfig.speedUpOfstream->open(pConfig.automatedFileName +"/config.sh",std::ofstream::app);
    
    if( remove("config.sh") != 0) perror( "WARNING: FAIL TO REMOVE CONFIG.SH" );
    else puts( "CONFIG.SH SUCCESSFULLY CLEARED");
    pConfig.speedUpOfstream->open("config.sh",std::ofstream::app);

    fileStream << "#!/bin/bash\n";
    
    OfstreamPath(pConfig);
    
    for (int approxMethod=0; approxMethod < 5; approxMethod++){
        for (int topologyType=1; topologyType < 4; topologyType += 2){
            switch (topologyType) {
                case 0:
                    pConfig.topology = "Triangular";
                    pConfig.refLevel = pConfig.ref2D;
                    pConfig.dim = 2;
                    break;
                case 1:
                    pConfig.topology = "Quadrilateral";
                    pConfig.refLevel = pConfig.ref2D;
                    pConfig.dim = 2;
                    break;
                case 2:
                    pConfig.topology = "Tetrahedral";
                    pConfig.refLevel = pConfig.ref3D;
                    pConfig.dim = 3;
                    break;
                case 3:
                    pConfig.topology = "Hexahedral";
                    pConfig.refLevel = pConfig.ref3D;
                    pConfig.dim = 3;
                    break;
                default:
                    DebugStop();
                    break;
            }
            
            switch (approxMethod) {
                case 0:
                    pConfig.approx = "Mixed";
                    pConfig.n=0;
                    break;
                case 1:
                    pConfig.approx = "Mixed";
                    pConfig.n=1;
                    break;
                case 2:
                    pConfig.approx = "Hybrid";
                    if(topologyType < 2){
                        pConfig.n=2;
                    } else{
                        pConfig.n=3;
                    }
                    break;
                case 3:
                    pConfig.approx = "HybridizedMixed";
                    pConfig.n=1;
                    break;
                case 4:
                pConfig.approx = "Hybrid2";
                if(topologyType < 2){
                        pConfig.n=2;
                    } else{
                        pConfig.n=3;
                    }
                break;
                default:
                    DebugStop();
                    break;
            }
            
            ManageOfstream(pConfig);
            if (pConfig.target.automated)
                GenerateThreadSpan(pConfig);
            else if (pConfig.target.timeEfficiency || pConfig.target.errorMeasurement)
                GenerateRefSpan(pConfig);
            else{
                DebugStop();
            }
        }
            
    }
    pConfig.speeedUpPath->close();
    pConfig.speedUpOfstream->close();
        
}

void Summarize(PreConfig &pConfig){
        
    std::ifstream *PathsIfs = new std::ifstream;
    std::ifstream *DataIfs = new std::ifstream;
        
    pConfig.stat.csv.resize(3);
    pConfig.stat.txt.resize(3);

    for (int i = 0; i < 3; i++){
        pConfig.stat.csv[i] = new std::ofstream;
        pConfig.stat.txt[i] = new std::ofstream;
    }
    if (pConfig.target.errorMeasurement){
        pConfig.stat.txt[3]=new std::ofstream;
    }

    PathsIfs->open("OfstreamPath.txt");
        std::string pathLine;
        int counter = 0;
        while (std::getline(*PathsIfs,pathLine))
        {
            counter++;
            std::string path;
            std::istringstream iss(pathLine);
            iss >> path;
            std::string fileName = GetDataFileName(path,pConfig);

            DataIfs->open(fileName);
            std::string dataLine;
            
            std::vector<std::string> avgName, statName;
            int nStatistics;
            
            if (pConfig.target.errorMeasurement){
                nStatistics = 2;
                statName.push_back("/L2Error");
                statName.push_back("/EnergyError");
            }else{
                nStatistics = 3;
                statName.push_back("/assem");
                statName.push_back("/solve");
                statName.push_back("/total");
            }
            
            for (int i = 0; i < nStatistics; i++){
                avgName.push_back(path + statName[i]);
                std::string csvheader,statType;
                if (pConfig.target.automated){
                    csvheader = "nThreads";
                } else {
                    csvheader = "nRef";
                    if (pConfig.target.errorMeasurement)
                        csvheader = csvheader + ",h,nDof,error";
                }
                
                if (!pConfig.target.errorMeasurement){
                    csvheader = csvheader + ",avg,cvar";
                    statType = "Avg";
                }else{
                    csvheader =  csvheader + ",rate";
                    statType = "Rate";
                }
                statType = statType + ".csv";
                pConfig.stat.csv[i]->open(avgName[i]+statType);
                *pConfig.stat.csv[i] << csvheader;
                std::string suffix = "";
                std::string fileName;

                if (pConfig.target.automated){
                    suffix = ",speedUp";
                    fileName = "SpeedUp.txt";
                }else if (pConfig.target.timeEfficiency){
                    fileName = "timeEfficiency.txt";
                }else if (pConfig.target.errorMeasurement)
                    //fileName = "errorMeasurement.txt";
                    fileName = ".txt";
                *pConfig.stat.csv[i] << suffix << std::endl;

                pConfig.stat.txt[i]->open(avgName[i] + fileName);
                *pConfig.stat.txt[i] << "coordinates{" << std::endl;
            }
            
            if (pConfig.target.errorMeasurement){
                for (int i = nStatistics; i < nStatistics+2; i++){
                    avgName.push_back(path + statName[i-2]+"_h.txt");
                    pConfig.stat.txt[i]->open(avgName[i]);
                    *pConfig.stat.txt[i] << "coordinates{" << std::endl;
                }
            }
            
            int lineCounter = 0;
            int nSpan = -1;
            int loopSize = -1;
            int iterCounter = -1;
            bool iterEnd = true;
            
            int lastSpan = -1;
            while (std::getline(*DataIfs,dataLine)){
        
                if (lineCounter == 0){
                    lineCounter++;
                    continue;
                }
                std::vector<std::string> cellVec;
                
                std::stringstream str(dataLine);
                std::string cell;
                while(getline(str, cell,',')){
                    cellVec.push_back(cell);
                }
                
                ReadSimFile(cellVec,pConfig);
                
                bool newSpanValue = false;
                
                if (pConfig.target.automated){
                    newSpanValue = (nSpan != pConfig.rSimulation.nThreads);
                }else if (pConfig.target.timeEfficiency || pConfig.target.errorMeasurement){
                    newSpanValue = (nSpan != pConfig.rSimulation.nRef);
                }

                if(newSpanValue && iterEnd){
                    if (pConfig.target.automated){
                        nSpan = pConfig.rSimulation.nThreads;
                    }else if (pConfig.target.timeEfficiency || pConfig.target.errorMeasurement){
                        nSpan = pConfig.rSimulation.nRef;
                    }
                    loopSize = pConfig.rSimulation.iterNum;
                    
                    if (loopSize < 1)
                        break;
                    
                    iterCounter = 1;
                    iterEnd = false;
                    pConfig.stat.timeVec.Resize(loopSize,3);
                    
                } else {
                    if (loopSize < 1)
                        break;
                    if (loopSize != 1 && pConfig.rSimulation.iterNum != -iterCounter){
                      std::cout << "loopsize " << loopSize << " pConfig.rSimulation.iterNum " << pConfig.rSimulation.iterNum << "iterCounter " << iterCounter << "iterEnd" << iterEnd << std::endl;
                            break;
                    }
                    if (pConfig.rSimulation.iterNum == 1-loopSize){
                        iterEnd = true;
                    } else {
                        
                        iterEnd = false;
                    }
                    iterCounter++;
                }
                
                int ipos = (iterCounter-1);
                
                if (!pConfig.target.errorMeasurement){
                    pConfig.stat.timeVec(ipos,0) = pConfig.rSimulation.assembleTime;
                    pConfig.stat.timeVec(ipos,1) = pConfig.rSimulation.solveTime;
                    pConfig.stat.timeVec(ipos,2) = pConfig.rSimulation.totalTime;
                        
                    if (iterEnd){
                        for (int i=0; i<3; i++){
                            pConfig.stat.avg[i]=0;
                            pConfig.stat.spu[i]=0;
                            pConfig.stat.cvar[i]=0;
                        }

                        for (int iter = 0; iter < loopSize; iter++){
                            for(int i=0; i < 3; i++)
                                pConfig.stat.avg[i] += pConfig.stat.timeVec(iter,i)/loopSize;
                        }
                        pConfig.stat.serialTime.Resize(3,0);
                        if (pConfig.target.automated){
                            if (pConfig.rSimulation.nThreads == 0){
                                for(int i=0; i < 3; i++){
                                    pConfig.stat.serialTime[i] = pConfig.stat.avg[i];
                                }
                            }
                            for (int i = 0; i < 3; i++){
                                pConfig.stat.spu[i] = pConfig.stat.serialTime[i]/pConfig.stat.avg[i];
                            }
                        }
                        
                        for (int iter = 0; iter < loopSize; iter++){
                            for(int i=0; i < 3; i++){
                                pConfig.stat.cvar[i] += pow(pConfig.stat.timeVec(iter,i)-pConfig.stat.avg[i],2.);
                            }
                        }
                        
                        for(int i=0; i < 3; i++){
                            pConfig.stat.cvar[i] = pow(pConfig.stat.cvar[i]/loopSize,0.5);
                        }
                    
                        if (pConfig.target.automated){
                            for (int i = 0; i < 3; i++){
                                *pConfig.stat.csv[i] << pConfig.rSimulation.nThreads << "," << pConfig.stat.avg[i] << "," << pConfig.stat.cvar[i] << "," << pConfig.stat.spu[i] << std::endl;
                                *pConfig.stat.txt[i] << "(" << pConfig.rSimulation.nThreads << "," << pConfig.stat.spu[i] << ")\n";
                            }
                        } else if (pConfig.target.timeEfficiency){
                            for (int i = 0; i < 3; i++){
                                *pConfig.stat.csv[i] << pConfig.rSimulation.nRef << "," << pConfig.stat.avg[i] << "," << pConfig.stat.cvar[i] << std::endl;
                                *pConfig.stat.txt[i] << "(" << pConfig.stat.avg[i] << ",\n";
                            }
                        }
                    
                        pConfig.stat.timeVec.Resize(0,0);
                        loopSize = -1;
                    
                    }
                }else {
                    std::vector<double> errorVec;
                    errorVec.push_back(pConfig.rSimulation.L2Error);
                    errorVec.push_back(pConfig.rSimulation.energyError);

                        
                    if (lastSpan != -1)
                    {
                        if (pConfig.rSimulation.nRef != lastSpan+1)
                            DebugStop();
                        lastSpan = pConfig.rSimulation.nRef;
                        if (pConfig.rSimulation.nRef > pConfig.rSimulation.maxRef)
                            DebugStop();
                        if(pConfig.rSimulation.nRef == pConfig.rSimulation.maxRef)
                            lastSpan = -1;
                        for (int i=0; i<nStatistics; i++){
                        pConfig.errorRate[i] = log10(errorVec[i]/pConfig.errorLog[i])/log10(0.5);
                        }
                    }else{
                        if (lineCounter == 1){
                            lastSpan = pConfig.rSimulation.nRef;
                            pConfig.errorLog.clear();
                            pConfig.errorLog.resize(2);
                            pConfig.errorRate.clear();
                            pConfig.errorRate.resize(2);
                        }
                    }
                    
                    pConfig.h = pConfig.href0 * pow(2,-pConfig.rSimulation.nRef);

                    for (int i = 0; i < nStatistics; i++){
                        *pConfig.stat.csv[i] << pConfig.rSimulation.nRef << "," << pConfig.h << "," << pConfig.rSimulation.nDof << "," << errorVec[i] << "," << pConfig.errorRate[i] <<std::endl;
                        *pConfig.stat.txt[i] << "," << errorVec[i] << ")\n";
                    }
                    
                    for (int i= nStatistics; i < nStatistics+2; i++){
                        *pConfig.stat.txt[i] << "(" << pConfig.h<< "," << errorVec[i-2] << ")\n";
                    }
                    
                    for (int i=0; i<nStatistics; i++){
                        pConfig.errorLog[i] = errorVec[i];
                        
                    }
                }
                lineCounter++;
            }
            
            DataIfs->close();
            for (int i = 0; i < nStatistics; i++){
                pConfig.stat.csv[i]->close();
                *pConfig.stat.txt[i] << "};" << std::endl;
                pConfig.stat.txt[i]->close();
            }
            if (pConfig.target.errorMeasurement){
                for (int i= nStatistics; i < nStatistics+2; i++){
                    *pConfig.stat.txt[i] << "};" << std::endl;
                    pConfig.stat.txt[i]->close();
                }
            }
        }
}

void GenerateThreadSpan(PreConfig &pConfig){
pConfig.tData.maxThreads = 24;
int iter0 = 10;
int iter;
    for (int nThreads=0; nThreads < pConfig.tData.maxThreads+1; nThreads+=2){
      if (nThreads == 0) iter = iter0;
      else iter = pConfig.stat.nLoops;
        for (int iterNum = 0; iterNum < iter; iterNum++){
            int iterCode = -999;
            if(iterNum == 0){
                iterCode = iter;
            }else {
                iterCode = -iterNum;
            }
            
            ManageOfstream(pConfig, nThreads,iterCode);

        }
    }
}

void GenerateRefSpan(PreConfig &pConfig){
    int ref0 = -1;
    switch (pConfig.dim) {
        case 2:
            pConfig.tData.maxRef = pConfig.ref2D;
            ref0 = 2;
            break;
        case 3:
            pConfig.tData.maxRef = pConfig.ref3D;
            ref0 = 1;
            break;
        default:
            DebugStop();
            break;
    }
    for (int nref=ref0; nref < pConfig.tData.maxRef+1; nref++){
        for (int iterNum = 0; iterNum < pConfig.stat.nLoops; iterNum++){
            int iterCode = -999;
            if(iterNum == 0){
                iterCode = pConfig.stat.nLoops;
            }else {
                iterCode = -iterNum;
            }
            
            ManageOfstream(pConfig,nref,iterCode);

        }
    }
}
std::string GetDataFileName(std::string &path, PreConfig &pConfig){
    std::string fileName;
    if (pConfig.target.automated){
        fileName = path + "speedUp.csv";
    } else if (pConfig.target.timeEfficiency){
            fileName = path + "timeEfficiency.csv";
    } else if (pConfig.target.errorMeasurement){
    
            fileName = path + "errorMeasurement.csv";
        
    }else{
        
        DebugStop();
    }
    return fileName;
}

void CountDataPerSetup(std::vector<std::vector<int>*> &dataPerSetup, PreConfig &pConfig){
    std::ifstream *PathsIfs = new std::ifstream;
    std::ifstream *DataIfs = new std::ifstream;
    
    std::vector<std::string> ofsPathFile;
    ofsPathFile.push_back("OfstreamPathError.txt");
    ofsPathFile.push_back("OfstreamPathTime.txt");

    for (int i = 0; i < 2; i++){
    
        PathsIfs->open(ofsPathFile[i]);
        std::string pathLine;
        int pathCounter = 0;
        
        while (std::getline(*PathsIfs,pathLine))
        {
            pathCounter++;
            (*(dataPerSetup[i])).resize(pathCounter);
            std::string dataLine;
            std::string path;
            std::istringstream iss(pathLine);
            iss >> path;
            
            std::string fileName;
            if (i == 0){
                fileName = path + "L2Error.txt";
            }else{
                fileName = path + "assemtimeEfficiency.txt";

            }
            
            DataIfs->open(fileName);
            int dataCounter = 0;
            while (std::getline(*DataIfs,dataLine)){
                dataCounter++;
            }
            DataIfs->close();
            (*(dataPerSetup[i]))[pathCounter-1] = dataCounter;
        }
        PathsIfs->close();
    }
    if (dataPerSetup[0]->size() != dataPerSetup[1]->size())
        DebugStop();
    
    for (int i = 0; i < dataPerSetup[0]->size(); i++){
        if ((*(dataPerSetup[0]))[i] != (*(dataPerSetup[1]))[i]){
            DebugStop();
        }
    }
}
