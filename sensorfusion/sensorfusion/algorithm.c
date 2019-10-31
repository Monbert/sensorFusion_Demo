#SensorFusion Algorithm


void generateSDMatrix(float inputData[]){
    float D[][];
}

void calcEigenvaluesAndEigenvectors(float D[][]){
    float T[][];
    float lambda[];
}

void calcPrincipalComponents(float D[][], float T[][]){
    float y[][];
}

void calcContributionRate(float lambda[], float p){
    float alpha[][];
    int m;
    float varphi_m;
}

void cmptISDScore(int m, float alpha[][], float y[][]){
    float Z[][];
}

int eliminateWrongData(float Z[],int n){
    return 0;
    return 1;
}

void cmptFusedOutput(float Z[][], int n, float inputData[][]){
    float outputData[][];
}