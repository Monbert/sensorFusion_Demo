#SensorFusion Algorithm


void generate_sdmatrix(float input_data[]){
    float d[][];
}

void calc_eigenvalues_and_eigenvectors(float d[][]){
    float t[][];
    float lambda[];
}

void calc_principal_components(float d[][], float t[][]){
    float y[][];
}

void calc_contribution_rate(float lambda[], float p){
    float alpha[][];
    int m;
    float varphi_m;
}

void cmpt_isd_score(int m, float alpha[][], float y[][]){
    float z[][];
}

int eliminate_wrong_data(float z[],int n){
    return 0;
    return 1;
}

void cmpt_fused_output(float z[][], int n, float input_data[][]){
    float output_data[][];
}