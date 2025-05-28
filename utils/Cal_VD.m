clear all
close all
warning off
clc
experiment=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]; 
for data_ind=[14]
    data=experiment(1,data_ind);
    [HIM,p_index,GT]=read_data(data);   
    [no_lines,no_rows,no_bands]=size(HIM); 
    img0=(ToVector(HIM))';% L*N
% VD
% % HFC
p_hfc=HFC(HIM,0.0001)
% % NWHFC
p_nwhfc=NWHFC(HIM,0.0001)
% % MOCA
groundtruth = 0;
% falseAlarm = [10^(-4) 10^(-3) 10^(-2) 10^(-1)];
falseAlarm = 10^(-4);
whitenChoice = 3;       % 1: residule analysis; 2: spectral correlation; 3: spectral and spatial correlation
baseChoice = 1;         % 1: SVD; 2: ATGP; 3: UNCLS; 4: UFCLS
noiseChoice = 3;        % 1: modified Gaussian;   2: Laplacian;   3: original Gaussian
[p_moca,endindex]=MOCA_set(HIM, groundtruth, falseAlarm, whitenChoice, baseChoice, noiseChoice);
p_moca
end