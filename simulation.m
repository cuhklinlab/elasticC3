clear
clc
close all

%input data (take setting 1 as an example)
load('data/simulated data/simu_rna1.mat');
load('data/simulated data/simu_atac1.mat');
load('data/simulated data/simu_rna_label1.mat');
load('data/simulated data/simu_atac_label1.mat');

i = 1;ind = ((i-1)*100+1):(i*100); 
X = atac_simu(:,ind); clu_X_truth = atac_clu_simu(:,i); % auxiliary data
Y = rna_simu(:,ind); clu_Y_truth = rna_clu_simu(:,i);  % target data

%initialize parameters
nclu_X = 2;nclu_Y = 2;nclu_Z = 3;
niter_outer0 = 2;niter_inter0 = 6;
alpha = 0.1; beta = 0;

%run elasticC3 algorithm   
[~,~,~,~,~, ~,~,~,clu_X_bes,clu_Y_bes] = elasticC3_main(X, Y, nclu_X, nclu_Y, nclu_Z, niter_outer0, niter_inter0, alpha, beta, i, 2); 
[TAB_X, TAB_Y, Eval_tab] = clu_eval(clu_X_truth, clu_Y_truth, clu_X_bes, clu_Y_bes);
       


