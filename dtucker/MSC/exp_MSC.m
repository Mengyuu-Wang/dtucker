clear all; 
close all; 
clc;

addpath('ClusteringMeasure', 'mylib', 'datasets');
addpath("\dtuckerrand")
addpath("\tensor_toolbox-v3.2\")
addpath("\data")
% load
load("MSRC.mat")
rX = [30 7 30 7 5];
lambda = 0.0005;
    
% load('yale.mat');
% rX = [11 15 15 11 3];
% lambda = 0.001;
%% Note: each column is an sample (same as in LRR)
V=length(X);
cls_num = length(unique(gt));
%data preparation...
for v=1:V
    [X{v}]=NormalizeData(X{v});
end
% Initialize...
%% parameter setting
K = length(X); N = size(X{1},2); %sample number
%%%%%%%%%%%%%%%%%%%%%%%% tune 
%%%%%%%%%%%%%%%%%%%%%%%% see MSC; for other algorithms
paras.R{1}=[2 2 2 2 2];
paras.lambda=lambda;
paras.rX   = rX;
        
% -------------------0。5------------------- clustering 
   tic;
   [S] = MSC(X,paras);
   Time=toc;
   for q=1:10
       C=SpectralClustering(S, cls_num);
       [Fi(q),Pi(q),Ri(q)] = compute_f(gt,C);
       [A1 nmi1(q) avgenti(q)] = compute_nmi(gt,C);    
       ACCi(q) = Accuracy(C,double(gt));
       if (min(gt)==0)
           [ARi(q),RIi(q),MIi(q),HIi(q)]=RandIndex(gt+1,C);
       else
           [ARi(q),RIi(q),MIi(q),HIi(q)]=RandIndex(gt,C);
       end       
    end
     
    F= mean(Fi); VF= std(Fi);
    P= mean(Pi); VP= std(Pi);
    Rx= mean(Ri); VR= std(Ri);
    nmi= mean(nmi1); Vnmi= std(nmi1);
    avgent= mean(avgenti); avgent= std(avgenti);
    AR= mean(ARi); VAR= std(ARi);
    ACC=mean(ACCi); VACC=std(ACCi);  
    
    fprintf('F-score: %.3f(%.3f)\n', F,VF);
    fprintf('P: %.3f(%.3f)\n', P,VP);    
    fprintf('R: %.3f(%.3f)\n', Rx,VR);
    fprintf('nmi:%.3f(%.3f)\n', nmi,Vnmi);
    fprintf('AR: %.3f(%.3f)\n', AR,VAR);
    fprintf('ACC: %.3f(%.3f)\n', ACC,VACC);  