addpath("data\")
load("Indian_pines.mat");  %% load("HSV_test.mat")
X = indian_pines;
clear indian_pines

sz = size(X);
maxiters = 50;
tol = 1e-6;
J1 = 400;
J2 = 4*J1;
R_1 = [3 3 5];
R_2 = R_1;

[a, b, cores, ticktock, error, ~] = dtucker_als(X, R_1, R_2,'maxiters',maxiters,'tol',tol);

%[~, ~, ~, ticktock, error, ~] = dtucker_Sampled_als(X, R_1, R_2, J1*ones(size(sz)) ,J2*ones(size(sz)) ,'maxiters',maxiters,'tol',tol);

%[~, ~, ~, ticktock, error, ~] = dtucker_SRFT_als(X, R_1, R_2, J1*ones(size(sz)) ,J2*ones(size(sz)) ,'maxiters',maxiters,'tol',tol);

%[~, ~, ~, ticktock, error, ~] = dtucker_preSRFT_als(X, R_1, R_2, J1*ones(size(sz)) ,J2*ones(size(sz)) ,'maxiters',maxiters,'tol',tol);
