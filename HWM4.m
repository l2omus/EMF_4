close all, clear all; clc; 

%--------------------------------------- Empirical Methods in Finance, Homework 3 -------------------------------------------------%
%----------------------------------------------------------------------------------------------------------------------------------%
%----------------------------------------------------------------------------------------------------------------------------------%
%---   Homework Assignment #3: Test of the Capital Asset Pricing Model                         ------------------------------------%
%---   HEC Lausanne, MScF                                                                      ------------------------------------%
%---   Empirical Methods in Finance, Pr. Erice Jondeau                                         ------------------------------------%
%---   AUTHORS:                                                                                ------------------------------------%
%---   Romain Pauli (09412099) | Ludovic Mojonnet (09413840) | Guillaume Nagy (09417304)       ------------------------------------%
%----------------------------------------------------------------------------------------------------------------------------------%

%% Data Mapping

%--------------------------------------- Section 0 -----------------------------------------------------
% This section is dedicated to import empirical data from
% the excel file and allocate them into arrays of variables.
%-------------------------------------------------------------------------------------------------------

% Import data from spreadsheet
data_stocks=xlsread('DATA_HW3.xls','Stocks');
data_market=xlsread('DATA_HW3.xls','Market');

rmrf = data_market(2:end, 2)/100;
smb = data_market(2:end,3)/100;
hml = data_market(2:end,4)/100;
rf = data_market(2:end,5)/100;

stocks_P = data_stocks(2:end,2:end);

N = length(stocks_P(1,:));
T = length(stocks_P(1:end-1,1));

for i=1:N                    %Create a matrix of returns for each firm and each period
    for t = 1:T
        stocks_R(t, i) = (stocks_P(t+1,i) - stocks_P(t,i))/stocks_P(t, i);
    end
    stocks_Z(:,i) = stocks_R(:,i) - rf; %Create the same matrix for excess returns
end

%% Fama-French Factors

%--------------------------------------- Section 1 -----------------------------------------------------
%
%
%-------------------------------------------------------------------------------------------------------

% Point 1a

% Regression time-series 
regressors_T = [ones(size(rmrf)) rmrf smb	hml];

for i=1:N
	reg_T(i) = ols(stocks_Z(:,i), regressors_T);
	alphas(i,1) = reg_T(i).beta(1);
    betas_m(i,1) = reg_T(i).beta(2);
    betas_smb(i,1) = reg_T(i).beta(3);
    betas_hml(i,1) = reg_T(i).beta(4);
    residuals(:,i) = reg_T(i).resid;
end

% Regression cross-sectional

stocks_mean_Z = mean(stocks_Z, 1)';

regressors_N = [ones(size(alphas)) betas_m betas_smb betas_hml];

reg_N = ols(stocks_mean_Z, regressors_N);
for p = 1:4
	psi_hat(p) = reg_T(i).beta(p);
end
clear p;
prt(reg_N);

% Point 1b

% Estimate covariance matrix
sigma_hat = zeros(N,N);
for t = 1:T
	sigma_hat = sigma_hat + residuals(t,:)'*residuals(t,:);
end
sigma_hat = sigma_hat/T;

mean_z = mean(rmrf);
var_z = var(rmrf);

wald_stat(1) = T*(1+mean_z^2/var_z)^(-1)*alphas'*inv(sigma_hat)*alphas;

wald_stat(2) = chi2inv(0.99, N);
wald_stat(3) = chi2inv(0.95, N);
wald_stat(4) = chi2inv(0.90, N)


