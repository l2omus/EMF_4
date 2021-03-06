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
data_stocks=xlsread('DATA_HW3.xls','Stocks'); % 
data_market=xlsread('DATA_HW3.xls','Market'); % 

rmrf = data_market(2:end, 2)/100;	% Create a vector with data and adjust for percentage
smb = data_market(2:end,3)/100;		% Create a vector with data and adjust for percentage
hml = data_market(2:end,4)/100;		% Create a vector with data and adjust for percentage
rf = data_market(2:end,5)/100;		% Create a vector with data and adjust for percentage

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


%% PCA

%--------------------------------------- Section 2 -----------------------------------------------------
%
%
%-------------------------------------------------------------------------------------------------------

for t=1:T
	stocks_R_cen(t,:) = stocks_R(t,:) - mean(stocks_R); 
	stocks_R_cs(t,:) = stocks_R_cen(t,:)./std(stocks_R);
end

stocks_R_cs_corr = corr(stocks_R_cs);
[eig_vec eig_val] = eig(stocks_R_cs_corr);

pc = (eig_vec'*stocks_R_cen')';


eig_val_v = diag(eig_val);
k = 0;
m = 0;
while m<0.45
	k = k+1;
	m = sum(eig_val_v(1:k))/sum(eig_val_v);
end
message = ['With ', num2str(k) , ' factors we can explain ', num2str(m*100), ' percents of the total variability in the excess returns matrix.'];
disp(message)

pc_crop = pc(:,1:k);

regressors_pc = [ones(length(pc_crop(:,1)),1) pc_crop];
for i=1:N
	reg_pc(i) = ols(stocks_R(:,i), regressors_pc);
	mu(i) = reg_pc(i).beta(1);
	gamma_1(i) = reg_pc(i).beta(2);
	gamma_2(i) = reg_pc(i).beta(3);
	gamma_3(i) = reg_pc(i).beta(4);
end

label = {'3M', 'ATT', 'AE', 'BO', 'CP', 'CHE','COCA', 'PDM', 'EXX','GE', 'HD', 'INT','IBUS', 'GPM', 'J&G','MCD', 'MNC', 'M$','NIKE', 'PZR', 'P&G','TRA', 'UT', 'UGP','VC', 'WM', 'WD',};
figure(1)
bar(gamma_1)
set(gca, 'XTick', 1:N, 'XTickLabel', label)

figure(2)
bar(gamma_2)
set(gca, 'XTick', 1:N, 'XTickLabel', label)

for i=1:N
	weights_1(i) = gamma_1(i)/sum(gamma_1);
end

portfolio_R_1 = (weights_1*stocks_R')';

for i=1:N
	eq_weights(i) = [1/N];
end

portfolio_R_eq = (eq_weights*stocks_R')';

corr_port_R_1_eq = corr(portfolio_R_eq, portfolio_R_1)

for i=1:N
	weights_2(i) = gamma_2(i)/sum(gamma_2);
end

portfolio_R_2 = (weights_2*stocks_R')';

corr_port_R_1_2 = corr(portfolio_R_2, portfolio_R_1)





