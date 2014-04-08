clear all; clc, close all; 

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

rmrf = data_market(:, 2)/100;
smb = data_market(:,3)/100;
hml = data_market(:,4)/100;
rf = data_market(:,5)/100;

stocks_P = data_stocks(2:end,2:end);