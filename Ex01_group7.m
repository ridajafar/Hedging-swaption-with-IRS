%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 1 Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUT data are stored as:
%   1. Market_data: table of ICAP quotes for the strip of IRS
%       Column #1: IRS maturity (year frac)
%       Column #2: MID rate
%   2. A set of scalar variable with self-explanatory names
% 
%    OUTPUT data are stored as:
%   3. ZC_curve: Table of ZC rates (cont. comp. 30/360)
%      Maturities are year fractions
%   4. A set of scalar variable with self-explanatory names
%   
%   5. Required functions' template:
%   5.1. ZC_bootstrap_IRS_only: Analytical bootstrap of ZC curve from IRS
%   function [ZC_curve]=ZC_bootstrap_IRS_only(IRS_data, freq)
%   5.2. s_fwd: IRS forward rate
%   function [ par_yield_fwd ] = s_fwd( t1, t2, freq, ZC_curve )
%   5.3. Swaption_Black_receiver: Price of swaptions with Black model
%   function [price, delta]=Swaption_Black_receiver(t1,t2,freq,sigma_black,strike,ZC_curve)
%   5.4. IRS_approx_fix_payer: Price of IRS approximated 
%   (we approximate the PV of the 
%   IRS as a short fixed rate bond, with coupons equal to the swap rate 
%   and face value equal to the IRS notional, plus a long cash position with 
%   same amount as the IRS notional)
%   function [price, delta]=IRS_approx_fix_payer(mat, freq, S,ZC_curve)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

%% Input data
market_data = [ 0.25 4.84; 0.50 5.10; 0.75 5.27; 1.00 5.43; 
                1.25 5.32; 1.50 5.22; 1.75 5.11; 2.00 5.00; 
                2.25 4.90; 2.50 4.80; 2.75 4.70; 3.00 4.60; 
                3.25 4.54; 3.50 4.48; 3.75 4.41; 4.00 4.35; 
                4.25 4.31; 4.50 4.27; 4.75 4.23; 5.00 4.20; 
                5.25 4.17; 5.50 4.15; 5.75 4.13; 6.00 4.11; 
                6.25 4.09; 6.50 4.07; 6.75 4.05; 7.00 4.03; 
                7.25 4.02; 7.50 4.01; 7.75 3.99; 8.00 3.98; 
                8.25 3.97; 8.50 3.97; 8.75 3.96; 9.00 3.95; 
                9.25 3.94; 9.50 3.94; 9.75 3.93; 10.0 3.92];
swaption_maturity = 3.0;
swaption_tenor = 5.0;
swaption_frequency = 4;     %Quarterly fixed coupons 
sigma_black = 0.3151;       %Black swaption volatility quoted by Refinitiv
irs_maturity = 5.0;         %Maturity of the IRS
irs_fixed_coupon_freq = 4;  %Quarterly fixed coupons 
irs_coupon_rate=4.20;       %Annual coupon rate

%% Q1:Bootstrap

[ZC_curve]=ZC_bootstrap_IRS_only(market_data,irs_fixed_coupon_freq);

%% Plot curve

figure
plot(ZC_curve(:,1),ZC_curve(:,2),'x','Linewidth', 2,'LineStyle','-')
%grid
set(gca, 'YGrid', 'on', 'XGrid', 'off')
title('Zero Coupon Curve')
xlabel('Time (year frac)')

% %Chart in percentage points
% yt = get(gca, 'ytick'); 
% % convert y to [0 100]
% yt100 = yt*100; 
% % convert to string
% % note the transpose so each number is on its own line
% ytstr = num2str(yt100'); 
% % make a cell string
% ytcell = cellstr(ytstr);  
% % get rid of leading/trailing spaces
% ytcell_trim = strtrim(ytcell); 
% % append % on there
% ytpercent = strcat(ytcell_trim, '%'); 
% set(gca, 'yticklabel', ytpercent); % make the change
% 

%% Q2: Strike of the ATM swaption

t1 = swaption_maturity;
t2 = t1 + swaption_tenor;
freq = swaption_frequency;
strike = s_fwd( t1, t2, freq, ZC_curve );


%% Price and analytic delta of the swaption

[sw_price, sw_delta]=Swaption_Black_receiver(t1,t2,freq,sigma_black,strike,ZC_curve);

%% Price and analytic delta of the swap 5y

mat = irs_maturity;
freq=irs_fixed_coupon_freq;
S=irs_coupon_rate/100;
[irs5y_price, irs5y_delta]=IRS_approx_fix_payer(mat, freq, S, ZC_curve);

%% Q3: MtM

% Actual portfolio exposure
sw_N = 100000000;                 % Swaption 2x5 receiver notional
irs5y_N = 34000000;               % IRS fixed rate receiver notional

MtM = sw_N*sw_price+irs5y_N*irs5y_price;


%% Q4: analytic delta_01

DV01_approx = 1e-4*(irs5y_delta*irs5y_N+sw_delta*sw_N);

%% Q5: DV01 by numerical differentiation

% Price and greeks of the swaption with 'upward shifted curves'
irs_curve_shift = 1e-2 * ones([40,1]);
mkt_data_up(:,1) = market_data(:,1);
mkt_data_up(:,2) =  market_data(:,2) + irs_curve_shift;

% Bootstrap ans pricing under the shifted curve
[ZC_curve_up] = ZC_bootstrap_IRS_only(mkt_data_up,irs_fixed_coupon_freq);

% Delta01 of the portfolio by numeric differentiation
[irs5y_price_DV01, ~] = IRS_approx_fix_payer(mat, freq, S, ZC_curve_up);
[sw_price_DV01, ~] = Swaption_Black_receiver(t1,t2,freq,sigma_black,strike,ZC_curve_up);

DV01 = sw_N*sw_price_DV01 + irs5y_N*irs5y_price_DV01 - MtM;


%% Q6: Bucketed delta_01 by numerical differentiation

gradient = gradient_comp(t1, t2, irs5y_N, sw_N, market_data, 1e-2, MtM, irs_maturity, irs_fixed_coupon_freq, swaption_frequency, S, sigma_black, strike);
DV01_grad = 0;

for i=[2 5 10]
    DV01_grad = DV01_grad + DV01_buck(i,[2 5 10],market_data,gradient);
end

% Check the accuracy of the numerical bucketed delta_01
check = DV01 - DV01_grad;

%%
gradient = gradient_comp(t1, t2, irs5y_N, sw_N, market_data, 1e-2, MtM, irs_maturity, irs_fixed_coupon_freq, swaption_frequency, S, sigma_black, strike);

weight_2=DV01_buck(2,[2 5 10],market_data,gradient);
mkt_data_up(:,1) = market_data(:,1);
mkt_data_up(:,2) = market_data(:,2) + weight_2;

[ZC_curve_up] = ZC_bootstrap_IRS_only(mkt_data_up,irs_fixed_coupon_freq);
[irs2y_price_up, ~] = IRS_approx_fix_payer(2, freq,5/100, ZC_curve_up);
[sw_price_up, ~] = Swaption_Black_receiver(t1,t2,freq,sigma_black,strike,ZC_curve_up);

MtM_swaption = sw_N*sw_price;
MtM_up_swaption = sw_N*sw_price_up;
DV01_swaption_2y = MtM_up_swaption-MtM_swaption;


[irs2y_price, ~] = IRS_approx_fix_payer(2, freq,5/100, ZC_curve);
MtM_swap=irs2y_price;
MtM_up_swap= irs2y_price_up;
DV01_swap_2y = 1e6*(MtM_up_swap-MtM_swap)

N_s2y=-DV01_swaption_2y/DV01_swap_2y







