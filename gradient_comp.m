function [gradient] = gradient_comp(t1, t2, irs_N, sw_N, market_data, shift, MtM, irs_maturity, irs_fixed_coupon_freq, swaption_frequency, S, sigma_black, strike)
% Approximated gradient of the price
%
% INPUTS:
% 
% t1:                       swaption maturity
% t2:                       swaption maturity + tenor
% irs_N:                    notional of the IRS
% sw_N:                     notional of the swaption
% market_data:              Market data
% shift:                    shift in the rates as percentage
% MtM:                      mark to market
% irs_maturity:             IRS maturity
% irs_fixed_coupon_freq:    frequency of IRS fixed coupon payments
% swaption_frequency:       frequency of swaption payments
% S:                        IRS coupon rate 
% sigma_black:              sigma 
% strike:                   strike 
%
% OUTPUTS:
%
% gradient:                 gradient of the price


gradient = zeros(length(market_data),1);

for i = 1:length(market_data)
    % Computnig increments
    incr = zeros(size(market_data));
    incr(i,2) = shift;

    % Peforming Bootstrap
    [ZC_curve_shift] = ZC_bootstrap_IRS_only(market_data + incr, 4);

    % Computing new prices
    [sw_price_shift, ~] = Swaption_Black_receiver(t1,t2,swaption_frequency,sigma_black,strike,ZC_curve_shift);
    [irs_price_shift, ~] = IRS_approx_fix_payer(irs_maturity, irs_fixed_coupon_freq, S, ZC_curve_shift);
    
    % Mark-to-Market with the new prices
    MtM_shift = irs_price_shift*irs_N + sw_price_shift*sw_N;

    % Approx gradient computation 
    gradient(i) = 1e-4 / (shift/100) * (MtM_shift - MtM);

end

end