function [ par_yield_fwd ] = s_fwd( t1, t2, freq, ZC_curve )
% Computes the IRS forward rate

% INPUT: 
% t1:               swaption maturity
% t2:               swaption end time (maturity + tenor)    
% freq:             frequency of coupons
% ZC_curve:         Table of ZC rates (cont. comp. 30/360)
%                   Maturities are year fractions
%
% OUTPUT:
% par_yield_fwd:    forward swap rate 

t = t1:1/freq:t2;

% Yearly Rates and Discounts needed
r = zeros(length(t),1);
B = zeros(length(t),1);
for i=1:length(t)
    r(i) = ZC_curve(ZC_curve(:,1)==t(i),2);
    B(i) = exp(-r(i)*t(i));
end

% Compute Forward swap rate 
par_yield_fwd = (1-B(end)/B(1))/(sum(B(2:end)/B(1))/freq);

end