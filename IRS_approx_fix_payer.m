function [price, delta]=IRS_approx_fix_payer(mat, freq, S, ZC_curve)
% Computes the approximated price of the IRS (we approximate the PV of the 
% IRS as a short fixed rate bond, with coupons equal to the swap rate 
% and face value equal to the IRS notional, plus a long cash position with 
% same amount as the IRS notional) 

% INPUT:
%
% mat:          IRS maturity
% freq:         frequency of coupons
% S:            IRS coupon
% ZC_curve:     Table of ZC rates (cont. comp. 30/360)
%               Maturities are year fractions
%
% OUTPUT:
%
% price:        price of IRS
% delta:        delta of IRS


B = exp(-ZC_curve(:,2).*ZC_curve(:,1));

c = S/freq*ones(length(B(1:freq*mat)),1);
c(end) = c(end) + 1;

% Price computation
price = -sum(c.*B(1:freq*mat)) + 1;

t = (1/freq:1/freq:mat)';
delta_y = 1; 
MacD = -sum(c.*t.*B(1:mat*freq))/sum(c.*B(1:mat*freq));

% Delta computation (in bps)
delta = - sum(c.*B(1:mat*freq)) .* delta_y .* MacD; % in bps

end