function [price, delta] = Swaption_Black_receiver(t1,t2,freq,sigma_black,strike,ZC_curve)
% Computes the price and analytic delta of the swaption

% INPUT:
%
% t1:               swaption maturity
% t2:               swaption end time (maturity + tenor)    
% freq:             frequency of coupons
% sigma_black:      sigma of the swaption
% strike:           strike of the swaption
% ZC_curve:         Table of ZC rates (cont. comp. 30/360)
%                   Maturities are year fractions
%
% OUTPUT:
%
% price:            price of the swaption
% delta:            delta of the swaption

% S0 computation 
par_yield_fwd = s_fwd( t1, t2, freq, ZC_curve );

% Price computation
d1 = log(par_yield_fwd/strike)/(sigma_black*sqrt(t1)) + 0.5*sigma_black*sqrt(t1);
d2 = log(par_yield_fwd/strike)/(sigma_black*sqrt(t1)) - 0.5*sigma_black*sqrt(t1);

B = exp(-ZC_curve(t1*freq:t2*freq,2).*ZC_curve(t1*freq:t2*freq,1));
BPV = sum(B(2:end)/B(1))/freq;

price = B(1)*BPV*(strike*cdf('Normal',-d2,0,1) - par_yield_fwd*cdf('Normal',-d1,0,1)); % Receiver

% Delta computation (in bps)
delta = -B(1)*BPV*cdf('Normal',-d1,0,1); % Receiver

end