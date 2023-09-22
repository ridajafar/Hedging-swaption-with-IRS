function [DV01_buck] = DV01_buck(bucket,all_buckets,market_data,gradient)
% DV01 of the bucket
% 
% INPUTS:
% bucket:       bucket in which to compute DV01 bucket
% all_buckets:  buckets
% market_data:  Column #1: IRS maturity (year frac)
%               Column #2: MID rate
% gradient:     gradient of the price
%
% OUTPUTS:
% DV01_buck:    bucket DV01

% function handle to compute weights for the bucket
if bucket == all_buckets(1)
    f_weights = @(t) max(0, (t > bucket).*(1-(t-bucket)./(all_buckets(2)-bucket)) + (t <= bucket) );
else
    if bucket == all_buckets(end)
        f_weights = @(t) max(0, (t < bucket).*((t-all_buckets(end-1))./(bucket-all_buckets(end-1))) + (t >= bucket) );
    else
        index = find(all_buckets == bucket);
        f_weights = @(t) max(0,(t >= bucket).*(1-(t-bucket)./(all_buckets(index+1)-bucket)) + (t < bucket).*((t-all_buckets(index-1))./(bucket-all_buckets(index-1))));
    end
end

% Weights computation
weights = f_weights(market_data(:,1));

% DV01 bucket computation
DV01_buck = sum(weights.*gradient);

end