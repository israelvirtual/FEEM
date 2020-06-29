function [ m ] = nonneg( m )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Converts matrix elements with values < 0 or NaN to 0

I = find(isnan(m));
m(I) = zeros(1, length(I));
I = find(m <= eps);
m(I) = zeros(1, length(I));

end

