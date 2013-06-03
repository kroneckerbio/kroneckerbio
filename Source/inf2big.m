function A = inf2big(A)
%inf2big converts infinities to just big numbers
%
%   A = inf2big(A)

A(A < -1e6) = -1e6;
A(A > 1e6)  =  1e6;