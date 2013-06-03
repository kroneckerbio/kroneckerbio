function A = nan2zero(A)
%nan2zero turns all the nans in matrix A into zeros

A(isnan(A)) = 0;