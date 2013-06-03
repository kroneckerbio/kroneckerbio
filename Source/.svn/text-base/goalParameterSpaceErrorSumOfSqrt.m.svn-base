function G = goalParameterSpaceErrorSumOfSqrt(F)

F = (F + F.')/2; % Symmetrize the FIM

lambdas = eig(F); % Compute eigenvalues

lambdas(lambdas < eps) = eps; % Floor ultra-small directions to eps

errors = 1 ./ sqrt(lambdas); % Compute parameter direction uncertainty

G = sum(sqrt(errors));