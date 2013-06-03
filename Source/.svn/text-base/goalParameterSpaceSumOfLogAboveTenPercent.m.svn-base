function G = goalParameterSpaceSumOfLogAboveTenPercent(F)

F = (F + F.')/2; % Symmetrize the FIM

lambdas = eig(F); % Compute eigenvalues

lambdas(lambdas < eps) = eps; % Floor ultra-small directions to eps

errors = 1 ./ sqrt(lambdas); % Compute parameter direction uncertainty

errors(errors < 0.1) = 0.1; % Floor errors at 10%

G = sum(log(errors));