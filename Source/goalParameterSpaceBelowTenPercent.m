function G = goalParameterSpaceBelowTenPercent(F)

nT = size(F, 1);

% Compute eigenvalues
lambdas = eig(F);

% Floor ultra-small directions to eps
lambdas(lambdas < eps) = eps;

% Compute cutoff
lambdacut = 100; % 1/(0.1)^2; that is, cutoff for ten percent

% Subtract 1 for each eigenvalue above threshhold
G = -sum(lambdas > lambdacut);

% Subtract a fraction that cannot sum above 1 to break ties
G = G - sum(lambdas(lambdas <= lambdacut) ./ (lambdacut*nT));
