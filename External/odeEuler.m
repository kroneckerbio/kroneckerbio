function [t, x] = odeEuler(der, t, ic, opts, u)

step = 0.1;
t0 = min(t);
tF = max(t);

nStep = ceil(tF / step);
nx = numel(ic);

tAll = zeros(1,nStep);
tAll(1) = t0;
xAll = zeros(nx,nStep);
xAll(:,1) = ic;

nonnegative = false(nx,1);
nonnegative(opts.NonNegative) = true;

for i = 2:nStep
    tAll(i)   = tAll(i-1) + step;
    xAll(:,i) = xAll(:,i-1) + der(tAll(i-1), xAll(:,i-1), u) * step;
    
    % Negatives are 0
    xAll(xAll(:,i) < 0 & nonnegative,i) = 0;
end

t = t.';
x = piecewiselinear(t, tAll, xAll.');