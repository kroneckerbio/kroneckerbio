function [ratio, dfdx_finite, dfdx_analytic] = fdiff(x0, f, dfdx)
% Compare an analytical derivative with a finite difference derivative
diff = 1e-7;

f0 = vec(f(x0));
nf = numel(f0);
nx = numel(x0);

dfdx_analytic = reshape(dfdx(x0), [nf,nx]);

dfdx_finite = zeros(nf,nx);
for i = 1:nx
    x = x0;
    x(i) = x(i) + diff;
    dfdx_finite(:,i) = (vec(f(x)) - f0) ./ diff;
end

dfdx_finite = dfdx_finite;
ratio = dfdx_finite ./ dfdx_analytic;
