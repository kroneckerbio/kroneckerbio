function [ratio, dfdx_finite, dfdx_analytic] = fdiff(x0, f, dfdx)
% Compare an analytical derivative with a finite difference derivative
diff = 1e-8;

f0 = vec(f(x0));
nf = numel(f0);
nx = numel(x0);

dfdx_analytic = reshape(dfdx(x0), [nf,nx]);

dfdx_finite = zeros(nf,nx);
for i = 1:nx
    x = x0;
    x(i) = x(i) + diff*1i;
    dfdx_finite(:,i) = imag(vec(f(x))) ./ diff;
end

ratio = dfdx_finite ./ dfdx_analytic;
