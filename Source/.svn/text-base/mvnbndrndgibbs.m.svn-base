function r = mvnbndrndgibbs(mu, V, lb, ub, r0, n, burn, skip)
%MVNBNDRNDGIBBS Sample from the bounded multivariate normal distribution
%   using a Gibbs sampler
%
%   r = mvnbndrndgibbs(mu, V, lb, ub, r0, n, burn, skip)
%
%   Matrix r is length(mu) by n. Each column is a random sample of the
%   bounded multivariate normal distribution with mean mu, variance V,
%   lower bounds lb, and upper bounds ub. Vectors mu, lb, and ub are of the
%   same length. V is a square matrix length(mu) by length(mu). Vector r0
%   indicates the starting location for the Gibbs sampler. Scalar burn
%   indicates the number of iterations to discard before saving samples.
%   Scalar skip indicates how many iterations to discard between saving
%   samples.
%
%   This Gibbs sampler takes 2*m steps per iteration. 1*m steps are
%   along the directions of the bounds (in random order). 1*m steps are
%   along the eigendirections of the variance (in random order). Each step
%   is a bounded univariate normal distribution.
%
%   According to the theory of Gibbs samplers, the distribution of samples
%   will eventually match the desired distribution. However, serial
%   correlation among the samples almost always is always a problem. This
%   particular implementation has been designed to bypass the strongest
%   sources of correlation and has been shown in practice to be very
%   reliable.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 8
    skip =[];
    if nargin < 7
        burn = [];
        if nargin < 6
            n = [];
            if nargin < 5
                r0 = [];
                if nargin < 4
                    ub = [];
                    if nargin < 3
                        lb = [];
                        if nargin < 2
                            V = [];
                        end
                    end
                end
            end
        end
    end
end

% Constants
m = numel(mu);

% Default inputs
if isempty(V)
    V = eye(m);
end
if isempty(lb)
    lb = -inf(m,1);
end
if isempty(ub)
    ub = inf(m,1);
end
if isempty(r0)
    r0 = mu;
end
if isempty(n)
    n = 1;
end
if isempty(burn)
    burn = 0;
end
if isempty(skip)
    skip = 0;
end

% Decompose V
Vinv = inv(V);
h2 = diag(Vinv).^-1;
c = zeros(m,m-1);
for i = 1:m
    c(i,:) = -h2(i) * Vinv(i,[1:i-1 i+1:m]);
end
h = sqrt(h2);
L = chol(V, 'lower');
Linv = inv(L);

% Zero mean
ylb = lb - mu;
yub = ub - mu;
y = r0 - mu;

% Stack bounds
B = [-eye(m); eye(m)];
b = [-lb; ub];
yb = b - B*mu;
zB = B*L;

% Burn as specified
for i = 1:burn
    y = boxIter(y);
    y = mvnIter(y);
end

% Main loop
r = zeros(m,n);
for i = 1:n
    for j = 1:skip+1
        y = boxIter(y);
        y = mvnIter(y); 
    end
    
    r(:,i) = mu + y; % Unzero mean
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Helper functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function y = boxIter(y)
        for k = randperm(m)
            condmean = c(k,:) * vec(y([1:k-1, k+1:m]));
            y(k) = normbndrnd(condmean, h(k), ylb(k), yub(k));
            %y(k) = randraw('normaltrunc', [ylb(k), yub(k), condmean, h(k)]);
        end
    end

    function y = mvnIter(y)
        z = Linv*y;
        for k = randperm(m)
            mk = [1:k-1, k+1:m]; % minus k: directions held constant
            zbk = yb - zB(:,mk)*vec(z(mk)); % conditional bounds for zk
            zlb = max(zbk(zB(:,k) < 0) ./ zB(zB(:,k) < 0,k));
            zub = min(zbk(zB(:,k) > 0) ./ zB(zB(:,k) > 0,k));
            z(k) = normbndrnd(0, 1, zlb, zub);
            %z(k) = randraw('normaltrunc', [zlb, zub, 0, 1]);
        end
        y = L*z;
    end

end