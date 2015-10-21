function r = normbndrnd(mu, sigma, lb, ub, n)
%NORMBNDRND Draw random samples from the bounded normal distribution. 
%   Each entry in r is a random number drawn
%   from the bounded univariate normal distribution with mean mu, standard
%   deviation sigma, lower bound lb, and upper bound ub. mu, sigma, lb, ub, and
%   n can be scalars or vectors of length m, where vectors must be the same dimensions and
%   scalars are expanded to match the vector length.
%
%   r = normbndrnd(mu, sigma, lb, ub, n)
%
%   m is the length of the longest vector of mu, sigma, lb, ub, and n. The
%   lengths of the arguments must be equal to m or 1.
%
% Inputs:
%   mu [ scalar double | 1 x m vector of doubles ]
%       Mean(s)
%   sigma [ scalar double | 1 x m vector of doubles ]
%       Standard deviation(s)
%   lb [ scalar double | 1 x m vector of doubles ]
%       Lower bound(s)
%   ub [ scalar double | 1 x m vector of doubles ]
%       Upper bound(s)
%   n [ scalar positive integer {1} | 1 x m vector of scalar positive integers ]
%       Desired number of random numbers to generate
%
% Outputs:
%   r [ n x m double matrix | 1 x m cell vector of ni x 1 double vectors ]
%       Random number drawn from corresponding mu, sigma, lb, and ub.
%       If n = 1 (the default), then r is a vector matching the length
%       of the inputs. If n is a non-unity scalar or a vector with the same
%       non-unity values, r is a matrix with rows sharing the corresponding n from
%       the r vector. If n is a vector with different values, r is a cell vector
%       with each element holding the n values for the corresponding n.
%
%   This function uses five different algorithms (chosen according to the
%   bounds) designed to efficiently draw random samples regardless of how
%   the bounds are positioned. It has been shown that these algorithms
%   contained are sufficient to draw from any bounded region of the normal
%   distribution.

% (c) 2015 Kevin Shi, David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if nargin < 5
    n = 1;
end

% Make sure n is valid
assert(all(n > 0 & n == floor(n)), 'normbndrnd:invalidn', 'n is not a vector of positive integers')

% Standardize to row vectors
mu    = vec(mu)';
sigma = vec(sigma)';
lb    = vec(lb)';
ub    = vec(ub)';
n     = vec(n)';

nmu    = length(mu);
nsigma = length(sigma);
nlb    = length(lb);
nub    = length(ub);
nn     = length(n);

% See if there is a vector of different n values
nDifferent = length(unique(n)) > 1;

% Make sure all arg vector lengths are consistent
lengths = [nmu, nsigma, nlb, nub, nn];
m = max(lengths);
lengthsValid = lengths == 1 | lengths == m;
assert(all(lengthsValid), 'normbndrnd:invalidArgs', 'Lengths of args not consistent')

% Standardize vectors of args
mu    = repmat(mu, m/nmu, 1);
sigma = repmat(sigma, m/nsigma, 1);
lb    = repmat(lb, m/nlb, 1);
ub    = repmat(ub, m/nub, 1);
n     = repmat(n, m/nn, 1);

% Preallocate output
if nDifferent
    r = cell(1, nn);
else % each entry has the same number of values
    r = zeros(n(1), nn);
end

% Generate random values for single set of args
for im = 1:m
    ri = normbndrndsingle(mu(im), sigma(im), lb(im), ub(im), n(im));
    if nDifferent
        r{im} = ri;
    else
        r(:,im) = ri;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Main Function %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

    function r = normbndrndsingle(mu, sigma, lb, ub, n)
        
        % Special case of sigma exactly equals zero
        if sigma == 0
            r = repmat(mu, n, 1);
            return
        end
        
        % Standardize the bounds
        lb = (lb - mu) / sigma;
        ub = (ub - mu) / sigma;
        
        % Computes the truncated normal distribution
        
        twoSideCutoff = 5;
        ratioCutoff = 1e52;
        expCutoff = 5;
        
        r = zeros(n,1);
        
        if lb >= 0 || ub <= 0
            % lb and ub are on the same side
            if ub <= 0
                % Reduce possible cases by making bounds larger than mu
                temp = lb;
                lb = -ub;
                ub = -temp;
                flip = true;
            else
                flip = false;
            end
            
            pdfRatio = exp(-0.5*(lb^2-ub^2));
            if pdfRatio < ratioCutoff
                for i = 1:n
                    r(i) = sideuniformreject();
                end
            else
                if lb < expCutoff
                    for i = 1:n
                        r(i) = halfnormalreject();
                    end
                else
                    %lambda = 0.5*(lb+sqrt(lb^2+4));
                    %c = exp(0.5*lambda^2 - lambda*lb);
                    %expmean = 1/lb;
                    for i = 1:n
                        r(i) = expreject();
                    end
                end
            end
        else
            % 0 is a member of lb:ub
            flip = false;
            if ub <= twoSideCutoff || lb <= -twoSideCutoff
                for i = 1:n
                    r(i) = normalreject();
                end
            else
                for i = 1:n
                    r(i) = miduniformreject();
                end
            end
        end
        
        if flip
            r = -r;
        end
        
        r = r*sigma + mu;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Helper Functions %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ri = normalreject()
            while true %dowhile
                ri = randn;
                if ri > lb && ri < ub
                    break
                end
            end
        end
        
        function ri = halfnormalreject()
            while true %dowhile
                ri = abs(randn);
                if ri > lb && ri < ub
                    break
                end
            end
        end
        
        function ri = miduniformreject()
            while true %dowhile
                ri = rand * (ub - lb) + lb;
                if rand <= bellcurve(ri);
                    break
                end
            end
        end
        
        function ri = sideuniformreject()
            while true %dowhile
                ri = rand * (ub - lb) + lb;
                if rand <= exp(-0.5*(ri^2-lb^2))
                    break
                end
            end
        end
        
        function ri = expreject()
            while true
                ri = -log(rand) / lb;
                %ri = exprnd(expmean); % dependent on Matlab statistics toolbox
                %if ri < (ub - lb) && rand <= (bellcurve(ri+lb) / (c * exp(-lambda*ri)))
                if ri <= (ub - lb) && rand <= bellcurve(ri)
                    break
                end
            end
            ri = ri + lb;
        end
        
        function y = bellcurve(x)
            y = exp(-0.5 .* x.^2);
        end
        
    end

end