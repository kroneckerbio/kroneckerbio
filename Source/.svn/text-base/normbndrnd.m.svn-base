function r = normbndrnd(mu, sigma, lb, ub, n)
%NORMBNDRND Draw random samples from the bounded normal distribution.
%
%   r = normbndrnd(mu, sigma, lb, ub, n)
%
%   Vector r has a length of n. Each entry in r is a random number drawn
%   from the bounded univarite normal distribution with mean mu, standard
%   deviation sigma, lower bound lb, and upper bound ub.
%
%   This function uses five different algorithms (chosen according to the
%   bounds) designed to efficiently draw random samples regardless of how
%   the bounds are positioned. It has been shown that these algorithms
%   contained are sufficient to draw from any bounded region of the normal
%   distribution.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if nargin < 5
    n = 1;
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
            expmean = 1/lb;
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
            ri = exprnd(expmean);
%            if ri < (ub - lb) && rand <= (bellcurve(ri+lb) / (c * exp(-lambda*ri)))
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