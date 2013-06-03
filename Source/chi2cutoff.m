function cutoff = chi2cutoff(dof, alpha)
%CHI2CUTOFF Chi-square cut-off
%   cutoff = chi2cutoff(dof, alpha)
%
%   The cut-off value cutoff is the value of chi-square with dof degrees of
%   freedom which has fraction alpha of its distribution to the right.
%
%   Mathematically: Solve alpha = chi2pvalue(cutoff,dof) for cutoff

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

cutoff = fzero(@(x)(chi2pvalue(x,dof) - alpha), 2*sqrt(2*dof)+dof);