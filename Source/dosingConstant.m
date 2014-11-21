function [d, dddq] = dosingConstant(ns, whichSeeds, whichTimes, whichControls)
%dosingConstant Easily create dosing functions for a series of constant
%   doses
%
%   [d, dddq] = dosingConstant(ns, whichSeeds, whichTimes, whichControls)
%
%   This returns a dosing function and its derivative (ready for plopping
%   into and experiment) where the dose is equal to q (part of the
%   experiment). The dose is given at each timepoint in whichTimes. The
%   seeds increased are given by whichSeeds and the control parameter that
%   applies to it is given by whichControls. Don't forget to set the value
%   of the dose (q) on the experiment!
%
%   Inputs
%   ns: [ nonnegative scalar ]
%       The number of seeds in the corresponding model. This is only used
%       to ensure that the output is the right size
%   whichSeeds: [ index vector ]
%       To which seeds will this dose be applied
%   whichTimes: [ nonnegative vector ]
%       At which times will this dose be applied
%   whichControls [ index vector {1:numel(whichSeeds)} ]
%       Optional: Which control parameters contain the value of the dose. This must
%       be the same length as whichSeeds as there will be a one-to-one
%       mapping of each control parameter indexed by whichControls to each
%       seed indexed by whichSeeds. If whichControls is missing, this
%       function assumes that the first values in q are the doses. 
%
%   Outputs
%   d: [ @(t,q) function handle returns nonnegative vector ns by 1 ]
%       Ready dosing function
%   dddq: [ @(t,q) function handle returns nonnegative vector ns by nq ]
%       Ready derivative of dosing function

% (c) 2014 David R Hagen
% This work is released under the MIT license.

if nargin < 4
    whichControls = vec(1:numel(whichSeeds));
end

whichTimes = vec(whichTimes);

d = @uniformDose;
dddq = @uniformDoseDerivative;

    function s = uniformDose(t,q)
        s = zeros(ns,1);
        if any(t == whichTimes)
            s(whichSeeds) = q(whichControls);
        end
    end

    function dsdq = uniformDoseDerivative(t,q)
        dsdq = zeros(ns,numel(q));
        if any(t == whichTimes)
            dsdq(whichSeeds,whichControls) = eye(numel(q));
        end
    end
end
