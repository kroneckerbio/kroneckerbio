function dos = doseList(m, amounts, times, receivers)
%doseList Easily create dosing functions for a series of doses at various
%   times to various
%
%   dos = doseList(m, amounts, times, receivers)
%
%   This returns a dose structure with each amount listed appleid at each
%   time listed to each seed listed. There are no dose control parameters
%   associated with this structure.
%
%   Inputs
%   m: [ Model ]
%       A KroneckerBio Model for which this dose will be applied. The dose
%       structure is valid for any Model which has the same number of
%       seeds
%   amounts: [ nonegative vector ]
%       Amount of each dose.
%   times: [ nonegative vector ]
%       Time at which a dose will be given.
%   receivers: [ index vector {1:numel(amount)} ]
%       The seed to whic hthe dose will be applied.
%
%   Outputs
%   dos
%       A valid dose structure

% (c) 2016 David R Hagen
% This work is released under the MIT license.

if nargin < 4
    receivers = [];
end

if isempty(receivers)
    receivers = ones(size(amounts));
end

% m
assert(is(m, 'Model'), 'KroneckerBio:doseList:m', 'm must be a Model')
m = keepfields(m, {'Type', 'ns'});

% amounts
assert(isnumeric(amounts) && all(amounts >= 0), 'KroneckerBio:doseList:amounts', 'amounts must be a nonegative vector')
amounts = vec(amounts);

% times
assert(isnumeric(times) && all(times >= 0), 'KroneckerBio:doseList:times', 'times must be a nonegative vector')

% receivers
ns = m.ns;
assert(isindex(receivers, ns), 'KroneckerBio:doseList:receviers', 'receivers must be an index vector into m.s')
receivers = vec(receivers);

% Replicate scalar inputs
if isscalar(amounts)
    amounts = zeros(size(times)) + amounts;
end
if isscalar(receivers)
    receivers = zeros(size(times)) + receivers;
end

dos = Dose(m, @d, times);

    function val = d(t, h)
        inds = times == t;
        val = sparse(receivers(inds), ones(nnz(inds),1), amounts(inds), ns,1);
    end
end
