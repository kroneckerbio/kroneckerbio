function dos = doseConstant(m, amounts, times, receivers)
%doseConstant Easily create dosing functions for a series of constant
%   doses
%
%   dos = dosingConstant(m, amounts, times, receivers)
%
%   This returns a dose structure with a dose applied in a given amount at
%   a given time. The amount of the dose will be a dose control parameter.
%   Set UseDoseControls to [] to disable optimization of these dose
%   amounts.
%
%   Inputs
%   m: [ Model ]
%       A KroneckerBio Model for which this dose will be applied. The dose
%       structure is valid for any Model which has the same number of
%       seeds
%   amounts: [ nonegative vector ]
%       The amount of each dose. If a dose is applied to different seeds,
%       this will be a vector.
%   times: [ nonegative vector ]
%       Times at which a dose will be given
%   receivers: [ index vector {1:numel(amount)} ]
%       If some seeds are not to be dosed, use this vector to indicate
%       which seeds are receiving the dose. If it is not provided, this
%       function will assume that the dose is given to the first
%       numel(amount) seeds.
%
%   Outputs
%   dos
%       A valid dose structure

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if nargin < 4
    receivers = [];
end

if isempty(receivers)
    receivers = vec(1:numel(amounts));
end

% m
assert(is(m, 'Model'), 'KroneckerBio:doseConstant:m', 'm must be a Model')
m = keepfields(m, {'Type', 'ns'});

% amounts
assert(isnumeric(amounts) && all(amounts >= 0), 'KroneckerBio:doseConstant:amounts', 'amounts must be a nonegative vector')
amounts = vec(amounts);
nh = numel(amounts);

% times
assert(isnumeric(times) && all(times >= 0), 'KroneckerBio:doseConstant:times', 'times must be a nonegative vector')
times = vec(unique2013a(times));

% receivers
ns = m.ns;
assert(isindex(receivers, ns), 'KroneckerBio:doseConstant:receviers', 'receivers must be an index vector into m.s')
receivers = vec(receivers);

dos = Dose(m, @d, times, amounts, @dddh, @d2ddh2);

    function val = d(t, h)
        val = zeros(ns,1);
        if any(t == times)
            val(receivers) = h;
        end
    end

    function val = dddh(t, h)
        val = sparse(ns,nh);
        if any(t == times)
            val(receivers,:) = eye(nh);
        end
    end

    function val = d2ddh2(t, h)
        val = sparse(ns*nh,nh);
    end
end
