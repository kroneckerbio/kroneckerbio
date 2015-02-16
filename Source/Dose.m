function dos = Dose(m, d, schedule, h, dddh, d2ddh2)
%Dosing Construct a dose object
%
%   dos = Dose(m, d, discontinuities, h, dddh, d2dh2)
%
%   Inputs
%   m: [ Model ]
%       A KroneckerBio Model for which this dose will be applied. The dose
%       structure is valid for any Model which has the same number of seeds
%   d: [ handle @(t,h) returns nonnegative matrix ns by numel(t) ]
%       A function that returns the values of the doses at any time t. Doses are applied via the seeds
%       to instantaneously update the state of the system. Must defined for
%       all t and return zeros(ns,1) for times at which a dose is not
%       given.
%   schedule: [ nonnegative vector ]
%       Any dosing times in the d function must be listed here in order to
%       ensure successful evaluation of d.
%   h: [ real vector nq {[]} ]
%       The values of the dose control parameters
%   dddh: [ handle @(t,h) returns real matrix ns by nh {[]} ]
%       Partial derivative of d with respect to the dose control parameters
%       h at time t
%   d2ddh2: [ handle @(t,h) returns ns*nh by nh {[]} ]
%       Partial derivative of dddh with respect to h at time t
%
%   Outputs
%   dos: [ dose struct scalar ]
%       A valid dose structure
%
%   If h is not provided or empty, appropriate empty values will be
%   supplied for the derivatives. If h is specified and nonempty, functions
%   requiring these derivatives will crash.

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if nargin < 6
    d2ddh2 = [];
    if nargin < 5
        dddh = [];
        if nargin < 4
            h = [];
        end
    end
end

% m
assert(is(m, 'Model'), 'KroneckerBio:Dose:m', 'm must be a Model')
m = keepfields(m, {'Type', 'ns'});

% d
assert(isfunction(d) && nargin(d) == 2, 'KroneckerBio:Dose:d', 'd must be a function handle acceptiong 2 arguments')

% discontinuities
assert(isnumeric(schedule) && all(schedule >= 0), 'KroneckerBio:Dosing:discontinuities', 'discontinuities must be a nonegative vector')
schedule = vec(unique(schedule));

% h
assert(isnumeric(h) && all(h > 0), 'KroneckerBio:Dosing:h', 'h must be a nonegative vector')
h = vec(h);
nh = numel(h);

% Derivatives
assert(isfunction(dddh) && nargin(dddh) == 2, 'KroneckerBio:Dosing:dddh', 'dddh must be a function handle acceptiong 2 arguments')
assert(isfunction(d2ddh2) && nargin(d2ddh2) == 2, 'KroneckerBio:Dosing:d2dh2', 'd2dh2 must be a function handle acceptiong 2 arguments')

% Build Dose
dos.Type = 'Dose';
dos.d = d;
dos.h = h;
dos.nh = nh;
dos.dddh = dddh;
dos.d2ddh2 = d2ddh2;
dos.discontinuities = schedule;
dos.Update = @update;

    function dos_out = update(h)
        assert(numel(h) == nh, 'KroneckerBio:Dose:Update:h', 'h must be a vector of length nh')
        dos_out = Dose(m, d, schedule, h, dddh, d2ddh2);
    end
end
