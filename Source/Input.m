function inp = Input(m, u, discontinuities, q, dudq, d2udq2)
%Input Construct an input object
%
%   inp = Input(m, u, discontinuities, q, dudq, d2udq2)
%
%   Inputs
%   m: [ Model ]
%       A KroneckerBio Model for which this input will be applied. The
%       input structure is valid for any Model which has the same number of
%       inputs
%   u: [ handle @(t,q) returns nonnegative matrix nu by numel(t) ]
%       A function handle that returns the value of the inputs at any time
%       t with associated control parameters q. This function should accept
%       t as a vector.
%   discontinuities: [ nonnegative vector {[]} ]
%       Any discontinuous times in the u function must be listed here in
%       order to ensure successful evaluation u.
%   q: [ real vector nq {[]} ]
%       The values of the input control parameters
%   dudq: [ handle @(t,q) returns real matrix nu by nq ]
%       A function handle that returns the derivative of each input with
%       respect to each input control parameter. This does not have to be
%       supplied but some Kronecker features will be disabled.
%   d2udq2: [ handle @(t,q) returns real matrix nu*nq by nq ]
%       Partial derivative of dudq with respect to q at time t
%
%   Outputs
%   inp: [ input struct scalar ]
%       A valid input structure
%
%   If q is not provided or empty, appropriate empty values will be
%   supplied for the derivatives. If q is specified and nonempty, functions
%   requiring these derivatives will crash.

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if nargin < 6
    d2udq2 = [];
    if nargin < 5
        dudq = [];
        if nargin < 4
            q = [];
            if nargin < 3
                discontinuities = [];
            end
        end
    end
end

% m
assert(is(m, 'Model'), 'KroneckerBio:Dose:m', 'm must be a Model')
m = keepfields(m, {'Type', 'u', 'nu'});

% u
assert(isfunction(u) && nargin(u) == 2, 'KroneckerBio:Input:u', 'u must be a function handle accepting 2 arguments')

% discontinuities
assert(isnumeric(discontinuities) && all(discontinuities >= 0), 'KroneckerBio:Input:discontinuities', 'discontinuities must be a nonegative vector')
discontinuities = vec(discontinuities);

% q
assert(isnumeric(q) && all(q >= 0), 'KroneckerBio:Input:q', 'q must be a nonegative vector')
q = vec(q);
nq = numel(q);

% Derivatives
assert(isfunction(dudq) && nargin(dudq) == 2, 'KroneckerBio:Input:dudq', 'dudq must be a function handle accepting 2 arguments')
assert(isfunction(d2udq2) && nargin(d2udq2) == 2, 'KroneckerBio:Input:d2udq2', 'd2udq2 must be a function handle accepting 2 arguments')

% Build input
inp.Type = 'Input';
inp.u = u;
inp.q = q;
inp.nq = nq;
inp.dudq = dudq;
inp.d2udq2 = d2udq2;
inp.discontinuities = discontinuities;
inp.Update = @update;

    function inp_out = update(q)
        assert(numel(q) == nq, 'KroneckerBio:Input:Update:q', 'q must be a vector of length nq')
        inp_out = Input(m, u, discontinuities, q, dudq, d2udq2);
    end
end
