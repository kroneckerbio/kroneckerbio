function con = SteadyStateExperiment(m, s, inp, dos, time_scale, name)
%SteadyStateExperiment Construct a KroneckerBio experimental conditions
%   structure describing a initial value problem first run to steady state
%
%   con = SteadyStateExperiment(m, tF, s, inp, dos, name)
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model for which these experiments will be run
%   s: [ nonnegative vector ns ]
%       Default = m.s
%       The values of the seed parameters
%   inp: [ input struct scalar | handle @(t) returns nonegative vector nu |
%          nonnegative vector nu ]
%       Default = m.u
%       The definition of the input species values
%   dos: [ dose struct scalar ]
%       Default = doseZero(m)
%       The definition of the dose amounts and schedule
%   time_scale: [ nonnegative scalar ]
%       Default = 10
%       The typical time scale for an observation. The steady state is
%       determined to be reached when the states are expected to change
%       less than the tolerance over this time scale.
%   name: [ string ]
%       Default = ''
%       An arbitrary name for the experiment
%
%   Outputs
%   con: [ experiment struct scalar ]
%       The KroneckerBio experimental conditions structure
%
%   For the meanings of the fields of con see "help experimentZero"
% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 6
    name = [];
    if nargin < 5
        time_scale = [];
        if nargin < 4
            dos = [];
            if nargin < 3
                inp = [];
                if nargin < 2
                    s = [];
                end
            end
        end
    end
end

if isempty(time_scale)
    time_scale = 10;
end
if isempty(s)
    s = m.s;
end
if isempty(inp)
    inp = inputConstant(m, m.u);
end
if isempty(dos)
    dos = doseZero(m);
end
if isempty(name)
    name = '';
end

% m
assert(isscalar(m) && is(m, 'Model'), 'KroneckerBio:Experiment:m', 'm must be a Model')
m = keepfields(m, {'Type', 's', 'u', 'ns', 'nu'});
nu = m.nu;

% s
assert(numel(s) == m.ns, 'KroneckerBio:Experiment:s', 's must a vector with length equal to m.ns')
s = vec(s);

% inp
if isnumeric(inp)
    assert(numel(inp) == m.nu, 'KroneckerBio:Experiment:inp', 'inp, when numeric, must have a length of m.nu')
    inp = inputConstant(m, inp);
end

assert(is(inp, 'Input'), 'KroneckerBio:Experiment:inp', 'inp must be an Input')

% dos
assert(is(dos, 'Dose'), 'KroneckerBio:Experiment:dos', 'dos must be a Dose')

% name
assert(ischar(name), 'KroneckerBio:Experiment:name', 'name must be a string')

% Build experiment
con.Type = 'Experiment:SteadyState';
con.Name = name;
con.nu = m.nu;
con.ns = m.ns;
con.nq = numel(inp.q);
con.nh = numel(dos.h);
con.s  = s;
con.q  = inp.q;
con.h  = dos.h;
[con.u,con.dudq,con.d2udq2] = getU(inp,m.nu);
con.d  = @(t)dos.d(t,dos.h);
con.dddh = @(t)dos.dddh(t,dos.h);
con.d2ddh2 = @(t)dos.d2ddh2(t,dos.h);
con.inp = inp;
con.dos = dos;
con.SteadyState = true;
con.Periodic = false;
con.Discontinuities = vec(unique([inp.discontinuities; dos.discontinuities]));
con.Update = @update;
con.private.TimeScale = time_scale;

    function con_out = update(s, q, h)
        con_out = SteadyStateExperiment(m, s, inp.Update(q), dos.Update(h), time_scale, name);
    end

end

function [u_t,dudq_t,d2udq2_t] = getU(inp,nu)

u_tq = inp.u;
dudq_tq = inp.dudq;
d2udq2_tq = inp.d2udq2;
q = inp.q;
clear inp
u_t = @(t) u_tq(t,q);

% Determine if u is vectorized, and fix if not
testut = u_t([1 2]);
if size(testut,2) == 1
    u_t = @ut_vectorized;
elseif size(testut,2) ~= 2
    error('u should return an nu-by-1 or nu-by-nt vector of input values')
end

dudq_t = @(t) dudq_tq(t,q);
d2udq2_t = @(t) d2udq2_tq(t,q);

    function u = ut_vectorized(t)
        nt = numel(t);
        u = zeros(nu,nt);
        for ti = 1:nt
            u(:,ti) = u_tq(t(ti),q);
        end
    end

end