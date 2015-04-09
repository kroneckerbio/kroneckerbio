function con = SteadyStateExperiment(m, s, inp, dos, time_scale, name)
%SteadyStateExperiment Construct a KroneckerBio experimental conditions
%   structure describing a initial value problem first run to steady state
%
%   con = Experiment(m, tF, s, steadyState, periodic, u, d,
%                    discontinuities, q, dudq, dddq, name)
%
%   The inputs to this function allow one to set all the variables that are
%   permitted on a KroneckerBio experimental conditions structure.
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
if nargin < 5
    name = [];
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
con.Type = 'Experiment:InitialValue';
con.Name = name;
con.nu = m.nu;
con.ns = m.ns;
con.nq = numel(inp.q);
con.nh = numel(dos.h);
con.s  = s;
con.q  = inp.q;
con.h  = dos.h;
con.u  = @(t)inp.u(t,inp.q);
try
    % See if u is vectorized
    test = con.u([1,2]);
    assert(size(test) == [nu,2])
catch
    % Nope, vectorize it
    con.u = @(t)vectorize_u(t, @(ti)inp.u(ti,inp.q));
end
con.dudq = @(t)inp.dudq(t,inp.q);
con.d2udq2 = @(t)inp.d2udq2(t,inp.q);
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

    function val = vectorize_u(t, f_u)
        nt = numel(t);
        val = zeros(nu,nt);
        for i = 1:nt
            val(:,i) = f_u(t(i));
        end
    end

    function con_out = update(s, q, h)
        con_out = InitialValueExperiment(m, s, inp.Update(q), dos.Update(h), name);
    end
end
