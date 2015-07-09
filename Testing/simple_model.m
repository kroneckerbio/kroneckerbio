function [m, con, obj, opts] = simple_model(simpleopts)
% Build simple model, experiment, and fitting procedure.
% Takes optional input objectiveFun, a string for the name of the objective
%   function to test with.

% Don't warn about extraneous input arguments:
%#ok<*INUSD,

if nargin < 1
    simpleopts = struct;
end

defaultopts.objectiveFun = 'objectiveWeightedSumOfSquares';
defaultopts.steadyState = false;

simpleopts = mergestruct(defaultopts, simpleopts);

objectiveFun = simpleopts.objectiveFun;
steadyState = simpleopts.steadyState;

if nargin < 1
    objectiveFun = [];
end

if isempty(objectiveFun)
    objectiveFun = 'objectiveWeightedSumOfSquares';
end

% Build model
m = LoadModelMassAction('Simple.txt');
m = addStatesAsOutputs(m);
m = AddOutput(m, 'all_ligand', 'ligand');
m = FinalizeModel(m);

ns = m.ns;
nu = m.nu;

defaultu = m.u;

% Input
q = [1;4;3];
nq = numel(q);
border1 = 2;
input = Input(m, @u, border1, q, @dudq, @d2udq2);

% Dose
h = [2;1;3];
nh = numel(h);
border2 = 3;
border3 = 10;
dose = Dose(m, @d, 0:border3, h, @dddh, @d2ddh2);

% Experiment
if steadyState
    con = experimentSteadyState(m, [], [], input, dose, [], 'SimpleExperiment');
else
    con = experimentInitialValue(m, [], input, dose, 'SimpleExperiment');
end

% Objective
sd = sdLinear(0.1, 1);
% Picked a few values near a simulation
% output, time, value
if steadyState
    values = [
        2, 1, 1;
        4, 2, 3;
        1, 4, 5;
        10, 3, 0;
        8, 6, 0;
        9, 5, 8.5;
        2, 2, 1;
        3, 4, 0.8;
        ];
else
    values = [ 
        2, 1, 13;
        4, 2, 3;
        1, 4, 3.4;
        10, 3, 17;
        8, 6, 0;
        9, 5, 3;
        2, 2, 16;
        3, 4, 8;
        ];
end

switch objectiveFun
    case 'objectiveWeightedSumOfSquares'
%       obj = objectiveWeightedSumOfSquares(values(:,1), values(:,2), sd, values(:,3), 'SimpleData');
        obs = observationLinearWeightedSumOfSquares(values(:,1), values(:,2), sd, 'SimpleData');
        obj = obs.Objective(values(:,3));
    case 'objectiveWeightedSumOfSquaresNonNeg'
        obj = objectiveWeightedSumOfSquaresNonNeg(values(:,1), values(:,2), sd, values(:,3), 'SimpleData');
    otherwise
        error('Error:simple_model:Objective function %s not recognized.', objectiveFun)
end

% Options
opts.RelTol = 1e-6;
opts.Verbose = 0;
opts.UseParams = [1;3;5];
opts.UseSeeds = [2;3];
opts.UseInputControls = [1;2];
opts.UseDoseControls = [1;3];

opts.AbsTol = GoodAbsTol(m, con, sd, opts);

    function val = u(t,q)
        if t < 0 % Steady state simulation case
            val = defaultu;
        elseif t < border1
            val = [t; q(1)];
        else
            val = [t; q(2)*q(3)+q(3)];
        end
    end

    function val = dudq(t, q)
        if t < 0 % Steady state simulation case
            val = sparse(nu, nq);
        elseif t < border1
            val = sparse([0, 0, 0; 1, 0, 0]);
        else
            val = sparse([0, 0, 0; 0, q(3), q(2)+1]);
        end
    end

    function val = d2udq2(t, q)
        if t < 0 % Steady state simulation case
            val = sparse(nu*nq, nq);
        elseif t < border1
            val = sparse(nu*nq,nq);
        else
            val = sparse(nu*nq,nq);
            val(4,3) = 1;
            val(6,2) = 1;
        end
    end

    function val = d(t, h)
        if floor(t) == t
            if t < 0 % Steady state simulation case
                val = zeros(nh, 1);
            elseif t < border2
                val = [h(1); 0; 0];
            elseif t <= border3
                val = [h(2)*h(3)+h(3); 0; 0];
            end
        else
            val = zeros(nh, 1);
        end
    end

    function val = dddh(t, h)
        if floor(t) == t
            if t < 0 % Steady state simulation case
                val = sparse(ns, nh);
            elseif t < border2
                val = [1, 0, 0; sparse(ns-1, nh)];
            elseif t <= border3
                val = [0, h(3), h(2)+1; sparse(ns-1, nh)];
            end
        else
            val = sparse(ns, nh);
        end
    end

    function val = d2ddh2(t, h)
        if floor(t) == t
            if t < 0
                val = sparse(ns*nh, nh);
            elseif t < border2
                val = sparse(ns*nh, nh);
            elseif t <= border3
                val = sparse(ns*nh, nh);
                val(4,3) = 1;
                val(7,2) = 1;
            end
        else
            val = sparse(ns*nh, nh);
        end
    end
end
