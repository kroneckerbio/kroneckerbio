function [m, con, obj, opts] = simple_model(objectiveFun)
% Build simple model, experiment, and fitting procedure.
% Takes optional input objectiveFun, a string for the name of the objective
%   function to test with.

if nargin < 1
    objectiveFun = [];
end

if isempty(objectiveFun)
    objectiveFun = 'objectiveWeightedSumOfSquares';
end

% Build model
m = LoadModelMassAction('Simple.txt');
m = addStatesAsOutputs(m);
m = FinalizeModel(m);

ns = m.ns;
nu = m.nu;

% Input
q = [1;4;3];
nq = numel(q);
border1 = 2;
input = Input(m, @u, border1, q, @dudq, @d2udq2);

% Dose
h = [2;1;3];
nh = numel(h);
border2 = 3;
% dose = Dose(m, @d, border2, h, @dddh, @d2ddh2);
dose = doseConstant(m, 3:5, 1:6);

% Experiment
tF = 6;
con = InitialValueExperiment(m, [], input, dose, 'SimpleExperiment');

% Objective
sd = sdLinear(0.1, 1);
values = [ % Picked a few values near a simulation
    2, 1, 13;
    4, 2, 3;
    1, 4, 3.4;
    10, 3, 17;
    8, 6, 0;
    9, 5, 3;
    2, 2, 16;
    3, 4, 8;
    ];

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
        if t < border1
            val = [t; q(1)];
        else
            val = [t; q(2)*q(3)+q(3)];
        end
    end

    function val = dudq(t, q)
        if t < border1
            val = sparse([0, 0, 0; 1, 0, 0]);
        else
            val = sparse([0, 0, 0; 0, q(3), q(2)+1]);
        end
    end

    function val = d2udq2(t, q)
        if t < border1
            val = sparse(nu*nq,nq);
        else
            val = sparse(nu*nq,nq);
            val(4,3) = 1;
            val(6,2) = 1;
        end
    end

    function val = d(t, h)
        if floor(t) == t
            if t < border2
                val = [h(1); 0; 0];
            else
                val = [h(2)*h(3)+h(3); 0; 0];
            end
        else
            val = zeros(nh,1);
        end
    end

    function val = dddh(t, h)
        if floor(t) == t
            if t < border2
                val = [1, 0, 0; sparse(ns-1,nh)];
            else
                val = [0, h(3), h(2)+1; sparse(ns-1,nh)];
            end
        else
            val = sparse(ns,nh);
        end
    end

    function val = d2ddh2(t, h)
        if floor(t) == t
            if t < border2
                val = sparse(ns*nh,nh);
            else
                val = sparse(ns*nh,nh);
                val(4,3) = 1;
                val(7,2) = 1;
            end
        else
            val = sparse(ns*nh,nh);
        end
    end
end
