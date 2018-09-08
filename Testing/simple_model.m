function [m, con, obj, opts, con_sscheck] = simple_model(simpleopts)
% Build simple model, experiment, and fitting procedure.
% Takes optional input objectiveFun, a string for the name of the objective
%   function to test with.

% Don't warn about extraneous input arguments:
%#ok<*INUSD,

if nargin < 1
    simpleopts = struct;
end

defaultopts.objectiveFun = 'observationLinearWeightedSumOfSquares';
defaultopts.steadyState = false;
defaultopts.kineticPrior = false;
defaultopts.logKineticPrior = false;
defaultopts.logSeedPrior = false;

simpleopts = mergestruct(defaultopts, simpleopts);

objectiveFun = simpleopts.objectiveFun;
steadyState = simpleopts.steadyState;
kineticPrior = simpleopts.kineticPrior;
logKineticPrior = simpleopts.logKineticPrior;
logSeedPrior = simpleopts.logSeedPrior;

if nargin < 1
    objectiveFun = [];
end

if isempty(objectiveFun)
    objectiveFun = 'observationLinearWeightedSumOfSquares';
end

% Build model
m = LoadModelMassAction('Simple.txt');
m = addStatesAsOutputs(m);
m = AddOutput(m, 'all_ligand', 'ligand'); % interpreted as regex?
m = FinalizeModel(m);

ns = m.ns;
nu = m.nu;

defaultu = m.u;

% Input
q = [1;4;3];
nq = numel(q);
border1 = 2;
input = Input(m, @u, border1, q, @dudq, @d2udq2);

% Basal input
basalborder1 = 3;
basalinput = Input(m, @basalu, basalborder1, q, @basaldudq, @basald2udq2);

% Dose
h = [2;1;3];
nh = numel(h);
border2 = 3;
border3 = 10;
dose = Dose(m, @d, 0:border3, h, @dddh, @d2ddh2);

% Experiment
if steadyState
    con = experimentSteadyState(m, [], basalinput, input, dose, [], 'SimpleExperiment');
    con_sscheck = experimentInitialValue(m, [], basalinput, [], 'SimpleExperiment');
else
    con = experimentInitialValue(m, [], input, dose, 'SimpleExperiment');
    con_sscheck = [];
end

% Objective
sd = sdLinear(0.1, 1);
% Picked a few values near a simulation
if steadyState
    
    % Initialize expected value array
    values = [
        % output, time, value
        2,  1, 0;
        4,  2, 0;
        1,  4, 0;
        10, 3, 0;
        8,  6, 0;
        9,  5, 0;
        2,  2, 0;
        3,  4, 0;
        ];
    
    % Set parameters
    kobj = [
        1.0965
        1.1320
        1.9421
        1.9561
        1.5752
        1.0598
        1.2348
        1.3532
        1.8212
        1.0154
        ];
    
    sobj = [
        1.0430
        1.1690
        1.6491
        ];
    
    qobj = [
        1.7317
        1.6477
        1.4509
        ];
    
    hobj = [
        1.5470
        1.2963
        1.7447
        ];
    
    % Simulate with new parameters to get expected values
    mobj = m.Update(kobj);
    conobj = con.Update(sobj, qobj, hobj);
    simobjopts.Verbose = 0;
    simobj = SimulateSystem(mobj, conobj, max(values(:,2)), simobjopts);
    
    % Fill in expected values with reasonable starting values close to correct values
    for vi = 1:size(values,1)
        val = simobj.y(values(vi,2), values(vi,1));
        values(vi,3) = val + normrnd(0, 0.01*val);
    end
    
else
    values = [ 
        2,  1, 13;
        4,  2, 3;
        1,  4, 3.4;
        10, 3, 17;
        8,  6, 0;
        9,  5, 3;
        2,  2, 16;
        3,  4, 8;
        ];
end

switch objectiveFun
    case 'observationLinearWeightedSumOfSquares'
        obs = observationLinearWeightedSumOfSquares(values(:,1), values(:,2), sd, 'SimpleData');
        obj = obs.Objective(values(:,3));
    case 'observationLogWeightedSumOfSquares'
        obs = observationLogWeightedSumOfSquares(values(:,1), values(:,2), sd, 'SimpleData');
        obj = obs.Objective(values(:,3));
    otherwise
        error('Error:simple_model:Objective function %s not recognized.', objectiveFun)
end

% Also include priors, if requested
if kineticPrior
    kprior = 1+rand(m.nk,1);
    Vlogkprior = rand(m.nk);
    Vlogkprior = Vlogkprior*Vlogkprior.'; % ensure symmetric
    obj_kprior = objectiveNormalPriorOnKineticParameters(kprior, Vlogkprior, 'SimpleKPrior');
    obj = [obj; obj_kprior];
end
if logKineticPrior
    kprior = 1+rand(m.nk,1);
    Vlogkprior = rand(m.nk);
    Vlogkprior = Vlogkprior*Vlogkprior.'; % ensure symmetric
    Vkprior = diag(kprior)*Vlogkprior*diag(kprior); % Multiply squared parameters in because objective function takes "non-normalized" inputs
    obj_logkprior = objectiveLogNormalPriorOnKineticParameters(kprior, Vkprior, 'SimpleLogKPrior');
    obj = [obj; obj_logkprior]; 
end
if logSeedPrior
    sprior = 1+rand(m.ns,1);
    Vlogsprior = rand(m.ns);
    Vlogsprior = Vlogsprior*Vlogsprior.'; % ensure symmetric
    Vsprior = diag(sprior)*Vlogsprior*diag(sprior); % Multiply squared parameters in because objective function takes "non-normalized" inputs
    obj_logsprior = objectiveLogNormalPriorOnSeedParameters(sprior, Vsprior, 'SimpleLogSPrior');
    obj = [obj; obj_logsprior]; 
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

    function val = basalu(t,q)
        if t < basalborder1
            val = t/basalborder1*defaultu/10*q(1).^2;
        else
            val = defaultu/10*q(1).^2;
        end
    end

    function val = dudq(t, q)
        if t < border1
            val = sparse([0, 0, 0; 1, 0, 0]);
        else
            val = sparse([0, 0, 0; 0, q(3), q(2)+1]);
        end
    end

    function val = basaldudq(t, q)
        val = zeros(nu,nq);
        if t < basalborder1
            val(:,1) = 2*t/basalborder1*defaultu/10*q(1);
        else
            val(:,1) = 2*defaultu/10*q(1);
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

    function val = basald2udq2(t, q)
        val = sparse(nu*nq, nq);
        linearinds = sub2ind([nu nq nq], (1:nu)', ones(nu,1), ones(nu,1));
        if t < basalborder1
            val(linearinds) = 2*t/basalborder1*defaultu/10;
        else
            val(linearinds) = 2*defaultu/10;
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
            else
                val = zeros(nh,1);
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
