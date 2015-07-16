function [m, con, obj, opts, eve] = michaelis_menten_model()
% Simple analytic model for a variety of tests

%% Build model
m = InitializeModelAnalytic('MichaelisMentenModel');

m = AddCompartment(m, 'solution', 3, 1);

m = AddInput(m, 'E', 'solution', 1);

m = AddState(m, 'S', 'solution', 'S0^2');
m = AddState(m, 'P', 'solution', 10);

m = AddParameter(m, 'Km', 10);
m = AddParameter(m, 'kcat', 2);

m = AddSeed(m, 'S0', 5);

m = AddReaction(m, 'rxn1', 'S', '', 'P', '', 'kcat*E*S/(Km+S)');

m = AddOutput(m, 'S', 'S');
m = AddOutput(m, 'P', 'P');
m = AddOutput(m, 'r', 'kcat*E*S/(Km+S)');

m = FinalizeModel(m);

%% Make experimenal condition and generate test data
if nargout > 1

    dos = doseConstant(m, 3, [2; 4; 6; 8]);
    u = @(t,q) repmat(q.^2,1,numel(t));
    dudq = @(t,q) 2.*q;
    d2udq2 = @(t,q) 2;
    inp = Input(m, u, [], 1, dudq, d2udq2); 
    con = experimentInitialValue(m, [], inp, dos);

    sd = sdLinear(0.1, 1);
    % Values approximately from simulation for k = [15;10], other parameters the same
    values = [
        1   1       24
        1   2.5     30
        1   4       37 
        2   3       14
        2   5       17
        2   8.5     23
        3   1.5     1.4
        3   3       1.5
        3   4.5     1.6
        ];
    obs = observationLinearWeightedSumOfSquares(values(:,1), values(:,2), sd, 'MichaelisMentenData');
    obj = obs.Objective(values(:,3));

    opts.Verbose = false;
    opts.RelTol = 1e-9;
    opts.Verbose = 0;
    opts.UseParams = [1;2];
    opts.UseSeeds = [1];
    opts.UseInputControls = [1];
    opts.UseDoseControls = [1];

    opts.AbsTol = GoodAbsTol(m, con, sd, opts);
    
    eve1 = eventDropsBelow(m, 3, 1);
    eve2 = eventDropsBelow(m, 1, 2);
    eve = [eve1;eve2];
    
end