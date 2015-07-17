function [m_kron, con, obj, opts, eve] = michaelis_menten_model()

% Build symbolic model
syms Km kcat S0 En S P

m.kNames = {'kcat';'Km'};
m.kSyms = [kcat;Km];
m.k = [10;2];

m.sNames = {'S0'};
m.sSyms = S0;
m.s = 5;

m.uNames = {'En'};
m.uSyms = En;
m.u = 1;
m.nu = length(m.u);

m.xNames = {'S';'P'};
m.xSyms = [S;P];
m.x0 = [S0.^2;10];
m.nx = length(m.xSyms);

r = kcat*En*S/(Km+S);
m.f = [-r;r];

m.vNames = {'solution'};
m.vSyms = sym('solution');
m.v = 1;
m.nv = 1;
m.dv = 3;
m.vxInd = ones(m.nx,1);
m.vuInd = ones(m.nu,1);

m.yNames = {};
m.yStrings = {};
m.y = sym([]);
yNames = {'S','P','r'};
yStrings = {'S','P','kcat*En*S/(Km+S)'};
m = AddOutputsToSymbolic(m, yNames, yStrings);

% Convert symbolic model to analytic model
m_kron = symbolic2PseudoKronecker(m);

if nargout > 1

    dos = doseConstant(m_kron, 3, [2; 4; 6; 8]);
    u = @(t,q) repmat(q.^2,1,numel(t));
    dudq = @(t,q) 2.*q;
    d2udq2 = @(t,q) 2;
    inp = Input(m_kron, u, [], 1, dudq, d2udq2); 
    con = experimentInitialValue(m_kron, [], inp, dos);

    sd = sdLinear(0.1, 1);
    % Values approximately from simulation for k = [15;10], other parameters the same
    values = [
        1   1       15
        1   2.5     12
        1   4       3 
        2   3       35
        2   5       48
        2   8.5     64
        3   1.5     11
        3   3       9.6
        3   4.5     9.5
        ];
    obs = observationLinearWeightedSumOfSquares(values(:,1), values(:,2), sd, 'MichaelisMentenData');
    obj = obs.Objective(values(:,3));

    opts.Verbose = false;
    opts.RelTol = 1e-6;
    opts.Verbose = 0;
    opts.UseParams = [1;2];
    opts.UseSeeds = [1];
    opts.UseInputControls = [1];
    opts.UseDoseControls = [1];

    opts.AbsTol = GoodAbsTol(m_kron, con, sd, opts);
    
    eve1 = eventDropsBelow(m, 3, 1);
    eve2 = eventDropsBelow(m, 1, 2);
    eve = [eve1;eve2];
    
end

end