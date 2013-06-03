function m = symbolic2MassAction(SymModel, opts)

if nargin < 2
    opts = [];
end

defaultOpts.Verbose = 0;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1, 0);

%% Extract symbolic components
name    = SymModel.Name;

nv      = SymModel.nv;
nk      = SymModel.nk;
ns      = SymModel.ns;
nq      = SymModel.nq;
nu      = SymModel.nu;
nx      = SymModel.nx;
nr      = SymModel.nr;

vSyms   = SymModel.vSyms;
vNames  = SymModel.vNames;
dv      = SymModel.dv;
v       = SymModel.v;

kSyms   = SymModel.kSyms;
kNames  = SymModel.kNames;
k       = SymModel.k;

sSyms   = SymModel.sSyms;
sNames  = SymModel.sNames;
s       = SymModel.s;

qSyms   = SymModel.qSyms;
qNames  = SymModel.qNames;
q       = SymModel.q;

uSyms   = SymModel.uSyms;
uNames  = SymModel.uNames;
vuInd   = SymModel.vuInd;
u       = SymModel.u;

xSyms   = SymModel.xSyms;
xNames  = SymModel.xNames;
vxInd   = SymModel.vxInd;
x0      = SymModel.x0;

rNames  = SymModel.rNames;
r       = SymModel.r;
S       = SymModel.S;
Su      = SymModel.Su;

%% Basic constants
vStrs = cell(nv,1);
for iv = 1:nv
    vStrs{iv} = char(vSyms(iv));
end

uNamesFull = cell(nu,1);
for iu = 1:nu
    uNamesFull{iu} = [vNames{vuInd(iu)} '.' uNames{iu}];
end

xNamesFull = cell(nx,1);
for ix = 1:nx
    xNamesFull{ix} = [vNames{vxInd(ix)} '.' xNames{ix}];
end

xuSyms = [uSyms; xSyms];
vxuInd = [vuInd; vxInd];
vxuNames = vNames(vxuInd); % TODO: remove the need for this variable
xuNamesFull = [uNamesFull; xNamesFull];

%% Sanity check
assert(all(vec(fix(S) == S)) && all(vec(abs(S) <= 2)), 'KroneckerBio:symbolic2MassAction:InvalidStoichiometry', 'Invalid stoichiometry matrix does not consist entirely of -2, -1, 0, 1, and 2')

%% Initialize model
% Basic structure
m = InitializeModel(name);

%% Loop over compartments
for iv = 1:nv
    m = AddCompartment(m, vNames{iv}, dv(iv), v(iv));
end

%% Loop over kinetic parameters
for ik = 1:nk
    m = AddParameter(m, kNames{ik}, k(ik));
end

%% Loop over seed parameters
for is = 1:ns
    m = AddSeed(m, sNames{is}, s(is));
end

%% Loop over inputs
% Each input control parameter q must be used by exactly one input species
qs_used = false(nq,1);

for iu = 1:nu
    qSyms_i = symvar(u(iu));
    qSyms_i = qSyms_i(find(qSyms_i ~= sym('t'))); % Always allow t to pass % ** MATLAB bug do not replace find with logical **
    nqu_i = numel(qSyms_i);

    % Ensure none have been used by other inputs
    q_indexes = lookup(qSyms_i, qSyms);
    assert(~any(qs_used(q_indexes)), 'KroneckerBio:symbolic2MassAction:MultipurposeControl', 'An input control parameter is used by multiple inputs')
    qs_used(q_indexes) = true;
    
    % Convert the symbolics into strings
    u_func        = symbolic2string('u', nu);
    dudq     = symbolic2string('dudq', nu,nq);
    
    % Replace input control parameters with systematic names
    for i = 1:nqu_i
        qStr = char(qSyms_i(i));
        sys_string = sprintf('q(%d)', i);
        u_func    = regexprep(u, qStr, sys_string, 0);
        dudq = regexprep(dudq, qStr, sys_string, 0);
    end
    
    % Convert to a function handle
    u_func = eval(['@(t,q) [' u_func ']']);
    dudq = eval(['@(t,q) inf2big(nan2zero(sparse([' dudq '])))']);
    
    m = AddInput(m, uNames{iu}, vNames{vuInd(iu)}, {u_func; dudq}, q(q_indexes));
end

%% Loop over states
dx0ds = double(jacobian(x0, sSyms));
x0c = double(x0 - dx0ds*sSyms);

for ix = 1:nx
    % Add constant value first
    if x0c(ix)
        values_i = {'', x0c(ix)};
    else
        values_i = cell(0,1);
    end
    
    % Append seed parameters
    values_i = [values_i; 
                sNames(vec(logical(dx0ds(ix,:)))), dx0ds(ix,dx0ds(ix,:) ~= 0)];

    m = AddState(m, xNames{ix}, vNames{vxInd(ix)}, values_i);
end

%% Loop over reaction rate equations
if verbose; fprintf('Parsing reaction symbolics to kronecker reactions...'); end

% Note second-order parameters that have the volume baked in
bakedk = zeros(nk,1);

for ir = 1:nr
    % Extract all variables used in this reaction rate
    all_var = symvar(r(ir));
    
    % Variables that are rate constants; index in kSyms
    k_pos = lookup(all_var, kSyms);
    k_pos = k_pos(k_pos ~= 0);
    
    % Number of microscopic reactions in this rate
    nMic = numel(k_pos);
    
    % Count of products and reactants in stoichiometry matrix
    % Double up on species that get twice as much
    uReac = [find(Su(:,ir) < 0); find(Su(:,ir) == -2)];
    uProd = [find(Su(:,ir) > 0); find(Su(:,ir) == 2)];
    xReac = [find(S(:,ir) < 0); find(S(:,ir) == -2)];
    xProd = [find(S(:,ir) > 0); find(S(:,ir) == 2)];
    reac = [uReac; xReac + nu];
    prod = [uProd; xProd + nu];
    nReac = numel(reac);
    nProd = numel(prod);
    assert(nReac <= 2 && nProd <= 2, 'KroneckerBio:symbolic2MassAction:InvalidStoichiometry', 'Reaction rate %i has an unacceptable stoichiometry vector', ir)
    
    % Loop over each microscopic reaction and add it
    for iMic = 1:nMic
        % Differentiate wrt the rate constant to return only the species
        der_mic = diff(r(ir), kSyms(k_pos(iMic)));
        
        % Extract those species involved
        xu_var_mic = symvar(der_mic);
        
        % Ensure second order or less with species
        xu_pos_mic = lookup(xu_var_mic, xuSyms);
        %xuPosMic = xuPosMic(xuPosMic ~= 0); % TODO: I think this line is a mistake, it simply ignores nonspecies in the rate equation, which is wrong
        n_xu_mic = numel(xu_pos_mic);
        assert(all(xu_pos_mic ~= 0) && n_xu_mic <= 2, 'KroneckerBio:symbolic2MassAction:InvalidReaction', 'Reaction rate %i does not represent a quadratic ODE', ir)
        
        % Switch one zeroth, first, or second order
        switch(n_xu_mic)
            case(0)
                % Zeroth-order reaction
                assert(der_mic == -1 || der_mic == 1, 'KroneckerBio:symbolic2MassAction:ZerothOrderCheck', 'Reaction rate %i has a zeroth-order term that failed a sanity check', ir)
                
                if der_mic == 1
                    % Add zeroth order synthesis forward
                    assert(nReac == 0, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a zeroth-order term that does not match the stoichiometry', ir)
                    
                    % Switch on number of products
                    switch(nProd)
                        case(0)
                            % Pointless reaction
                            m = AddReaction(m, rNames{ir}, '', '', '', '', '', kNames{k_pos(iMic)});
                        case(1)
                            m = AddReaction(m, rNames{ir}, '', '', '', xuNamesFull{prod}, '', kNames{k_pos(iMic)});
                        case(2)
                            m = AddReaction(m, rNames{ir}, '', '', '', xuNamesFull{prod(1)}, xuNamesFull{prod(2)}, kNames{k_pos(iMic)});
                    end
                elseif der_mic == -1
                    % Add zeroth order synthesis reverse
                    assert(nProd == 0, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a zeroth-order term that does not match the stoichiometry', ir)
                    
                    % Switch on number of "reactants"
                    switch(nReac)
                        case(0)
                            % Pointless reaction
                            m = AddReaction(m, rNames{ir}, '', '', '', '', '', kNames{k_pos(iMic)});
                        case(1)
                            m = AddReaction(m, rNames{ir}, '', '', '', xuNamesFull{reac}, '', kNames{k_pos(iMic)});
                        case(2)
                            m = AddReaction(m, rNames{ir}, '', '', '', xuNamesFull{reac(1)}, xuNamesFull{reac(2)}, kNames{k_pos(iMic)});
                    end
                end
            case(1)
                % Derivative check for first-order reaction
                der_mic = diff(der_mic, xuSyms(xu_pos_mic));
                
                % Switch on actual first order or self second-order
                if der_mic == -1 || der_mic == 1
                    % True first-order reaction
                    
                    if der_mic == 1
                        % Add first order forward
                        assert(nReac == 1, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a single species term that does not match the stoichiometry', ir)
                        
                        % Switch on number of products
                        switch(nProd)
                            case(0)
                                m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic}, '', '', '', kNames{k_pos(iMic)});
                            case(1)
                                m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic}, '', xuNamesFull{prod}, '', kNames{k_pos(iMic)});
                            case(2)
                                m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic}, '', xuNamesFull{prod(1)}, xuNamesFull{prod(2)}, kNames{k_pos(iMic)});
                        end
                    elseif der_mic == -1
                        % Add first order reverse
                        assert(nProd == 1, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a single species term that does not match the stoichiometry', ir)
                        
                        % Switch on number of "reactants"
                        switch(nReac)
                            case(0)
                                m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic}, '', '', '', kNames{k_pos(iMic)});
                            case(1)
                                m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic}, '', xuNamesFull{reac}, '', kNames{k_pos(iMic)});
                            case(2)
                                m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic}, '', xuNamesFull{reac(1)}, xuNamesFull{reac(2)}, kNames{k_pos(iMic)});
                        end
                    end
                else
                    % Species reacts with itself
                    der_mic = diff(der_mic, xuSyms(xu_pos_mic));
                    allv = symvar(der_mic);
                    
                    % Process compartment volume
                    vPosMic = lookup(allv, vSyms);
                    vPosMic = vPosMic(vPosMic ~= 0);
                    if numel(vPosMic) == 0
                        % Compartment is baked in
                        vPosMic = find(strcmp(vxuNames{xu_pos_mic}, vNames), 1);
                        
                        % If we are going to change the parameter it must
                        % not be used by multiple compartments
                        assert(bakedk(k_pos(iMic)) == 0 || bakedk(k_pos(iMic)) == vPosMic || v(bakedk(k_pos(iMic))) == v(vPosMic), 'KroneckerBio:symbolic2MassAction:BakedReactionReuse', 'Reaction rate %i uses a second-order term that has the compartment volume baked into the kinetic parameter, but this parameter is shared between compartments with different volumes')
                        
                        % Claim this parameter for this compartment
                        if k_pos(iMic) == 0
                            bakedk(k_pos(iMic)) = vPosMic;
                        
                            % Replace the old parameter value
                            m = AddParameter(m, kNames{k_pos(iMic)}, k(k_pos(iMic)) * v(vPosMic));
                        end
                    elseif numel(vPosMic) == 1
                        % Compartment is written in the equation
                        der_mic = der_mic * vSyms(vPosMic);
                        assert(der_mic == 2 || der_mic == -2, 'KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', ir)
                        assert(vxuSyms(xu_pos_mic) == vSyms(vPosMic), 'KroneckerBio:symbolic2MassAction:CompartmentMembership', 'Reaction rate %i uses a compartment that neither reactant is a member of', ir)
                    else
                        error('KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', ir)
                    end
                    
                    assert(der_mic == 2 || der_mic == -2, 'KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', ir)
                    
                    if der_mic == 2
                        % Add second order forward
                        assert(nReac == 2, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a single species term that does not match the stoichiometry', ir)
                        
                        % Switch on number of products
                        switch(nProd)
                            case(0)
                                m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic}, xuNamesFull{xu_pos_mic}, '', '', kNames{k_pos(iMic)});
                            case(1)
                                m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic}, xuNamesFull{xu_pos_mic}, xuNamesFull{prod}, '', kNames{k_pos(iMic)});
                            case(2)
                                m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic}, xuNamesFull{xu_pos_mic}, xuNamesFull{prod(1)}, xuNamesFull{prod(2)}, kNames{k_pos(iMic)});
                        end
                    elseif der_mic == -2
                        % Add second order reverse
                        assert(nProd == 2, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a single species term that does not match the stoichiometry', ir)
                        
                        % Switch on number of "reactants"
                        switch(nReac)
                            case(0)
                                m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic}, xuNamesFull{xu_pos_mic}, '', '', kNames{k_pos(iMic)});
                            case(1)
                                m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic}, xuNamesFull{xu_pos_mic}, xuNamesFull{reac}, '', kNames{k_pos(iMic)});
                            case(2)
                                m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic}, xuNamesFull{xu_pos_mic}, xuNamesFull{reac(1)}, xuNamesFull{reac(2)}, kNames{k_pos(iMic)});
                        end
                    end
                end
            case(2)
                % Derivative check for second-order reaction
                der_mic = diff(der_mic, xuSyms(xu_pos_mic(1)));
                der_mic = diff(der_mic, xuSyms(xu_pos_mic(2)));
                allv = symvar(der_mic);
                
                % Process compartment volume
                vPosMic = lookup(allv, vSyms);
                vPosMic = vPosMic(vPosMic ~= 0);
                if numel(vPosMic) == 0
                    % Compartment is baked in
                    vPosMic = lookup(vxuNames(xu_pos_mic), vNames);
                    
                    % If we are going to change the parameter it must
                    % not be used by multiple compartments
                    assert(bakedk(k_pos(iMic)) == 0 || any(bakedk(k_pos(iMic)) == vPosMic) || any(v(bakedk(k_pos(iMic))) == v(vPosMic)), 'KroneckerBio:symbolic2MassAction:BakedReactionReuse', 'Reaction rate %i uses a second-order term that has the compartment volume baked into the kinetic parameter, but this parameter is shared between compartments with different volumes')
                    
                    % Claim this parameter for this compartment
                    if k_pos(iMic) == 0
                        % Use the freest compartment
                        [unused, ind] = max(d(vPosMic));
                        vPosMic = vPosMic(ind);
                        bakedk(k_pos(iMic)) = vPosMic;
                        
                        % Replace the old parameter value
                        m = AddParameter(m, kNames{k_pos(iMic)}, k(k_pos(iMic)) * v(vPosMic));
                    end
                elseif numel(vPosMic) == 1
                    % Compartment is written in the equation
                    der_mic = der_mic * vSyms(vPosMic);
                    assert(der_mic == 2 || der_mic == -2, 'KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', ir)
                    assert(vxuSyms(xu_pos_mic) == vSyms(vPosMic), 'KroneckerBio:symbolic2MassAction:CompartmentMembership', 'Reaction rate %i uses a compartment that neither reactant is a member of', ir)
                else
                    error('KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', ir)
                end
                
                assert(der_mic == 1 || der_mic == -1, 'KroneckerBio:symbolic2MassAction:SecondOrderCheck', 'Reaction rate %i has a double species term that failed a sanity check', ir)
                
                if der_mic == 1
                    % Add second order forward
                    assert(nReac == 2, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a double species term that does not match the stoichiometry', ir)
                    
                    % Switch on number of products
                    switch(nProd)
                        case(0)
                            m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic(1)}, xuNamesFull{xu_pos_mic(2)}, '', '', kNames{k_pos(iMic)});
                        case(1)
                            m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic(1)}, xuNamesFull{xu_pos_mic(2)}, xuNamesFull{prod}, '', kNames{k_pos(iMic)});
                        case(2)
                            m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic(1)}, xuNamesFull{xu_pos_mic(2)}, xuNamesFull{prod(1)}, xuNamesFull{prod(2)}, kNames{k_pos(iMic)});
                    end
                elseif der_mic == -1
                    % Add second order reverse
                    assert(nProd == 2, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a double species term that does not match the stoichiometry', ir)
                    
                    % Switch on number of "reactants"
                    switch(nReac)
                        case(0)
                            m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic(1)}, xuNamesFull{xu_pos_mic(2)}, '', '', kNames{k_pos(iMic)});
                        case(1)
                            m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic(1)}, xuNamesFull{xu_pos_mic(2)}, xuNamesFull{reac}, '', kNames{k_pos(iMic)});
                        case(2)
                            m = AddReaction(m, rNames{ir}, '', xuNamesFull{xu_pos_mic(1)}, xuNamesFull{xu_pos_mic(2)}, xuNamesFull{reac(1)}, xuNamesFull{reac(2)}, kNames{k_pos(iMic)});
                    end
                end
                
        end
    end
end
if verbose; fprintf('done.\n'); end

%% Finalize model
m = FinalizeModel(m);

%% %%%%%%%%%%%%%%%%%%%%%%%
%%%% Helper functions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
    function string_rep = symbolic2string(variable_name, varargin)
        dimensions = [varargin{:}];
        if numel(dimensions) == 1
            dimensions = [dimensions, 1];
        end
        
        if any(dimensions == 0)
            string_rep = ['zeros([' num2str(dimensions) '])'];
        else
            string_rep = evalc(['disp(' variable_name ')']);
        end
    end

end