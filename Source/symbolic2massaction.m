function m = symbolic2massaction(symbolic, opts)
% Inputs:
%   symbolic [ Model.Symbolic struct ]
%   opts [ options struct ]
% Outputs:
%   m [ Model.MassActionAmount struct ]
%
% Note: mass action kroneckerbio models now require unique species names.
%   TODO: implement converter that fixes models that rely on unique
%   compartment.species

if nargin < 2
    opts = [];
end

% Default options
opts_.Verbose = 0;

opts = mergestruct(opts_, opts);

verbose = logical(opts.Verbose);

% Sanity check
if ~strcmpi(symbolic.Type, 'Model.Symbolic')
    error('symbolic2analytic: input model is not a Model.Symbolic')
end

%% Fail fast sanity checks to make sure model is mass action (incomplete)
% No rules - warn if rules are found and ignore them
if symbolic.nz > 0
    warning('symbolic2massaction: Rules found, ignoring.')
end

% Only up to 2 reactants and products per reaction allowed
nr     = symbolic.nr;
rNames = symbolic.rNames;
r      = symbolic.r;
for i = 1:nr
    reactants = r{i,1};
    products = r{i,2};
    
    if length(reactants) > 2 || length(products) > 2
        error('symbolic2massaction: reaction %s has > 2 reactants or products', rNames)
    end
end

% Assume state initial conditions are linear combinations of seeds
%   Always true if converting from SBML or SimBio because IC's are identically seeds

% Check that reaction rates are valid later

%% Add symbolic components to model
if verbose; fprintf('Adding components to mass action model...'); end

% Initialize model
m = InitializeModelMassActionAmount(symbolic.Name);

% Compartments
nv     = symbolic.nv;
vNames = symbolic.vNames;
vIDs   = symbolic.vIDs;
v      = symbolic.v; % doubles
dv     = symbolic.dv;
for i = 1:nv
    m = AddCompartment(m, vNames{i}, dv(i), v(i));
end

% States
nx      = symbolic.nx;
xNames  = symbolic.xNames;
xIDs    = symbolic.xIDs;
xvNames = symbolic.xvNames;
x       = symbolic.x; % cell array of strings
for i = 1:nx
    m = AddState(m, xNames{i}, xvNames{i}, strrep(x{i}, '"', '')); % remove quotes since mass action finalizer needs it
end

% Inputs
nu      = symbolic.nu;
uNames  = symbolic.uNames;
uIDs    = symbolic.uIDs;
uvNames = symbolic.uvNames;
u       = symbolic.u; % doubles
for i = 1:nu
    m = AddInput(m, uNames{i}, uvNames{i}, u(i));
end

% Seeds
ns      = symbolic.ns;
sNames  = symbolic.sNames;
sIDs    = symbolic.sIDs;
s       = symbolic.s; % doubles
for i = 1:ns
    m = AddSeed(m, sNames{i}, s(i));
end

% Parameters
nk     = symbolic.nk;
kNames = symbolic.kNames;
kIDs   = symbolic.kIDs;
k      = symbolic.k; % doubles
for i = 1:nk
    m = AddParameter(m, kNames{i}, k(i));
end

if verbose; fprintf('done.\n'); end

%% Reactions
if verbose; fprintf('Parsing and validating mass action rate forms...'); end
nr     = symbolic.nr;
rNames = symbolic.rNames;
rIDs   = symbolic.rIDs;
r      = symbolic.r; % cell matrix

% Get names/ids for symbolic substitutions
xuNames  = [xNames; uNames];
xuvNames = [xvNames; uvNames];
xuIDs    = [xIDs; uIDs];
allNames = [xuNames; vNames; kNames; sNames];
allIDs   = [xuIDs; vIDs; kIDs; sIDs];

% Full species names
% Note/TODO: currently hacked to just use short names
% xuNamesFull = strcat(vxuNames, '.', xuNames); % compartment.species
xuNamesFull = xuNames;

% Make symbolic variables
vSyms = sym(vIDs);
xuSyms = sym(xuIDs);
kSyms = sym(kIDs);

% Note second-order parameters that have the volume baked in
kBaked = zeros(nk,1);

% Process reactions
for i = 1:nr
    
    name = rNames{i};
    
    reactants = r{i,1};
    products = r{i,2};
    rate = r{i,3};
    
    nReac = length(reactants);
    nProd = length(products);
    
    % Convert names to ids and symbolic vars so they can be manipulated symbolically
    for j = 1:nReac
        reactants{j} = name2id(quoteInvalid(reactants{j}), xuNames, xuIDs, xuvNames);
    end
    for j = 1:nProd
        products{j} = name2id(quoteInvalid(products{j}), xuNames, xuIDs, xuvNames);
    end
    reactants = sym(reactants);
    products = sym(products);
    
    rate = sym(name2id(rate, allNames, allIDs, xuvNames));
    rate = evaluateExternalFunctions(rate, allIDs); % for resolving "power" function
    
    % Extract variables in the rate - Note: vars cell array of strings is in
    %   alphabetic order
    vars = symvar(rate);
    
    % Get positions of rate params and species
    kPos = lookup(vars, kSyms);
    kPos(kPos==0) = [];
    reac = lookup(reactants, xuSyms);
    prod = lookup(products, xuSyms);
    
    % Get number of microscopic reactions in this rate (1 or 2)
    nMic = length(kPos);
    
    % Loop over each microscopic reaction and add it
    for iMic = 1:nMic
        % Differentiate wrt the rate constant to return only the species
        derMic = diff(rate, kSyms(kPos(iMic)));
        
        % Extract those species involved
        xuMicSyms = symvar(derMic);
        
        % Ensure second order or less with species
        xuMicPos = lookup(xuMicSyms, xuSyms);
        %xuPosMic = xuPosMic(xuPosMic ~= 0); % TODO: I think this line is a mistake, it simply ignores nonspecies in the rate equation, which is wrong
        nxuMic = numel(xuMicPos);
        assert(all(xuMicPos ~= 0) && nxuMic <= 2, 'KroneckerBio:symbolic2MassAction:InvalidReaction', 'Reaction rate %i does not represent a quadratic ODE', i)
        
        % Switch one zeroth, first, or second order
        switch(nxuMic)
            case(0)
                % Zeroth-order reaction
                assert(derMic == -1 || derMic == 1, 'KroneckerBio:symbolic2MassAction:ZerothOrderCheck', 'Reaction rate %i has a zeroth-order term that failed a sanity check', i)
                
                if derMic == 1
                    % Add zeroth order synthesis forward
                    assert(nReac == 0, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a zeroth-order term that does not match the stoichiometry', i)
                    
                    % Switch on number of products
                    switch(nProd)
                        case(0)
                            % Pointless reaction
                            m = AddReaction(m, name, '', '', '', '', kNames{kPos(iMic)});
                        case(1)
                            m = AddReaction(m, name, '', '', xuNamesFull{prod}, '', kNames{kPos(iMic)});
                        case(2)
                            m = AddReaction(m, name, '', '', xuNamesFull{prod(1)}, xuNamesFull{prod(2)}, kNames{kPos(iMic)});
                    end
                elseif derMic == -1
                    % Add zeroth order synthesis reverse
                    assert(nProd == 0, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a zeroth-order term that does not match the stoichiometry', i)
                    
                    % Switch on number of "reactants"
                    switch(nReac)
                        case(0)
                            % Pointless reaction
                            m = AddReaction(m, name, '', '', '', '', kNames{kPos(iMic)});
                        case(1)
                            m = AddReaction(m, name, '', '', xuNamesFull{reac}, '', kNames{kPos(iMic)});
                        case(2)
                            m = AddReaction(m, name, '', '', xuNamesFull{reac(1)}, xuNamesFull{reac(2)}, kNames{kPos(iMic)});
                    end
                end
            case(1)
                % Derivative check for first-order reaction
                derMic = diff(derMic, xuSyms(xuMicPos));
                
                % Switch on actual first order or self second-order
                if derMic == -1 || derMic == 1
                    % True first-order reaction
                    
                    if derMic == 1
                        % Add first order forward
                        assert(nReac == 1, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a single species term that does not match the stoichiometry', i)
                        
                        % Switch on number of products
                        switch(nProd)
                            case(0)
                                m = AddReaction(m, name, xuNamesFull{xuMicPos}, '', '', '', kNames{kPos(iMic)});
                            case(1)
                                m = AddReaction(m, name, xuNamesFull{xuMicPos}, '', xuNamesFull{prod}, '', kNames{kPos(iMic)});
                            case(2)
                                m = AddReaction(m, name, xuNamesFull{xuMicPos}, '', xuNamesFull{prod(1)}, xuNamesFull{prod(2)}, kNames{kPos(iMic)});
                        end
                    elseif derMic == -1
                        % Add first order reverse
                        assert(nProd == 1, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a single species term that does not match the stoichiometry', i)
                        
                        % Switch on number of "reactants"
                        switch(nReac)
                            case(0)
                                m = AddReaction(m, name, xuNamesFull{xuMicPos}, '', '', '', kNames{kPos(iMic)});
                            case(1)
                                m = AddReaction(m, name, xuNamesFull{xuMicPos}, '', xuNamesFull{reac}, '', kNames{kPos(iMic)});
                            case(2)
                                m = AddReaction(m, name, xuNamesFull{xuMicPos}, '', xuNamesFull{reac(1)}, xuNamesFull{reac(2)}, kNames{kPos(iMic)});
                        end
                    end
                else
                    % Species reacts with itself
                    derMic = diff(derMic, xuSyms(xuMicPos));
                    allv = symvar(derMic);
                    
                    % Process compartment volume
                    vPosMic = lookup(allv, vSyms);
                    vPosMic = vPosMic(vPosMic ~= 0);
                    if numel(vPosMic) == 0
                        % Compartment is baked in
                        vPosMic = find(strcmp(xuvNames{xuMicPos}, vNames), 1);
                        
                        % If we are going to change the parameter it must
                        % not be used by multiple compartments
                        assert(kBaked(kPos(iMic)) == 0 || kBaked(kPos(iMic)) == vPosMic || v(kBaked(kPos(iMic))) == v(vPosMic), 'KroneckerBio:symbolic2MassAction:BakedReactionReuse', 'Reaction rate %i uses a second-order term that has the compartment volume baked into the kinetic parameter, but this parameter is shared between compartments with different volumes')
                        
                        % Claim this parameter for this compartment
                        if kPos(iMic) == 0
                            kBaked(kPos(iMic)) = vPosMic;
                            
                            % Replace the old parameter value
                            m = AddParameter(m, kNames{kPos(iMic)}, k(kPos(iMic)) * v(vPosMic));
                        end
                    elseif numel(vPosMic) == 1
                        % Compartment is written in the equation
                        derMic = derMic * vSyms(vPosMic);
                        assert(derMic == 2 || derMic == -2, 'KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', i)
                        assert(vxuSyms(xuMicPos) == vSyms(vPosMic), 'KroneckerBio:symbolic2MassAction:CompartmentMembership', 'Reaction rate %i uses a compartment that neither reactant is a member of', i)
                    else
                        error('KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', i)
                    end
                    
                    assert(derMic == 2 || derMic == -2, 'KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', i)
                    
                    if derMic == 2
                        % Add second order forward
                        assert(nReac == 2, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a single species term that does not match the stoichiometry', i)
                        
                        % Switch on number of products
                        switch(nProd)
                            case(0)
                                m = AddReaction(m, name, xuNamesFull{xuMicPos}, xuNamesFull{xuMicPos}, '', '', kNames{kPos(iMic)});
                            case(1)
                                m = AddReaction(m, name, xuNamesFull{xuMicPos}, xuNamesFull{xuMicPos}, xuNamesFull{prod}, '', kNames{kPos(iMic)});
                            case(2)
                                m = AddReaction(m, name, xuNamesFull{xuMicPos}, xuNamesFull{xuMicPos}, xuNamesFull{prod(1)}, xuNamesFull{prod(2)}, kNames{kPos(iMic)});
                        end
                    elseif derMic == -2
                        % Add second order reverse
                        assert(nProd == 2, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a single species term that does not match the stoichiometry', i)
                        
                        % Switch on number of "reactants"
                        switch(nReac)
                            case(0)
                                m = AddReaction(m, name, xuNamesFull{xuMicPos}, xuNamesFull{xuMicPos}, '', '', kNames{kPos(iMic)});
                            case(1)
                                m = AddReaction(m, name, xuNamesFull{xuMicPos}, xuNamesFull{xuMicPos}, xuNamesFull{reac}, '', kNames{kPos(iMic)});
                            case(2)
                                m = AddReaction(m, name, xuNamesFull{xuMicPos}, xuNamesFull{xuMicPos}, xuNamesFull{reac(1)}, xuNamesFull{reac(2)}, kNames{kPos(iMic)});
                        end
                    end
                end
            case(2)
                % Derivative check for second-order reaction
                derMic = diff(derMic, xuSyms(xuMicPos(1)));
                derMic = diff(derMic, xuSyms(xuMicPos(2)));
                allv = symvar(derMic);
                
                % Process compartment volume
                vPosMic = lookup(allv, vSyms);
                vPosMic = vPosMic(vPosMic ~= 0);
                if numel(vPosMic) == 0
                    % Compartment is baked in
                    vPosMic = lookup(xuvNames(xuMicPos), vNames);
                    
                    % If we are going to change the parameter it must
                    % not be used by multiple compartments
                    assert(kBaked(kPos(iMic)) == 0 || any(kBaked(kPos(iMic)) == vPosMic) || any(v(kBaked(kPos(iMic))) == v(vPosMic)), 'KroneckerBio:symbolic2MassAction:BakedReactionReuse', 'Reaction rate %i uses a second-order term that has the compartment volume baked into the kinetic parameter, but this parameter is shared between compartments with different volumes')
                    
                    % Claim this parameter for this compartment
                    if kPos(iMic) == 0
                        % Use the freest compartment
                        [~, ind] = max(d(vPosMic));
                        vPosMic = vPosMic(ind);
                        kBaked(kPos(iMic)) = vPosMic;
                        
                        % Replace the old parameter value
                        m = AddParameter(m, kNames{kPos(iMic)}, k(kPos(iMic)) * v(vPosMic));
                    end
                elseif numel(vPosMic) == 1
                    % Compartment is written in the equation
                    derMic = derMic * vSyms(vPosMic);
                    assert(derMic == 2 || derMic == -2, 'KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', i)
                    assert(vxuSyms(xuMicPos) == vSyms(vPosMic), 'KroneckerBio:symbolic2MassAction:CompartmentMembership', 'Reaction rate %i uses a compartment that neither reactant is a member of', i)
                else
                    error('KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', i)
                end
                
                assert(derMic == 1 || derMic == -1, 'KroneckerBio:symbolic2MassAction:SecondOrderCheck', 'Reaction rate %i has a double species term that failed a sanity check', i)
                
                if derMic == 1
                    % Add second order forward
                    assert(nReac == 2, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a double species term that does not match the stoichiometry', i)
                    
                    % Switch on number of products
                    switch(nProd)
                        case(0)
                            m = AddReaction(m, name, xuNamesFull{xuMicPos(1)}, xuNamesFull{xuMicPos(2)}, '', '', kNames{kPos(iMic)});
                        case(1)
                            m = AddReaction(m, name, xuNamesFull{xuMicPos(1)}, xuNamesFull{xuMicPos(2)}, xuNamesFull{prod}, '', kNames{kPos(iMic)});
                        case(2)
                            m = AddReaction(m, name, xuNamesFull{xuMicPos(1)}, xuNamesFull{xuMicPos(2)}, xuNamesFull{prod(1)}, xuNamesFull{prod(2)}, kNames{kPos(iMic)});
                    end
                elseif derMic == -1
                    % Add second order reverse
                    assert(nProd == 2, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a double species term that does not match the stoichiometry', i)
                    
                    % Switch on number of "reactants"
                    switch(nReac)
                        case(0)
                            m = AddReaction(m, name, xuNamesFull{xuMicPos(1)}, xuNamesFull{xuMicPos(2)}, '', '', kNames{kPos(iMic)});
                        case(1)
                            m = AddReaction(m, name, xuNamesFull{xuMicPos(1)}, xuNamesFull{xuMicPos(2)}, xuNamesFull{reac}, '', kNames{kPos(iMic)});
                        case(2)
                            m = AddReaction(m, name, xuNamesFull{xuMicPos(1)}, xuNamesFull{xuMicPos(2)}, xuNamesFull{reac(1)}, xuNamesFull{reac(2)}, kNames{kPos(iMic)});
                    end
                end
                
        end
        
    end
end

if verbose; fprintf('done.\n'); end


end

function expr = quoteInvalid(expr)
% Wrap expression (single identifier) in quotes if it contains invalid
% characters ('\W')
if regexp(expr, '\W') % wrap in quotes if invalid chars present in state name
    expr = ['"', expr, '"'];
end
end

function rOut = evaluateExternalFunctions(rIn, ids)
% Evaluate symbolic functions/pull in functions defined in path
%   Necessary for "power" and other MathML function translation
% Initialize symbolic variables
syms(ids{:});

% Evaluate the expressions to remove function calls
rOut = eval(rIn);

% Clear the symbolic variables
clear(ids{:})
end