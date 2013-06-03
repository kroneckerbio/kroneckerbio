function m = symbolic2MassAction(SymModel, opts)

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1, 0);

%% Extract symbolic components
name       = SymModel.Name;

vSyms     = SymModel.vSyms;
vNames    = SymModel.vNames;
dv        = SymModel.dv;
vValues   = SymModel.vValues;

xuSyms    = SymModel.xuSyms;
xuNames   = SymModel.xuNames;
vxuNames  = SymModel.vxuNames;
isu       = SymModel.isu;
xu0       = SymModel.xu0;

kSyms     = SymModel.kSyms;
kNames    = SymModel.kNames;
k         = SymModel.k;

rNames    = SymModel.rNames;
rSyms     = SymModel.r;
S         = SymModel.S;

%% Basic constants
% Basic constants
nv = numel(vNames);
nxu = numel(xuNames);
nk = numel(kNames);
nEqu = numel(rSyms);

xuNamesFull = strcat(vxuNames, '.', xuNames);

%% Sanity check
assert(all(vec(fix(S) == S)) && all(vec(abs(S) <= 2)), 'KroneckerBio:symbolic2MassAction:InvalidStoichiometry', 'Invalid stoichiometry matrix does not consist entirely of -2, -1, 0, 1, and 2')

%% Initialize model
% Basic structure
m = InitializeModel(name);

%% Loop over compartments
for iv = 1:nv
    m = AddCompartment(m, vNames{iv}, dv(iv), '', vValues(iv));
end

%% Loop over species
for ixu = 1:nxu
    m = AddSpecies(m, xuNames{ixu}, vxuNames{ixu}, xu0(ixu), '', isu(ixu));
end

%% Loop over parameters
for ik = 1:nk
    m = AddParameter(m, kNames{ik}, k(ik));
end

%% Loop over reaction rate equations
if verbose; fprintf('Parsing reaction symbolics to kronecker reactions...'); end

% Note second-order parameters that have the volume baked in
bakedk = zeros(nk,1);

for iEqu = 1:nEqu
    % Extract all variables used in this reaction rate
    allvar = symvar(rSyms(iEqu));
    
    % Variables that are rate constants; index in kSyms
    kPosEqu = lookup(allvar, kSyms);
    kPosEqu = kPosEqu(kPosEqu ~= 0);
    
    % Number of microscopic reactions in this rate
    nMic = numel(kPosEqu);
    
    % Count of products and reactants in stoichiometry matrix
    % Double up on species that get twice as much
    reac = [find(S(:,iEqu) < 0); find(S(:,iEqu) == -2)];
    prod = [find(S(:,iEqu) > 0); find(S(:,iEqu) == 2)];
    nreac = numel(reac);
    nprod = numel(prod);
    assert(nreac <= 2 && nprod <= 2, 'KroneckerBio:symbolic2MassAction:InvalidStoichiometry', 'Reaction rate %i has an unacceptable stoichiometry vector', iEqu)
    
    % Loop over each microscopic reaction and add it
    for iMic = 1:nMic
        % Differentiate wrt the rate constant to return only the species
        derMic = diff(rSyms(iEqu), kSyms(kPosEqu(iMic)));
        
        % Extract those species involved
        allxu = symvar(derMic);
        
        % Ensure second order or less with species
        xuPosMic = lookup(allxu, xuSyms);
        xuPosMic = xuPosMic(xuPosMic ~= 0);
        nxuMic = numel(xuPosMic);
        assert(all(xuPosMic ~= 0) && nxuMic <= 2, 'KroneckerBio:symbolic2MassAction:InvalidReaction', 'Reaction rate %i does not represent a quadratic ODE', iEqu)
        
        % Switch one zeroth, first, or second order
        switch(nxuMic)
            case(0)
                % Zeroth-order reaction
                assert(derMic == -1 || derMic == 1, 'KroneckerBio:symbolic2MassAction:ZerothOrderCheck', 'Reaction rate %i has a zeroth-order term that failed a sanity check', iEqu)
                
                if derMic == 1
                    % Add zeroth order synthesis forward
                    assert(nreac == 0, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a zeroth-order term that does not match the stoichiometry', iEqu)
                    
                    % Switch on number of products
                    switch(nprod)
                        case(0)
                            % Pointless reaction
                            m = AddReaction(m, rNames{iEqu}, '', '', '', '', '', kNames{kPosEqu(iMic)});
                        case(1)
                            m = AddReaction(m, rNames{iEqu}, '', '', '', xuNamesFull{prod}, '', kNames{kPosEqu(iMic)});
                        case(2)
                            m = AddReaction(m, rNames{iEqu}, '', '', '', xuNamesFull{prod(1)}, xuNamesFull{prod(2)}, kNames{kPosEqu(iMic)});
                    end
                elseif derMic == -1
                    % Add zeroth order synthesis reverse
                    assert(nprod == 0, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a zeroth-order term that does not match the stoichiometry', iEqu)
                    
                    % Switch on number of "reactants"
                    switch(nreac)
                        case(0)
                            % Pointless reaction
                            m = AddReaction(m, rNames{iEqu}, '', '', '', '', '', kNames{kPosEqu(iMic)});
                        case(1)
                            m = AddReaction(m, rNames{iEqu}, '', '', '', xuNamesFull{reac}, '', kNames{kPosEqu(iMic)});
                        case(2)
                            m = AddReaction(m, rNames{iEqu}, '', '', '', xuNamesFull{reac(1)}, xuNamesFull{reac(2)}, kNames{kPosEqu(iMic)});
                    end
                end
            case(1)
                % Derivative check for first-order reaction
                derMic = diff(derMic, xuSyms(xuPosMic));
                
                % Switch on actual first order or self second-order
                if derMic == -1 || derMic == 1
                    % True first-order reaction
                    
                    if derMic == 1
                        % Add first order forward
                        assert(nreac == 1, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a single species term that does not match the stoichiometry', iEqu)
                        
                        % Switch on number of products
                        switch(nprod)
                            case(0)
                                m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic}, '', '', '', kNames{kPosEqu(iMic)});
                            case(1)
                                m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic}, '', xuNamesFull{prod}, '', kNames{kPosEqu(iMic)});
                            case(2)
                                m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic}, '', xuNamesFull{prod(1)}, xuNamesFull{prod(2)}, kNames{kPosEqu(iMic)});
                        end
                    elseif derMic == -1
                        % Add first order reverse
                        assert(nprod == 1, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a single species term that does not match the stoichiometry', iEqu)
                        
                        % Switch on number of "reactants"
                        switch(nreac)
                            case(0)
                                m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic}, '', '', '', kNames{kPosEqu(iMic)});
                            case(1)
                                m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic}, '', xuNamesFull{reac}, '', kNames{kPosEqu(iMic)});
                            case(2)
                                m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic}, '', xuNamesFull{reac(1)}, xuNamesFull{reac(2)}, kNames{kPosEqu(iMic)});
                        end
                    end
                else
                    % Species reacts with itself
                    derMic = diff(derMic, xuSyms(xuPosMic));
                    allv = symvar(derMic);
                    
                    % Process compartment volume
                    vPosMic = lookup(allv, vSyms);
                    vPosMic = vPosMic(vPosMic ~= 0);
                    if numel(vPosMic) == 0
                        % Compartment is baked in
                        vPosMic = find(strcmp(vxuNames{xuPosMic}, vNames), 1);
                        
                        % If we are going to change the parameter it must
                        % not be used by multiple compartments
                        assert(bakedk(kPosEqu(iMic)) == 0 || bakedk(kPosEqu(iMic)) == vPosMic || vValues(bakedk(kPosEqu(iMic))) == vValues(vPosMic), 'KroneckerBio:symbolic2MassAction:BakedReactionReuse', 'Reaction rate %i uses a second-order term that has the compartment volume baked into the kinetic parameter, but this parameter is shared between compartments with different volumes')
                        
                        % Claim this parameter for this compartment
                        if kPosEqu(iMic) == 0
                            bakedk(kPosEqu(iMic)) = vPosMic;
                        
                            % Replace the old parameter value
                            m = AddParameter(m, kNames{kPosEqu(iMic)}, k(kPosEqu(iMic)) * vValues(vPosMic));
                        end
                    elseif numel(vPosMic) == 1
                        % Compartment is written in the equation
                        derMic = derMic * vSyms(vPosMic);
                        assert(derMic == 2 || derMic == -2, 'KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', iEqu)
                        assert(vxuSyms(xuPosMic) == vSyms(vPosMic), 'KroneckerBio:symbolic2MassAction:CompartmentMembership', 'Reaction rate %i uses a compartment that neither reactant is a member of', iEqu)
                    else
                        error('KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', iEqu)
                    end
                    
                    assert(derMic == 2 || derMic == -2, 'KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', iEqu)
                    
                    if derMic == 2
                        % Add second order forward
                        assert(nreac == 2, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a single species term that does not match the stoichiometry', iEqu)
                        
                        % Switch on number of products
                        switch(nprod)
                            case(0)
                                m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic}, xuNamesFull{xuPosMic}, '', '', kNames{kPosEqu(iMic)});
                            case(1)
                                m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic}, xuNamesFull{xuPosMic}, xuNamesFull{prod}, '', kNames{kPosEqu(iMic)});
                            case(2)
                                m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic}, xuNamesFull{xuPosMic}, xuNamesFull{prod(1)}, xuNamesFull{prod(2)}, kNames{kPosEqu(iMic)});
                        end
                    elseif derMic == -2
                        % Add second order reverse
                        assert(nprod == 2, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a single species term that does not match the stoichiometry', iEqu)
                        
                        % Switch on number of "reactants"
                        switch(nreac)
                            case(0)
                                m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic}, xuNamesFull{xuPosMic}, '', '', kNames{kPosEqu(iMic)});
                            case(1)
                                m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic}, xuNamesFull{xuPosMic}, xuNamesFull{reac}, '', kNames{kPosEqu(iMic)});
                            case(2)
                                m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic}, xuNamesFull{xuPosMic}, xuNamesFull{reac(1)}, xuNamesFull{reac(2)}, kNames{kPosEqu(iMic)});
                        end
                    end
                end
            case(2)
                % Derivative check for second-order reaction
                derMic = diff(derMic, xuSyms(xuPosMic(1)));
                derMic = diff(derMic, xuSyms(xuPosMic(2)));
                allv = symvar(derMic);
                
                % Process compartment volume
                vPosMic = lookup(allv, vSyms);
                vPosMic = vPosMic(vPosMic ~= 0);
                if numel(vPosMic) == 0
                    % Compartment is baked in
                    vPosMic = lookup(vxuNames(xuPosMic), vNames);
                    
                    % If we are going to change the parameter it must
                    % not be used by multiple compartments
                    assert(bakedk(kPosEqu(iMic)) == 0 || any(bakedk(kPosEqu(iMic)) == vPosMic) || any(vValues(bakedk(kPosEqu(iMic))) == vValues(vPosMic)), 'KroneckerBio:symbolic2MassAction:BakedReactionReuse', 'Reaction rate %i uses a second-order term that has the compartment volume baked into the kinetic parameter, but this parameter is shared between compartments with different volumes')
                    
                    % Claim this parameter for this compartment
                    if kPosEqu(iMic) == 0
                        % Use the freest compartment
                        [unused ind] = max(d(vPosMic));
                        vPosMic = vPosMic(ind);
                        bakedk(kPosEqu(iMic)) = vPosMic;
                        
                        % Replace the old parameter value
                        m = AddParameter(m, kNames{kPosEqu(iMic)}, k(kPosEqu(iMic)) * vValues(vPosMic));
                    end
                elseif numel(vPosMic) == 1
                    % Compartment is written in the equation
                    derMic = derMic * vSyms(vPosMic);
                    assert(derMic == 2 || derMic == -2, 'KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', iEqu)
                    assert(vxuSyms(xuPosMic) == vSyms(vPosMic), 'KroneckerBio:symbolic2MassAction:CompartmentMembership', 'Reaction rate %i uses a compartment that neither reactant is a member of', iEqu)
                else
                    error('KroneckerBio:symbolic2MassAction:FirstOrderCheck', 'Reaction rate %i has a single species term that failed a sanity check', iEqu)
                end
                
                assert(derMic == 1 || derMic == -1, 'KroneckerBio:symbolic2MassAction:SecondOrderCheck', 'Reaction rate %i has a double species term that failed a sanity check', iEqu)
                
                if derMic == 1
                    % Add second order forward
                    assert(nreac == 2, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a double species term that does not match the stoichiometry', iEqu)
                    
                    % Switch on number of products
                    switch(nprod)
                        case(0)
                            m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic(1)}, xuNamesFull{xuPosMic(2)}, '', '', kNames{kPosEqu(iMic)});
                        case(1)
                            m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic(1)}, xuNamesFull{xuPosMic(2)}, xuNamesFull{prod}, '', kNames{kPosEqu(iMic)});
                        case(2)
                            m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic(1)}, xuNamesFull{xuPosMic(2)}, xuNamesFull{prod(1)}, xuNamesFull{prod(2)}, kNames{kPosEqu(iMic)});
                    end
                elseif derMic == -1
                    % Add second order reverse
                    assert(nprod == 2, 'KroneckerBio:symbolic2MassAction:StoichiometryMismatch', 'Reaction rate %i has a double species term that does not match the stoichiometry', iEqu)
                    
                    % Switch on number of "reactants"
                    switch(nreac)
                        case(0)
                            m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic(1)}, xuNamesFull{xuPosMic(2)}, '', '', kNames{kPosEqu(iMic)});
                        case(1)
                            m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic(1)}, xuNamesFull{xuPosMic(2)}, xuNamesFull{reac}, '', kNames{kPosEqu(iMic)});
                        case(2)
                            m = AddReaction(m, rNames{iEqu}, '', xuNamesFull{xuPosMic(1)}, xuNamesFull{xuPosMic(2)}, xuNamesFull{reac(1)}, xuNamesFull{reac(2)}, kNames{kPosEqu(iMic)});
                    end
                end
                
        end
    end
end
if verbose; fprintf('done.\n'); end

%% Finalize model
m = FinalizeModel(m);
