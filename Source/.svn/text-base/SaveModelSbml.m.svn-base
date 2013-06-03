function mSimbio = SaveModelSbml(m, file)

%% Convert Kronecker model to SimBiology
% New model
mSimbio = sbiomodel(m.Name);

%% Add Compartments
vIsConstant = true(m.nv,1);
for iv = 1:m.nv
    v = m.Compartments(iv);
    
    % Determine if this comparment is constant
    vConstant = 0;
    for iExpr = 1:numel(v.Expressions)
        if ~isempty(v.Expressions{iExpr})
            vIsConstant(iv) = false;
            vConstant = v.Values(iExpr);
            break
        end
    end
    
    % Add compartment
    mSimbio.addcompartment(v.Name, vConstant, 'ConstantCapacity', vIsConstant(iv));
end

%% Add Compartment Rules
% Combine B1 and B2 for easier species indexing
B = zeros(m.nv,numel(m.Species));
B(:,~m.isu) = m.B1;
B(:,m.isu)  = m.B2;

for iv = 1:find(~vIsConstant')
    v = m.Compartments(iv);
    
    % Initialize rule
    rule = '';
    nEle = 0;
    
    % Add species elements
    ele = find(B(iv,:));
    for iEle = 1:numel(ele)
        % Increment the number of elements
        nEle = nEle + 1;
        if nEle > 1
            % Add a plus sign before all but the first element
            rule = [rule '+'];
        end
        % Add the element
        if B(iv,ele(iEle)) == 1
            % There is no need to multiply by 1
            rule = sprintf('%s%s.%s', rule, m.Species(ele(iEle)).Compartment, m.Species(ele(iEle)).Name);
        else
            % Multiply by the value
            rule = sprintf('%s%g*%s.%s', rule, B(iv,ele(iEle)), m.Species(ele(iEle)).Compartment, m.Species(ele(iEle)).Name);
        end
    end
    
    % Add rule
    if m.b(iv)
        % There is a constant to be added
        mSimbio.addrule(sprintf('%s = %g+%s', v.Name, full(m.b(iv)), rule), 'repeatedAssignment');
    else
        % There is no constant
        mSimbio.addrule(sprintf('%s = %s', v.Name, rule), 'repeatedAssignment');
    end
end

%% Add Species
xuIsVarInput = false(numel(m.Species),1);
for ixu = 1:numel(m.Species)
    xu = m.Species(ixu);
    
    % Find compartment object
    vObj = sbioselect(mSimbio, 'Type', 'compartment', 'Name', xu.Compartment);
    assert(numel(vObj) == 1, 'KroneckerBio:SaveModelSbml:MissingCompartment', 'SimBiology could not find compartment %s for species %s', xu.Compartment, xu.Name)
    
    if ~xu.IsInput
        % Add state species
        value = xu.Value;
    else
        % Determine if this input is variable input
        if isempty(regexp(func2str(xu.Value.Function), '@\(t,q\)repmat\([-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?,1,numel\(t\)\)', 'once'))
            % Input is variable
            xuIsVarInput(ixu) = true;
            value = 0;
        else
            % The function was originally a constant; extract it
            match = regexp(func2str(xu.Value.Function), '[-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?', 'match', 'once');
            value = str2double(match);
        end
    end
    
    % Add input species
    vObj.addspecies(xu.Name, value, 'BoundaryCondition', xu.IsInput);
end

%% Add Input Rules
for ixu = find(xuIsVarInput')
    u = m.Species(ixu);
    
    mSimbio.addrule(sprintf('%s = %s', u.Name, u.Value.Function), 'repeatedAssignment');
end

%% Add Outputs
% Combine C1 and C2 for easier species indexing
C = zeros(m.ny,numel(m.Species));
isu = cat(1, m.Species.IsInput);
C(:,~isu) = m.C1;
C(:,isu)  = m.C2;

for iy = 1:m.ny
    y = m.Outputs(iy);
    
    % Initialize rule
    rule = '';
    nEle = 0;
    
    % Add species elements
    ele = find(C(iy,:));
    for iEle = 1:numel(ele)
        % Increment the number of elements
        nEle = nEle + 1;
        if nEle > 1
            % Add a plus sign before all but the first element
            rule = [rule '+'];
        end
        % Add the element
        if C(iy,ele(iEle)) == 1
            % There is no need to multiply by 1
            rule = sprintf('%s%s.%s', rule, m.Species(ele(iEle)).Compartment, m.Species(ele(iEle)).Name);
        else
            % Multiply by the value
            rule = sprintf('%s%g*%s.%s', rule, C(iy,ele(iEle)), m.Species(ele(iEle)).Compartment, m.Species(ele(iEle)).Name);
        end
    end
    
    % Add rule
    if m.c(iy)
        % There is a constant to be added
        mSimbio.addrule(sprintf('%s = %g+%s', y.Name, full(m.c(iy)), rule), 'repeatedAssignment');
    else
        % There is no constant
        mSimbio.addrule(sprintf('%s = %s', y.Name, rule), 'repeatedAssignment');
    end
end

%% Add Parameters
for ik = 1:m.nk
    k = m.Parameters(ik);
    mSimbio.addparameter(k.Name, k.Value);
end

%% Add Reactions
for ir = 1:n.nr
    r = m.Reactions(ir);
    
    mSimbio.addreaction(r.Reactants, r.Products);
end

%% Save SimBiology model to file
sbmlexport(mSimbio, file);