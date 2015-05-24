function SaveModel(filename, m)
%SaveModel writes the Kronecker mass action model to a Kronecker mass
%   action model file
%
%   SaveModel(filename, m)
%
%   Inputs
%   filename: [ string ]
%       The path to the file that will be written. Existing files are
%       overwritten without prompt.
%   m: [ model struct scalar ]
%       The Kronecker model to be written.

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Open file
fid = fopen(filename, 'w+t');

%% Save all compartments
% Write compartment header
fprintf(fid, '%% Compartments');
if ~isempty(m.Name)
    % Write model name
    if ~isempty(regexp(m.Name, '[\s,]', 'once'))
        % Enclose name in quotes if it contains operators
        m.Name = ['"' m.Name '"'];
    end
    fprintf(fid, ' %s', m.Name);
end
fprintf(fid, '\n');

% Write comparments
for iv = 1:m.nv
    v = m.Compartments(iv);
    
    % Write name and dimension
    if ~isempty(regexp(v.Name, '^[%#]|[\s,]', 'once'))
        % Enclose name in quotes if it contains operators
        v.Name = ['"' v.Name '"'];
    end
    fprintf(fid, '%s %i', v.Name, v.Dimension);
    
    % Write compartment volume regular expressions
    nExpr = numel(v.Expressions);
    for iExpr = 1:nExpr
        fprintf(fid, ' ');
        
        % Write expressions if one exists
        if ~isempty(v.Expressions{iExpr})
            fprintf(fid, '%s=', v.Expressions{iExpr});
        end
        
        % Write value
        fprintf(fid, '%g', v.Values(iExpr));
        
        % Write comma if there are more expressions,
        if iExpr == nExpr
            break
        end
        fprintf(fid, ',');
    end
    
    % Write newline
    fprintf(fid, '\n');
end

fprintf(fid, '\n');

%% Save all parameters
% Write parameter header
fprintf(fid, '%% Parameters\n');

% Write parameters
for ik = 1:m.nk
    k = m.Parameters(ik);
    
    % Write Name and Value
    if ~isempty(regexp(k.Name, '^[%#]|[\s,]', 'once'))
        % Enclose name in quotes if it contains operators
        k.Name = ['"' k.Name '"'];
    end
    fprintf(fid, '%s %g\n', k.Name, k.Value);
end

fprintf(fid, '\n');

%% Save all seeds
% Write parameter header
fprintf(fid, '%% Seeds\n');

% Write parameters
for ik = 1:m.nk
    s = m.Seeds(ik);
    
    % Write Name and Value
    if ~isempty(regexp(s.Name, '^[%#]|[\s,]', 'once'))
        % Enclose name in quotes if it contains operators
        s.Name = ['"' s.Name '"'];
    end
    fprintf(fid, '%s %g\n', k.Name, k.Value);
end

fprintf(fid, '\n');

%% Save all inputs
% Write species header
fprintf(fid, '%% Species\n');

% Write species
for iu = 1:numel(m.Inputs)
    u = m.Species(iu);
    
    % Write Compartment.Name
    name = [u.Compartment '.' u.Name];
    if ~isempty(regexp(name, '^[%#]|[\s,]', 'once'))
        % Enclose name in quotes if it contains operators
        v.Name = ['"' name '"'];
    end
    fprintf(fid, '%s', name);
    
    % Write function
    if isempty(regexp(func2str(u.Value.Function), '@\(t,q\)repmat\([-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?,1,numel\(t\)\)', 'once'))
        % Input is not a constant
        fprintf(fid, ' %s', func2str(u.Function));
    else
        % The function was originally a constant; extract it
        match = regexp(func2str(u.Function), '[-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?', 'match', 'once');
        fprintf(fid, ' %s', match);
    end
    
    % If there are parameters, write them
    nPara = numel(u.Parameters);
    for iPara = 1:nPara
        % Write parameter
        fprintf(fid, ' %g', u.Parameters(iPara));
        
        % Write comma if there are more parameters
        if iExpr == nExpr
            break
        end
        fprintf(fid, ',');
    end
    
    % Write newline
    fprintf(fid, '\n');
end

fprintf(fid, '\n');

%% Save all states
% Write species header
fprintf(fid, '%% Species\n');

% Write species
for ix = 1:numel(m.States)
    x = m.States(ix);
    
    % Write Compartment.Name
    name = [x.Compartment '.' x.Name];
    if ~isempty(regexp(name, '^[%#]|[\s,]', 'once'))
        % Enclose name in quotes if it contains operators
        v.Name = ['"' name '"'];
    end
    fprintf(fid, '%s', name);
    
    % Write value only if nonzero
    if x.InitialValue ~= 0
        fprintf(fid, ' %g', x.InitialValue);
    end
    
    % Write newline
    fprintf(fid, '\n');
end

fprintf(fid, '\n');

%% Save all reactions
% Write reaction header
fprintf(fid, '%% Reactions\n');

% Write reactions
for ir = 1:m.nr
    r = m.Reactions(ir);
    
    % Write Reactants
    if ~isempty(r.Reactants{1})
        if ~isempty(regexp(r.Reactants{1}, '^[%#]|[\s,]', 'once'))
            % Enclose name in quotes if it contains operators
            r.Reactants{1} = ['"' r.Reactants{1} '"'];
        end
        fprintf(fid, [r.Reactants{1} ' ']);
    else
        fprintf(fid, '0 ');
    end
    if ~isempty(r.Reactants{2})
        if ~isempty(regexp(r.Reactants{2}, '^[%#]|[\s,]', 'once'))
            % Enclose name in quotes if it contains operators
            r.Reactants{2} = ['"' r.Reactants{2} '"'];
        end
        fprintf(fid, [r.Reactants{2} ' ']);
    else
        fprintf(fid, '0 ');
    end
    
    % Write Products
    if ~isempty(r.Products{1})
        if ~isempty(regexp(r.Products{1}, '^[%#]|[\s,]', 'once'))
            % Enclose name in quotes if it contains operators
            r.Products{1} = ['"' r.Products{1} '"'];
        end
        fprintf(fid, [r.Products{1} ' ']);
    else
        fprintf(fid, '0 ');
    end
    if ~isempty(r.Products{2})
        if ~isempty(regexp(r.Products{2}, '^[%#]|[\s,]', 'once'))
            % Enclose name in quotes if it contains operators
            r.Products{2} = ['"' r.Products{2} '"'];
        end
        fprintf(fid, [r.Products{2} ' ']);
    else
        fprintf(fid, '0 ');
    end
    
    % Write Parameter
    if ~isempty(regexp(r.Parameter, '^[%#]|[\s,]', 'once'))
        % Enclose name in quotes if it contains operators
        r.Parameter = ['"' r.Parameter '"'];
    end
    fprintf(fid, r.Parameter);
    
    % Write Name
    if ~isempty(r.Name)
        if ~isempty(regexp(r.Name, '[\s,]', 'once'))
            % Enclose name in quotes if it contains operators
            r.Name = ['"' r.Name '"'];
        end
        fprintf(fid, [' 0 ' r.Name]);
    end
    
    % Write newline
    fprintf(fid, '\n');
end

%% Save all outputs
% Write output header
fprintf(fid, '%% Outputs\n');

% Write outptus
for iy = 1:m.ny
    y = m.Outputs(iy);
    
    % Write name
    if ~isempty(regexp(y.Name, '^[%#]|[\s,]', 'once'))
        % Enclose name in quotes if it contains operators
        y.Name = ['"' y.Name '"'];
    end
    fprintf(fid, y.Name);
    
    % Write output value regular expressions
    nExpr = numel(y.Expressions);
    for iExpr = 1:nExpr
        fprintf(fid, ' ');
        
        % Write expressions if one exists
        if ~isempty(y.Expressions(iExpr))
            fprintf(fid, '%s=', y.Expressions{iExpr});
        end
        
        % Write value
        fprintf(fid, '%g', y.Values(iExpr));
        
        % Write comma if there are more expressions
        if iExpr == nExpr
            break
        end
        fprintf(fid, ',');
    end
    
    % Write newline
    fprintf(fid, '\n');
end

fprintf(fid, '\n');

%% Close file
fclose(fid);
