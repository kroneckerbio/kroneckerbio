function SaveModel(m, filename)
%SaveModel writes the Kronecker mass action model to a Kronecker mass
%   action model file
%
%   SaveModel(m, filename)
%
%   Inputs
%   m: [ model struct scalar ]
%       The Kronecker model to be written.
%   filename: [ string ]
%       The path to the file that will be written. Existing files are
%       overwritten without prompt.

% (c) 2017 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if ~m.Ready
    m = FinalizeModel(m);
end

%% Open file
fid = fopen(filename, 'w+t');
keep_me_open = onCleanup(@()fclose(fid));

%% Save all compartments
% Write compartment header
fprintf(fid, '%% Compartments');
if ~isempty(m.Name)
    % Write model name
    fprintf(fid, ' %s', enquoted(m.Name));
end
fprintf(fid, '\n');

% Write comparments
for iv = 1:m.nv
    v = m.Compartments(iv);
    
    % Write name and dimension
    fprintf(fid, '%s %i', enquoted(v.Name), v.Dimension);
    
    % Write compartment size contributors
    n_expr = size(v.Size, 1);
    for i_expr = 1:n_expr
        % Write contributor expression
        if isempty(v.Size{i_expr, 1})
            % Empty name is constant expression
            fprintf(fid, ' %g', v.Size{i_expr, 2});
        elseif v.Size{i_expr, 2} == 1
            % Value of 1 is the default
            fprintf(fid, ' %s', enquoted(v.Size{i_expr, 1}));
        else
            fprintf(fid, ' %s=%g', enquoted(v.Size{i_expr, 1}), v.Size{i_expr, 2});
        end
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
    fprintf(fid, '%s %g\n', enquoted(k.Name), k.Value);
end

fprintf(fid, '\n');

%% Save all seeds
% Write parameter header
fprintf(fid, '%% Seeds\n');

% Write parameters
for is = 1:m.ns
    s = m.Seeds(is);
    
    % Write Name and Value
    fprintf(fid, '%s %g\n', enquoted(s.Name), s.Value);
end

fprintf(fid, '\n');

%% Save all inputs
% Initialize to dummy value to ensure new block runs on first iteration
current_compartment = nan;

% Write species
for iu = 1:numel(m.Inputs)
    u = m.Inputs(iu);
    
    % Start new block if compartment changes
    if ~isequaln(u.Compartment, current_compartment)
        if iu ~= 1; fprintf(fid, '\n'); end
        fprintf(fid, '%% Inputs %s\n', enquoted(u.Compartment));
        current_compartment = u.Compartment;
    end
    
    % Write Name
    fprintf(fid, '%s', enquoted(u.Name));
    
    % Write default value only if nonzero
    if u.DefaultValue ~= 0
        fprintf(fid, ' %g', u.DefaultValue);
    end
    
    % Write newline
    fprintf(fid, '\n');
end

fprintf(fid, '\n');

%% Save all states
% Initialize to dummy value to ensure new block runs on first iteration
current_compartment = nan;

% Write species
for ix = 1:numel(m.States)
    x = m.States(ix);
    
    % Start new block if compartment changes
    if ~isequaln(x.Compartment, current_compartment)
        if ix ~= 1; fprintf(fid, '\n'); end
        fprintf(fid, '%% States %s\n', enquoted(x.Compartment));
        current_compartment = x.Compartment;
    end
    
    % Write Name
    fprintf(fid, '%s', enquoted(x.Name));
    
    % Write initial condition seeds
    n_expr = size(x.InitialValue, 1);
    for i_expr = 1:n_expr
        % Write contributor expression
        if isempty(x.InitialValue{i_expr, 1})
            % Empty name is constant expression
            fprintf(fid, ' %g', x.InitialValue{i_expr, 2});
        elseif x.InitialValue{i_expr, 2} == 1
            % Value of 1 is the default
            fprintf(fid, ' %s', enquoted(x.InitialValue{i_expr, 1}));
        else
            fprintf(fid, ' %s=%g', enquoted(x.InitialValue{i_expr, 1}), x.InitialValue{i_expr, 2});
        end
    end
    
    % Write newline
    fprintf(fid, '\n');
end

fprintf(fid, '\n');

%% Save all reactions
% Initialize to dummy value to ensure new block runs on first iteration
current_compartment = nan;

% Write reactions
for ir = 1:m.nr
    r = m.Reactions(ir);
    
    % Start new block if compartment changes
    if ~isequaln(r.Compartment, current_compartment)
        if ir ~= 1; fprintf(fid, '\n'); end
        fprintf(fid, '%% Reactions %s\n', enquoted(r.Compartment));
        current_compartment = r.Compartment;
    end
    
    % Reactions with more than 2 products are handled specially
    large_scale = numel(r.Products) > 2;
    
    if large_scale
        fprintf(fid, ', ');
    end
    
    % Write reactants
    if numel(r.Reactants) == 0
        fprintf(fid, '0 0');
    elseif numel(r.Reactants) == 1
         fprintf(fid, '%s 0', enquoted(r.Reactants{1}));
    elseif numel(r.Reactants) == 2
        fprintf(fid, '%s %s', enquoted(r.Reactants{1}), enquoted(r.Reactants{2}));
    end
    
    % Write Products
    if numel(r.Products) == 0
        fprintf(fid, ' 0 0');
    elseif numel(r.Products) == 1
         fprintf(fid, ' %s 0', enquoted(r.Products{1}));
    elseif numel(r.Products) == 2
        fprintf(fid, ' %s %s', enquoted(r.Products{1}), enquoted(r.Products{2}));
    else
        fprintf(fid, [' ', strjoin(cellfun(@enquoted, r.Products, 'UniformOutput', false), ' ')]);
    end
    
    % Write parameter
    fprintf(fid, ' %s', enquoted(r.Parameter{1}));
    
    % Write parameter scale if present
    if r.Parameter{2} ~= 1
        fprintf(fid, '=%g', r.Parameter{2});
    end
    
    % Write name if present
    if ~isempty(r.Name)
        fprintf(fid, ' 0 %s', enquoted(r.Name));
    end
    
    % Write newline
    fprintf(fid, '\n');
end

fprintf(fid, '\n');

%% Save all outputs
% Write output header
fprintf(fid, '%% Outputs\n');

% Write outptus
for iy = 1:m.ny
    y = m.Outputs(iy);
    
    % Write name
    fprintf(fid, enquoted(y.Name));
    
    % Write output value contributors
    n_expr = size(y.Expression, 1);
    for i_expr = 1:n_expr
        % Write contributor expression
        if isempty(y.Expression{i_expr, 1})
            % Empty name is constant expression
            fprintf(fid, ' %g', y.Expression{i_expr, 2});
        elseif y.Expression{i_expr, 2} == 1
            % Value of 1 is the default
            fprintf(fid, ' %s', enquoted(y.Expression{i_expr, 1}));
        else
            fprintf(fid, ' %s=%g', enquoted(y.Expression{i_expr, 1}), y.Expression{i_expr, 2});
        end
    end
    
    % Write newline
    fprintf(fid, '\n');
end
end

function string = enquoted(string)
% Add quotes if string contains characters that need escaping
if ~isempty(regexp(string, '[%#\s,]', 'once'))
    string = ['"' string '"'];
end
end
