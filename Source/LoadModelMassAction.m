function m = LoadModelMassAction(files)
%LoadModelMassAction Load a model from the standard Kronecker model file
%   format
%
%   m = LoadModelMassAction(filenames)
%
%   Inputs:
%   files: [ cell vector of strings | string ]
%       The names of the files to be loaded
%
%   Outputs:
%   m: [ Model.MassActionKronecker scalar ]
%       The mass action Kronecker model constructed from the supplied files
%
%   The files are read in the order in which they appear in the files
%   vector.
%
%   If compartments, species, outputs, or parameters appear multiple times
%   throughout the files, the last instance by that name is the one that
%   will define that component's properties. All others by that same name
%   will be silently overwritten.

% (c) 2016 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Initialize the model
m = InitializeModelMassActionAmount('');

% Standardize files a cell vector of strings
if ischar(files)
    files = {files};
end
n_files = numel(files);

% Check that files exist
for i = 1:n_files
    assert(exist(files{i}, 'file') == 2, 'KroneckerBio:LoadModelMassAction:FileNotFound', 'File %s was not found', files{i})
end

%% Loop over each file
decimal_regex = '^[+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?$';

for i_file = 1:n_files
    % Open file
    fid = fopen(files{i_file});
    
    % Initialize
    mode = 0;
    line_number = 0;
    
    % Loop over all lines in file
    while true % break on -1
        % Get next line
        line = fgetl(fid);
        line_number = line_number + 1;
        
        % Break if end of file is encountered
        if line == -1
            break
        end
        
        % Strip leading and trailing whitespace
        line = strtrim(line);
        
        % Strip comments ('#' until end of line)
        %   Ignore '#' inside double-quotes '"'
        comment_inds = regexp(line, '#(?=([^"]*"[^"]*")*[^"]*$)');
        if ~isempty(comment_inds)
            line = line(1:comment_inds(1)-1);
        end
        
        % Skip if line is empty
        %   Comment lines starting with '#' are skipped in combination with above
        if isempty(line)
            continue
        end
        
        % Check if mode switching line is encountered
        % First character is a "%"
        is_switch_mode = ~isempty(regexp(line, '^%', 'once'));
        
        if is_switch_mode
            [match, payload] = regexp(line, '^%\s*\w+', 'match', 'split', 'once');
            payload = payload{2}; % Only want the end, not the starting empty string
            match = regexp(match, '\w+', 'match', 'once');
            
            if strcmpi(match, 'compartments')
                % Switch mode
                mode = 1;
                
                % Extract model name
                model_name = regexp(payload, '".*"|[^\s]*', 'match', 'once');
                model_name = regexp(model_name, '[^"]*', 'match', 'once'); % Strip quotes
                m = RenameModel(m, model_name);
            elseif strcmpi(match, 'inputs')
                % Switch mode
                mode = 2;
                
                % Extract default compartment from payload
                compartment = regexp(payload, '"[^"]*"|[^."\s]+', 'match');
                assert(numel(compartment) == 1, 'KroneckerBio:LoadModelMassAction:SpeciesCompartments', 'Line %i in %s specified the beginning of an Input section but exactly one compartment was not provided: %s', line_number, files{i_file}, line)
                
                % Standardize compartment as a string
                if isempty(compartment)
                    compartment = '';
                else
                    compartment = compartment{1};
                end
            elseif strcmpi(match, 'seeds')
                % Switch mode
                mode = 3;
            elseif strcmpi(match, 'states')
                % Switch mode
                mode = 4;
                
                % Extract default compartment from payload
                compartment = regexp(payload, '"[^"]*"|[^."\s]+', 'match');
                assert(numel(compartment) == 1, 'KroneckerBio:LoadModelMassAction:SpeciesCompartments', 'Line %i in %s specified the beginning of an States section but exactly one compartment was not provided: %s', line_number, files{i_file}, line)
                
                % Standardize compartment as a string
                if isempty(compartment)
                    compartment = '';
                else
                    compartment = compartment{1};
                end
                
            elseif strcmpi(match, 'outputs')
                % Switch mode
                mode = 5;
            elseif strcmpi(match, 'parameters')
                % Switch mode
                mode = 6;
            elseif strcmpi(match, 'reactions')
                % Switch mode
                mode = 7;
            else
                error('KroneckerBio:LoadModelMassAction:UnrecognizedSectionHeader', 'Line %i in %s had an unrecognized section header: %s', line_number, files{i_file}, line)
            end
        else
            % Process entries
            
            % Tokenize
            tokens = vec(regexp(line, '"[^"]*"|[^"\s]+','match'));
            
            % Strip quotes
            for i_tok = 1:numel(tokens)
                tokens{i_tok} = regexp(tokens{i_tok}, '[^"]*', 'match', 'once');
            end
            
            if mode == 1
                %% Compartments
                assert(numel(tokens) >= 3, 'KroneckerBio:LoadModelMassAction:InsufficientCompartmentInformation', 'Line %i in %s does not provide the name, dimension, and size required for a compartment: %s', line_number, files{i_file}, line)
                
                % Extract name
                assert(~isempty(regexp(tokens{1}, '^[^.]+$', 'once')), 'KroneckerBio:LoadModelMassAction:InvalidCompartmentName', 'Line %i in %s has an invalid compartment name: %s', line_number, files{i_file}, line)
                name = tokens{1};
                
                % Extract dimension
                assert(~isempty(regexp(tokens{2}, '^0|1|2|3$', 'once')), 'KroneckerBio:LoadModelMassAction:InvalidCompartmentDimension', 'Line %i in %s has an invalid compartment dimension: %s', line_number, files{i_file}, line)
                dim = str2double(tokens{2});
                
                % Extract value/size
                tokens = tokens(3:end);
                exprs = extract_values(tokens);
                
                % Add compartment
                m = AddCompartment(m, name, dim, exprs);
            elseif mode == 2
                %% Inputs
                % Extract name
                assert(~isempty(regexp(tokens{1}, '^([^.]+\.)?[^.]+$', 'once')), 'KroneckerBio:LoadModelMassAction:InvalidSpeciesName', 'Line %i in %s has an invalid species name: %s', line_number, files{i_file}, line)
                name = tokens{1};
                
                % Second token is value
                if numel(tokens) >= 2
                    assert(~isempty(regexp(tokens{2}, decimal_regex, 'once')), 'KroneckerBio:LoadModelMassAction:InvalidInputValue', 'Line %i in %s has an input with a value that is not a number: %s', line_number, files{i_file}, line)
                    value = str2double(tokens{2});
                else
                    value = 0;
                end
                
                % Add Input
                m = AddInput(m, name, compartment, value);
            elseif mode == 3
                %% Seeds
                % Extract name
                assert(~isempty(regexp(tokens{1}, '^[^.]+$', 'once')), 'KroneckerBio:LoadModelMassAction:InvalidParameterName', 'Line %i in %s has an invalid seed name: %s', line_number, files{i_file}, line)
                name = tokens{1};
                
                % Second token is the number
                value = eval(tokens{2});
                
                % Add seed
                m = AddSeed(m, name, value);
            elseif mode == 4
                %% States
                % Extract name
                assert(~isempty(regexp(tokens{1}, '^([^.]+\.)?[^.]+$', 'once')), 'KroneckerBio:LoadModelMassAction:InvalidName', 'Line %i in %s has an invalid species name: %s', line_number, files{i_file}, line)
                name = tokens{1};
                
                % Process seeds
                tokens = tokens(2:end);
                seeds = extract_values(tokens);
                
                % Add State
                m = AddState(m, name, compartment, seeds);
            elseif mode == 5
                %% Outputs
                % Extract name
                assert(~isempty(regexp(tokens{1}, '^[^.]+$', 'once')), 'KroneckerBio:LoadModelMassAction:InvalidName', 'Line %i in %s has an invalid output name: %s', line_number, files{i_file}, line)
                name = tokens{1};
                
                % Process outputs
                tokens = tokens(2:end);
                exprs = extract_values(tokens);
                
                % Add output
                m = AddOutput(m, name, exprs);
            elseif mode == 6
                %% Parameters
                % Extract parameter name
                assert(~isempty(regexp(tokens{1}, '^[^.]+$', 'once')), 'KroneckerBio:LoadModelMassAction:InvalidName', 'Line %i in %s has an invalid parameter name: %s', line_number, files{i_file}, line)
                name = tokens{1};
                
                % Second token is value
                assert(numel(tokens) >= 2, 'KroneckerBio:LoadModelMassAction:MissingParameterValue', 'Line %i in %s is missing its parameter value: %s', line_number, files{i_file}, line)
                value = eval(tokens{2});
                
                % Add parameter
                m = AddParameter(m, name, value);
            elseif mode == 7
                %% Reactions
                if ~strcmp(tokens{1}, ',')
                    % Add empty reverse parameter if not specified
                    if numel(tokens) < 6
                        tokens  = [tokens; {''}];
                    end
                    
                    parameters = cell(2,2);
                    for i = 1:size(parameters,1)
                        % Split between the equal sign
                        elements = regexp(tokens{i+4}, '=', 'split');
                        assert(numel(elements) <= 2, 'KroneckerBio:LoadModelMassAction:InvalidExpression', 'Line %i in %s has more than one equal sign in an expression: %s', line_number, files{i_file}, line)
                        
                        if strcmp(elements{1}, '0')
                            % It is an empty reaction
                            parameters(i,:) = {'', 0};
                        elseif numel(elements) == 2
                            % It is a parameter and a modifier
                            parameters(i,:) = {elements{1}, eval(elements{2})};
                        else
                            % It is a single parameter
                            parameters(i,:) = {elements{1}, 1};
                        end
                    end
                    
                    reactants = tokens(1:2);
                    reactants = reactants(~strcmp(reactants, '0'));
                    products = tokens(3:4);
                    products = products(~strcmp(products, '0'));
                    
                    names = tokens(7:end);
                    
                    % This is a normal line
                    m = AddReaction(m, names, reactants, products, parameters(1,:), parameters(2,:));
                else
                    % Split between the equal sign
                    parameter = regexp(tokens{end}, '=', 'split');
                    assert(numel(parameter) <= 2, 'KroneckerBio:LoadModelMassAction:InvalidExpression', 'Line %i in %s has more than one equal sign in an expression: %s', line_number, files{i_file}, line)
                    
                    if numel(parameter) == 2
                        % It is a parameter and a modifier
                        parameter = {parameter{1}, eval(parameter{2})};
                    else
                        % It is a single parameter
                        parameter = {parameter{1}, 1};
                    end
                    
                    reactants = tokens(2:3);
                    reactants = reactants(~strcmp(reactants, '0'));
                    products = tokens(4:end-1);
                    products = products(~strcmp(products, '0'));
                    
                    % This is the special syntax for a multi-product line
                    m = AddReaction(m, [], reactants, products, parameter);
                end
            else
                error('KroneckerBio:LoadModelMassAction:SectionNotSpecified', 'Line %i in %s occurs before any section header: %s', line_number, files{i_file}, line)
            end
        end
    end
    
    % Close the model file
    fclose(fid);
end

%% Finalize Model
m = FinalizeModel(m);

%% Helper function
    function exprs = extract_values(tokens)
        % Extract components of linear combination in cell matrix form
        %
        % Inputs
        %   tokens [ 1 x nTokens cell vector of strings ]
        %       Input tokens to extract
        %
        % Output:
        %   exprs [ nExpr x 2 cell matrix of (string,double) rows ]
        %       Output cell matrix used in Add* functions
        n_expr = numel(tokens);
        exprs = cell(n_expr,2);
        for i_expr = 1:n_expr
            % Split between the equal sign
            elements = regexp(tokens{i_expr}, '=', 'split');
            assert(numel(elements) <= 2, 'KroneckerBio:LoadModelMassAction:InvalidExpression', 'Line %i in %s has more than one equal sign in an expression: %s', line_number, files{i_file}, line)
            
            if numel(elements) == 2
                % It is a expression and a modifier
                exprs(i_expr,:) = {elements{1}, str2double(elements{2})};
            elseif numel(elements) == 1 && ~isempty(regexp(elements{1}, '^\d', 'once'))
                % It is a constitutive value
                exprs(i_expr,:) = {'',  str2double(elements{1})};
            else
                % It is a single expression
                exprs(i_expr,:) = {elements{1}, 1};
            end
        end
    end
end
