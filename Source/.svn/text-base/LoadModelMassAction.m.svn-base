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

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Initialize the model
m = InitializeModel('');

% Standardize files a cell vector of strings
if ischar(files)
    files = {files};
end
nFiles = numel(files);

% Check that files exist
for i = 1:nFiles
    assert(exist(files{i}, 'file') == 2, 'KroneckerBio:LoadModelMassAction:FileNotFound', 'File %s was not found', files{i})
end

%% Loop over each file
for iFile = 1:nFiles
    % Open file
    fid = fopen(files{iFile});
    
    % Initialize
    mode = 0;
    lineNumber = 0;
    unitsRegexp = 'M|m|L';
    
    % Loop over all lines in code
    while true % break on -1
        % Get next line
        line = fgetl(fid);
        lineNumber = lineNumber + 1;
        
        % Break if end of file is encountered
        if line == -1
            break
        end
        
        % Skip if line is empty or begins with a "#"
        if isempty(line) || ~isempty(regexp(line, '^\s*#|^\s*$', 'once'))
            continue
        end
        
        % Check if mode switching line is encountered
        % First character is a "%"
        isSwitchMode = ~isempty(regexp(line, '^\s*%', 'once'));
        
        if isSwitchMode
            % Switch mode
            [match payload] = regexp(line, '^\s*%*\s*\w*', 'match', 'split', 'once');
            payload = payload{2}; % Only want the end, not the starting empty string
            match = regexp(match, '\w*', 'match');
            if strcmpi(match, 'compartments')
                % Switch mode
                mode = 1;
                
                % Extract model name
                modelName = regexp(payload, '".*"|[^.,\s]*', 'match', 'once');
                modelName = regexp(modelName, '[^"]*', 'match', 'once'); % Strip quotes
                m = RenameModel(m, modelName);
            elseif strcmpi(match, 'species')
                % Switch mode
                mode = 2;
                
                % Extract default compartment from payload
                defaultCompartment = regexp(payload, '".*"|[^.,\s]*', 'match');
                assert(numel(defaultCompartment) <= 1, 'KroneckerBio:LoadModelMassAction:MultipleSpeciesCompartments', 'Line %i in %s specified multiple compartments for a single section of species; a species section can have at most one default compartment: %s', lineNumber, files{iFile}, line)
                
                % Standardize compartment as a string
                if isempty(defaultCompartment)
                    defaultCompartment = '';
                else
                    defaultCompartment = defaultCompartment{1};
                end
            elseif strcmpi(match, 'outputs')
                % Switch mode
                mode = 3;
            elseif strcmpi(match, 'parameters')
                % Switch mode
                mode = 4;
            elseif strcmpi(match, 'reactions')
                % Switch mode
                mode = 5;

                % Extract default compartments from payload
                defaultCompartment = regexp(payload, '".*"|[^.,\s]*', 'match');
            elseif strcmpi(match, 'units')
                mode = 6;
            else
                error('KroneckerBio:LoadModelMassAction:UnrecognizedSectionHeader', 'Line %i in %s had an unrecognized section header: %s', lineNumber, files{iFile}, line)
            end
        else
            % Process entries
            if mode == 1
                % Compartments
                tokens = vec(regexp(line, ',|".*"|[^\s,]*','match'));
                assert(numel(tokens) >= 3, 'KroneckerBio:LoadModelMassAction:InsufficientCompartmentInformation', 'Line %i in %s does not provide the name, dimension, and volume required for a compartment: %s', lineNumber, files{iFile}, line)
                
                % Strip quotes
                for iTok = 1:numel(tokens)
                    tokens{iTok} = regexp(tokens{iTok}, '[^"]*', 'match', 'once');
                end
                
                % Extract compartment name
                assert(~isempty(regexp(tokens{1}, '^[^\s,.]*$', 'once')), 'KroneckerBio:LoadModelMassAction:InvalidName', 'Line %i in %s has an invalid compartment name: %s', lineNumber, files{iFile}, line)
                name = tokens{1};
                
                % Extract dimension
                assert(~isempty(regexp(tokens{2}, '^\d$', 'once')), 'KroneckerBio:LoadModelMassAction:InvalidCompartmentDimension', 'Line %i in %s has an invalid compartment dimension: %s', lineNumber, files{iFile}, line)
                dim = str2double(tokens{2});
                
                % Extract expressions into seperate cells
                expressions = cell(0,1);
                nExpr = 0;
                for i = 3:numel(tokens)
                    assert(~strcmp(tokens{i}, ','), 'KroneckerBio:LoadModelMassAction:UnexpectedToken', 'Line %i in %s has an unexpected comma: %s', lineNumber, files{iFile}, line)
                    
                    % Stop if this is the last token or the next token is
                    % not a comma
                    if i == numel(tokens) || ~strcmp(tokens{i+1}, ',')
                        expressions = tokens(3:i);
                        nExpr = numel(expressions);
                        break
                    else
                        % Remove the comma and continue extraction
                        tokens(i+1) = [];
                    end
                end
                
                % Parse each expression into regexp and value
                values = zeros(nExpr,1);
                units = cell(nExpr,1);
                for i = 1:nExpr
                    % Split between the equal sign
                    tokens = regexp(expressions{i}, '=', 'split');
                    assert(numel(tokens) <= 2, 'KroneckerBio:LoadModelMassAction:InvalidExpression', 'Line %i in %s has more than one equal sign in an expression: %s', lineNumber, files{iFile}, line)
                    
                    if numel(tokens) == 1
                        if ~isempty(regexp(tokens{1}, '^\d', 'once'))
                            % If it starts with a digit, it is a
                            % constituitive volume
                            expressions{i} = '';
                            
                            % Match the number, the rest are units
                            [value unit] = regexp(tokens{1}, '^[-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?', 'match', 'split', 'once');
                            values(i) = str2double(value);
                            units{i} = unit{2};
                        else
                            % Otherwise it is a regular expression with
                            % value equal to the default of 1
                            expressions{i} = tokens{1};
                            values(i) = 1;
                            units{i} = '';
                        end
                    else%if numel(tokens) == 2
                        % Store expression
                        expressions{i} = tokens{1};
                        
                        % Match the number, the rest are units
                        [value unit] = regexp(tokens{2}, '^[-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?', 'match', 'split', 'once');
                        values(i) = str2double(value);
                        units{i} = unit{2}; % The first had better be zero
                    end
                end
                
                % Add compartment
                m = AddCompartment(m, name, dim, expressions, values, units);
            elseif mode == 2
                % Species
                tokens = vec(regexp(line, ',|".*"|[^\s,]*','match'));
                
                % Strip quotes
                for iTok = 1:numel(tokens)
                    tokens{iTok} = regexp(tokens{iTok}, '[^"]*', 'match', 'once');
                end
                
                % Extract species name
                assert(~isempty(regexp(tokens{1}, '^[^\s,]*\.?[^\s,]*$', 'once')), 'KroneckerBio:LoadModelMassAction:InvalidName', 'Line %i in %s has an invalid species name: %s', lineNumber, files{iFile}, line)
                name = tokens{1};
                if isempty(regexp(name, '\.', 'once'))
                    % Compartment is defualt compartment
                    assert(~isempty(defaultCompartment), 'KroneckerBio:LoadModelMassAction:SpeciesWithNoCompartment', 'Line %i in %s has a species that is not associated with a compartment: %s', lineNumber, files{iFile}, line)
                    compartment = defaultCompartment;
                else
                    % Compartment is specified within name
                    compartment = regexp(name, '^[^.]*', 'match', 'once');
                    name = regexp(name, '[^.]*$', 'match', 'once');
                end
                
                % Extract expressions into seperate cells
                value = cell(0,1);
                nVal = 0;
                for i = 2:numel(tokens)
                    assert(~strcmp(tokens{i}, ','), 'KroneckerBio:LoadModelMassAction:UnexpectedToken', 'Line %i in %s has an unexpected comma: %s', lineNumber, files{iFile}, line)
                    
                    % Stop if this is the last token or the next token is
                    % not a comma
                    if i == numel(tokens) || ~strcmp(tokens{i+1}, ',')
                        value = tokens(2:i);
                        nVal = numel(value);
                        break
                    else
                        % Remove the comma and continue extraction
                        tokens(i+1) = [];
                    end
                end
                
                % Extract whether or not species is an input
                if numel(tokens) > 1 + nVal
                    % Input is specified
                    if strcmpi(tokens{2+nVal},'true') || strcmpi(tokens{2+nVal},'t') || strcmp(tokens{2+nVal},'1')
                        isInput = true;
                    elseif strcmpi(tokens{2+nVal},'false') || strcmpi(tokens{2+nVal},'f') || strcmp(tokens{2+nVal},'0')
                        isInput = false;
                    else
                        error('KroneckerBio:LoadModelMassAction:InvalidIsInput', 'Line %i in %s does not properly specify whether or not the species is an input, must be false/true or 0/1: %s', lineNumber, files{iFile}, line)
                    end
                else
                    % Input was not specified, default to state
                    isInput = false;
                end
                
                % Process values
                if isInput
                    % Parse input by getting each input string ready for
                    % parsing
                    for i = 1:numel(value)
                        % Find first substring matching a unit. Expression
                        % is everything before, units are everything after.
                        [expression unit] = regexp(value{i}, unitsRegexp, 'split', 'match', 'once');
                        expression = expression{1};
                        % Only keep the units on the first expression
                        if i == 1
                            units = unit;
                        end
                        if ~isempty(regexp(expression, '^[-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?$', 'match', 'once'))
                            % Expression is a constant
                            value = str2double(expression);
                            break
                        else
                            % Expression is something more complicated
                            if expression(1) == '@'
                                % It is already a function handle
                                value{i} = expression;
                            else
                                % It needs to be made into a function handle
                                value{i} = ['@(t,q)' expression];
                            end
                        end
                    end
                    
                    % Process parameters
                    parameters = zeros(0,1);
                    if numel(tokens) > 2+nVal
                        % There are parameters
                        for i = 3+nVal:numel(tokens)
                            assert(~strcmp(tokens{i}, ','), 'KroneckerBio:LoadModelMassAction:UnexpectedToken', 'Line %i in %s has an unexpected comma: %s', lineNumber, files{iFile}, line)
                            % Stop if this is the last token or the next token is
                            % not a comma
                            if i == numel(tokens) || strcmp(tokens{i}, ',')
                                parameters = tokens(3+nVal:i);
                                break
                            else
                                % Remove the comma and continue extraction
                                tokens(i+1) = [];
                            end
                        end
                        parameters = str2double(parameters);
                    end
                else
                    if isempty(value)
                        % Default initial value of 0
                        value = 0;
                        units = '';
                        parameters = [];
                    else
                        % Match the number, the rest are units
                        [value units] = regexp(value{1}, '^[-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?', 'match', 'split', 'once');
                        value = str2double(value);
                        units = units{2};
                        parameters = [];
                    end
                end
                
                
                % Add Species
                m = AddSpecies(m, name, compartment, value, units, isInput, parameters);
            elseif mode == 3
                % Outputs
                tokens = vec(regexp(line, ',|".*"|[^\s,]*','match'));
                
                % Strip quotes
                for iTok = 1:numel(tokens)
                    tokens{iTok} = regexp(tokens{iTok}, '[^"]*', 'match', 'once');
                end
                
                % Extract output name
                assert(~isempty(regexp(tokens{1}, '^[^\s,.]*$', 'once')), 'KroneckerBio:LoadModelMassAction:InvalidName', 'Line %i in %s has an invalid output name: %s', lineNumber, files{iFile}, line)
                name = tokens{1};
                
                % Extract expressions into seperate cells
                expressions = cell(0,1);
                nExpr = 0;
                for i = 2:numel(tokens)
                    assert(~strcmp(tokens{i}, ','), 'KroneckerBio:LoadModelMassAction:UnexpectedToken', 'Line %i in %s has an unexpected comma: %s', lineNumber, files{iFile}, line)
                    
                    % Stop if this is the last token or the next token is
                    % not a comma
                    if i == numel(tokens) || ~strcmp(tokens{i+1}, ',')
                        expressions = tokens(2:i);
                        nExpr = numel(expressions);
                        break
                    else
                        % Remove the comma and continue extraction
                        tokens(i+1) = [];
                    end
                end
                
                % Parse each expression into regexp and value
                values = zeros(nExpr,1);
                units = cell(nExpr,1);
                for i = 1:nExpr
                    % Split between the equal sign
                    tokens = regexp(expressions{i}, '=', 'split');
                    assert(numel(tokens) <= 2, 'KroneckerBio:LoadModelMassAction:InvalidExpression', 'Line %i in %s has more than one equal sign in an expression: %s', lineNumber, files{iFile}, line)
                    
                    if numel(tokens) == 1
                        if ~isempty(regexp(tokens{1}, '^\d', 'once'))
                            % If it starts with a digit, it is a
                            % constituitive volume
                            expressions{i} = '';
                            
                            % Match the number, the rest are units
                            [value unit] = regexp(tokens{1}, '^[-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?', 'match', 'split', 'once');
                            values(i) = str2double(value);
                            units{i} = unit{2};
                        else
                            % Otherwise it is a regular expression with
                            % value equal to the default of 1
                            expressions{i} = tokens{1};
                            values(i) = 1;
                            units{i} = '';
                        end
                    else%if numel(tokens) == 2
                        % Store expression
                        expressions{i} = tokens{1};
                        
                        % Match the number, the rest are units
                        [value unit] = regexp(tokens{2}, '^[-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?', 'match', 'split', 'once');
                        values(i) = str2double(value);
                        units{i} = unit{2}; % The first had better be zero
                    end
                end
                
                % Add output
                m = AddOutput(m, name, expressions, values, units);
            elseif mode == 4
                % Parameters
                tokens = vec(regexp(line, ',|".*"|[^\s,]*','match'));
                
                % Strip quotes
                for iTok = 1:numel(tokens)
                    tokens{iTok} = regexp(tokens{iTok}, '[^"]*', 'match', 'once');
                end
                
                % Extract parameter name
                assert(~isempty(regexp(tokens{1}, '^[^\s,.]*$', 'once')), 'KroneckerBio:LoadModelMassAction:InvalidName', 'Line %i in %s has an invalid compartment name: %s', lineNumber, files{iFile}, line)
                name = tokens{1};
                
                % Match the number, the rest are units
                [value units] = regexp(tokens{2}, '^[-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?', 'match', 'split', 'once');
                value = str2double(value);
                units = units{2}; % The first had better be zero
                
                % Add parameter
                m = AddParameter(m, name, value, units);
            elseif mode == 5
                % Reactions
                tokens = vec(regexp(line, ',|".*"|[^\s,]*','match'));
                
                % Strip quotes
                for iTok = 1:numel(tokens)
                    tokens{iTok} = regexp(tokens{iTok}, '[^"]*', 'match', 'once');
                end
                
                % Extract names
                name = cell(0,1);
                if numel(tokens) > 6
                    name{1,1} = tokens{7};
                end
                if numel(tokens) >= 8 && strcmp(tokens{8}, ',')
                    assert(numel(tokens) >= 9, 'KroneckerBio:LoadModelMassAction:UnexpectedToken', 'Line %i in %s contains an unsatisfied comma: %s', lineNumber, files{iFile}, line)
                    assert(numel(tokens) <= 9 || tokens{10} ~= ',', 'KroneckerBio:LoadModelMassAction:TooManyReactionNames', 'Line %i in %s supplies more than two reaction names: %s', lineNumber, files{iFile}, line)
                    name{2,1} = tokens{9};
                end
                
                m = AddReaction(m, name, defaultCompartment, tokens{1}, tokens{2}, tokens{3}, tokens{4}, tokens{5}, tokens{6});
            elseif mode == 6
                % Units
            else
                error('KroneckerBio:LoadModelMassAction:SectionNotSpecified', 'Line %i in %s occurs before any section header: %s', lineNumber, files{iFile}, line)
            end
        end
    end
    
    % Close the model file
    fclose(fid);
end

%% Finalize Model
m = FinalizeModel(m);