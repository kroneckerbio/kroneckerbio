function u = completeInputFunction(inputs)

nu = numel(inputs);

% Complete function
if nu == 0
    u = @(t,q)zeros(0,numel(t));
else
    % Extract input values
    inputs = cat(1, inputs.Value);
    
    % Stack each input function
    uEach = {inputs.Function}.';
    
    % Stack the parameters in a cell vector
    q = {inputs.Parameters}.';
    u = @(t)inputStack(t);
end

    function val = inputStack(t)
        % Evaluate each individual input function
        val = cellfun(@feval, uEach, repmat({t}, nu,1), q, 'UniformOutput', false);
        
        % Collapse the cell array into a matrix
        val = cat(1, val{:});
    end
end