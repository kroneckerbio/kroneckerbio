function u = completeInputFunction(inputs)
% This is a helper function that turns a bunch of input functions from
% different inputs (m.Inputs.Function) into a single function (m.u). This
% function is in /Source rather than /Source/private because it needs to be
% in the path to load models that are saved.

nu = numel(inputs);

% Complete function
if nu == 0
    u = @(t,q)zeros(0,numel(t));
else
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