function func = fixInputValueFunction(func)
% Standardize the input function as @(t,q)

if nargin < 1
    func = [];
end

if isempty(func)
    func = 0;
end

% TODO: deal with the derivatives
if iscell(func); func = func{1}; end

if isnumeric(func) && isscalar(func) && func >= 0
    % Function is a numeric scalar
    func = eval(sprintf('@(t,q)repmat(%g, 1,numel(t))', func));
elseif isa(func, 'function_handle')
    func = func;
elseif ~isempty(regexp(func, '^[+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?$', 'match', 'once'))
    % Expression is a string representation of a scalar
    func = eval(['@(t,q)repmat(' func ', 1,numel(t))']);
elseif ischar(func)
    if func(1) =='@'
        % Already a function handle
        func = eval(func);
    else
        % Needs the handle
        func = ['@(t,q)' func '+row(t).*0'];
        func = eval(func);
    end
else
    error('KroneckerBio:Input:Function', 'The input function was invalid')
end
