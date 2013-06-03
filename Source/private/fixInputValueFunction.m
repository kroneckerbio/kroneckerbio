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
    func = eval(sprintf('@(t,q)repmat(%g, 1,numel(t))', func));
elseif isa(func, 'function_handle')
    func = func;
elseif ischar(func)
    if func(1) =='@'
        % Already a function handle
        func = eval(func);
    else
        % Needs the handle
        func = ['@(t,q)' func];
        func = eval(func);
    end
else
    error('KroneckerBio:Input:Function', 'The input function was invalid')
end
