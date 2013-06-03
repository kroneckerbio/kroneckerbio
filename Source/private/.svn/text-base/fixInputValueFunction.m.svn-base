function [func handled] = fixInputValueFunction(value, units, v0)
% TODO: why did I put in v0?
if nargin < 2
    units = [];
    if nargin < 1
        value = [];
    end
end

if isempty(value)
    value = 0;
end

% TODO: units
if ~isempty(units); error('Units not yet implemented, tell David if you would like to volunteer'); end

% TODO: deal with the derivatives
if iscell(value); value = value{1}; end

handled = true;

if isnumeric(value)
    func = eval(sprintf('@(t,q)repmat(%g, 1,numel(t))', value));
elseif isa(value, 'function_handle')
    func = value;
elseif ischar(value)
    if value(1) =='@'
        % Already a function handle
        func = eval(value);
    else
        % Needs the handle
        value = ['@(t,q)' value];
        func = eval(value);
    end
else
    error('KroneckerBio:InputValueFunction', 'The class of the input value function was invalid')
end
