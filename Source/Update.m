function obj = Update(obj, varargin)
assert(isscalar(obj), 'KroneckerBio:Update', 'This function only accepts scalar objects')

if any(regexp(obj.Type, '^Model'))
    % Update a model
    if nargin <= 1
        for i = 1:numel(obj)
            obj(i) = obj.Update(obj(i).k);
        end
    elseif nargin <= 2
        if ~iscell(varargin{1})
            % Encapsulate
            varargin{1} = varargin(1);
        end
        for i = 1:numel(obj)
            obj(i) = obj.Update(varargin{1}{i});
        end
    end
elseif any(regexp(obj.Type, '^Experiment'))
    if nargin <= 1
        obj = obj.Update(obj.s, obj.q, obj.h);
    elseif nargin <= 2
        obj = obj.Update(varargin{1}, obj.q, obj.h);
    elseif nargin <= 3
        obj = obj.Update(varargin{1}, varargin{2}, obj.h);
    elseif nargin <= 4
        obj = obj.Update(varargin{1}, varargin{2}, varargin{3});
    end
else
    error('KroneckerBio:Update', 'That object type cannot be updated with Update at this time.')
end
