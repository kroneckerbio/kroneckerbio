function obj = Update(obj, varargin)
assert(isscalar(obj), 'KroneckerBio:Update', 'This function only accepts scalar objects')

if any(regexp(obj.Type, '^Model'))
    % Update a model
    if nargin <= 1
        for i = 1:numel(obj)
            obj(i) = obj.Update(obj(i).k, obj(i).s, obj(i).q);
        end
    elseif nargin <= 2
        if ~iscell(varargin{1})
            % Encapsulate
            varargin{1} = varargin(1);
        end
        for i = 1:numel(obj)
            obj(i) = obj.Update(varargin{1}{i}, obj(i).s, obj(i).q);
        end
    elseif nargin <= 3
        if ~iscell(varargin{1})
            % Encapsulate
            varargin{1} = varargin(1);
        end
        if ~iscell(varargin{2})
            % Encapsulate
            varargin{2} = varargin(2);
        end
        for i = 1:numel(obj)
            obj(i) = obj.Update(varargin{1}{i}, varargin{2}{i}, obj(i).q);
        end
    elseif nargin <= 4
        if ~iscell(varargin{1})
            % Encapsulate
            varargin{1} = varargin(1);
        end
        if ~iscell(varargin{2})
            % Encapsulate
            varargin{2} = varargin(2);
        end
        if ~iscell(varargin{3})
            % Encapsulate
            varargin{3} = varargin(3);
        end
        for i = 1:numel(obj)
            obj(i) = obj.Update(varargin{1}{i}, varargin{2}{i}, varargin{3}{i});
        end
    end
elseif any(regexp(obj.Type, '^Experiment'))
    if nargin <= 1
        obj = obj.Update(obj.s, obj.q);
    elseif nargin <= 2
        obj = obj.Update(varargin{1}, obj.q);
    elseif nargin <= 3
        obj = obj.Update(varargin{1}, varargin{2});
    end
else
    error('KroneckerBio:Update', 'That object type cannot be updated with Update at this time.')
end
