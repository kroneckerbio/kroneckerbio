function states = growStatesAnalytic(states, nx)

if nargin < 2
    nx = 0;
    if nargin < 1
        states = [];
    end
end

current = numel(states);
add = nx - current;
if add > 0 || isempty(states)
    % Double length
    add = max(current,add);
    states = [states; struct('Name', cell(add,1), 'ID', cell(add,1), 'Compartment', cell(add,1), 'InitialValue', cell(add,1))];
end