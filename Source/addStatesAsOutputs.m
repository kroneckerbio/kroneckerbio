function m = addStatesAsOutputs(m)
%addStatesAsOutputs A quick and dirty helper script that adds one output
%   for each state currently in the model which tracks the concentration of
%   that state
%
%   m = addStatesAsOutputs(m)

% (c) 2015 David R Hagen
% This work is released under the MIT license.

all_states = [m.States; m.add.States];

full_names = strcat(vec({all_states.Compartment}), '.', vec({all_states.Name}));
full_names_regex = strcat('^', vec({all_states.Name}), '$');

for i = 1:numel(full_names)
    m = AddOutput(m, full_names{i}, full_names_regex{i});
end
