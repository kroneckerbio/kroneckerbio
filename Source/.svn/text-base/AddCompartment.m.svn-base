function m = AddCompartment(m, name, dimension, expressions, values, units)
%AddCompartment Add a compartment to a KroneckerBio model
%
%   m = AddCompartment(m, name, dimension, expressions, values, units)
%
%   Compartments hold species and have volume. The volume is a function of
%   species values.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the compartment will be added
%   name: [ string ]
%       A name for the compartment
%   dimension: [ 0 1 2 3 ]
%       The dimensionality of the compartment. Example: the cytoplasm would
%       be 3, the plasma membrane or the nuclear membrane would be 2. This
%       is only used to determine which compartment's volume plays a part
%       in the rate of a bimolecular reaction between compartments.
%   expressions: [ string | cell array of strings ]
%       Each species full name will be substring matched against each
%       regular expression provided here. These species contribute to the
%       volume. An empty string indicates that the compartment has a
%       constant component to its volume.
%   values: [ nonnegative scalar | nonnegative vector numel(expressions)
%             {0} ]
%       Optional
%       Each entry in this vector is associated with one of the
%       expressions. This value tells how much a species will contribute
%       when it matches the corresponding expressions. If a species is
%       matched by multiple expressions, the last expression to match
%       overrides all others.
%   units: [ '' ]
%       This is currently not implemented.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new compartment added.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 6
    units = [];
end

% Validate compartment name
assert(isempty(regexp(name, '\.|,', 'once')), 'KroneckerBio:AddCompartment:CompartmentName', 'Compartment names cannot contain "." or ","')

% Increment counter
nv = m.add.nv + 1;
m.add.nv = nv;
m.add.Compartments = growCompartments(m.add.Compartments, m.add.nv);

% Add item
m.add.Compartments(nv).Name = fixName(name);
m.add.Compartments(nv).Expressions = fixExpression(expressions);
m.add.Compartments(nv).Values = fixCompartmentValues(values, fixUnits(units));
m.add.Compartments(nv).Dimension = dimension;

m.Ready = false;