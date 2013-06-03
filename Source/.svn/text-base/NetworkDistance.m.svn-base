function net = NetworkDistance(m, maxn, reflect)
%NetworkDistance Compute the distance that any two species in a network are
%   from each other.
%
%   net = NetworkDistance(m, maxn, reflect)
%
%   The distance is computed on model m. The network net is a matrix nx by
%   nx. The matrix tells the distance between the species in the columns
%   and the species in the rows. For example, a 1 at position (3,2) means
%   that the value of species 2 affects the ODE of species 3. If there are
%   irreversible reactions, such that the amount of 3 affects 2 but the
%   amount of 2 never affects the ODE of 3, the matrix will not be
%   symmetric. Therefore, the distance from 3 to 2 may be different from
%   the distance from 2 to 3.
%
%   The positive integer maxn puts a bound on how deep to calculate the
%   distance. For example, setting maxn to 1 will only find the adjacent
%   species. The maximum value is the default nx, because no two species
%   can be more than nx reactions away from each other and still be
%   connected.
%
%   For species that never affect each other, or are more than maxn away
%   from each other, the value in net is 0.
%
%   The logical reflect makes the distances symmetric by setting the
%   distances in both directions equal to their minimum. Thus, 2 is
%   connected to 3 by a distance of 1 even if 2 is not a factor in 3's ODE.
%
%   Note: Reactions with parameters equal to zero are the same as reactions
%   that do not exist. If you want to consider the network distance of the
%   actual topology, set your parameters to non-zero values.
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model
%   maxn: [ positive integer scalar ]
%       A bound on how deep to calculate the network distance. States that
%       are more than this integer away, will be returned as 0.
%   reflect: [ logical scalar ]
%       Symmetrize the distance matrix by taking the minimum distance
%       between two states.
%
%   Outputs
%   net: [ nonnegative integer matrix nx by nx ]
%       The distances between each state (columns) and each other state
%       (rows)

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% TODO: add ability to include dfdu as well
% Clean up inputs
assert(nargin >= 1, 'KroneckerBio:NetworkDistance:TooFewInputs', 'NetworkDistance requires at least 1 input arguments')
if nargin < 3
    reflect = [];
    if nargin < 2
        maxn = [];
    end
end

assert(isscalar(m), 'KroneckerBio:NetworkDistance:MoreThanOneModel', 'The model structure must be scalar')

% Constants
nx = m.nx;
nu = m.nu;

% Defaults
if isempty(maxn)
    maxn = nx;
end
if isempty(reflect)
    reflect = false;
end

% Fetch dfdx
dfdx = m.dfdx(rand,rand(nx,1),rand(nu,1));

% Convert to ones to get interaction matrix
dfdx = spones(dfdx);
net = dfdx;

% Recursively multiply to get increasing distance
for i = 2:maxn
    % These are the old, dominant indexes
    inds = (net ~= 0);
    
    % Compute the new interactions
    next = dfdx*net;
    next = spones(next)*i;
    
    % Paste the dominate indexes on top of the new ones
    next(inds) = net(inds);
    net = next;
end

% Create two-way interaction matrix if requested
if reflect
    net = min(net,net');
end

