function extractSpecies(this)

import Util.cellstr2str

%% Trim m.* components to only those actually added
this.Compartments = this.Compartments(1:this.nv);
this.Seeds        = this.Seeds(1:this.ns);
this.Parameters   = this.Parameters(1:this.nk);
this.Inputs       = this.Inputs(1:this.nu);
this.States       = this.States(1:this.nx);
this.Reactions    = this.Reactions(1:this.nr);
this.Rules        = this.Rules(1:this.nz);
this.Outputs      = this.Outputs(1:this.ny);

nx = this.nx;
nu = this.nu;
nxu = nx + nu;

%% Extract names
v_names = vec({this.Compartments.Name});
k_names = vec({this.Parameters.Name});
s_names = vec({this.Seeds.Name});
u_names = vec({this.Inputs.Name});
x_names = vec({this.States.Name});
r_names = vec({this.Reactions.Name});
z_names = vec({this.Rules.Name});
y_names = vec({this.Outputs.Name});
xu_names = [x_names; u_names];

% Make list of all compartment.species in model
x_full_names = vec(strcat({this.States.Compartment}, '.', {this.States.Name}));
u_full_names = vec(strcat({this.Inputs.Compartment}, '.', {this.Inputs.Name}));
xu_full_names = [x_full_names; u_full_names];

% Make species - compartment index mapping
vx_names = vec({this.States.Compartment});
vxInd = lookupmember(vx_names, v_names);
vu_names = vec({this.Inputs.Compartment});
vuInd = lookupmember(vu_names, v_names);

%% Error for repeated components
[~, ia, ~] = unique(v_names);
v_repeated = v_names(setdiff(1:this.nv, ia));
assert(isempty(v_repeated), 'KroneckerBio:FinalizeModel:RepeatCompartment', 'Compartment %s not unique', cellstr2str(v_repeated))

[~, ia, ~] = unique(k_names);
k_repeated = k_names(setdiff(1:this.nk, ia));
assert(isempty(k_repeated), 'KroneckerBio:FinalizeModel:RepeatParameter', 'Parameter %s not unique', cellstr2str(k_repeated))

[~, ia, ~] = unique(s_names);
s_repeated = s_names(setdiff(1:this.ns, ia));
assert(isempty(s_repeated), 'KroneckerBio:FinalizeModel:RepeatSeed', 'Seed %s not unique', cellstr2str(s_repeated))

[~, ia, ~] = unique(xu_full_names);
xu_repeated = xu_full_names(setdiff(1:nxu, ia));
assert(isempty(xu_repeated), 'KroneckerBio:FinalizeModel:RepeatSpecies', 'Species %s not unique', cellstr2str(xu_repeated))

[~, ia, ~] = unique(y_names);
y_repeated = y_names(setdiff(1:this.ny, ia));
assert(isempty(y_repeated), 'KroneckerBio:FinalizeModel:RepeatOutput', 'Output %s not unique', cellstr2str(y_repeated))

[~, ia, ~] = unique(z_names);
z_repeated = z_names(setdiff(1:this.nz, ia));
assert(isempty(z_repeated), 'KroneckerBio:FinalizeModel:RepeatRule', 'Rule %s not unique', cellstr2str(z_repeated))

% Reactions are checked for uniqueness of names here, not content
% [~, ia, ~] = unique(r_names);
% r_repeated = r_names(setdiff(1:this.nr, ia));
% assert(isempty(r_repeated), 'KroneckerBio:FinalizeModel:RepeatReactionName', 'Reaction name %s not unique', cellstr2str(r_repeated))

% Make logical vector of which species names are unique
unique_xu_names = false(nxu,1);
for ixu = 1:nxu
    unique_xu_names(ixu) = ~ismember(xu_names{ixu}, [xu_names(1:ixu-1); xu_names(ixu+1:end)]);
end
unique_x_names = unique_xu_names(1:nx);
unique_u_names = unique_xu_names(nx+1:nxu);

% Check for uniquness of names in entire model
% A name can be duplicated within species, but cannot otherwise be duplicated
% Species full names cannot be duplicated, though
assert_unique_name([v_names; s_names; k_names; unique(u_names); unique(x_names); z_names])
assert_unique_name(xu_full_names)

% Only unique species names can be referred to by their unqualified names
% Unique species names and full names point to the same symbols
% ambiguous_names = [u_names(~unique_u_names); x_names(~unique_x_names)];
all_names = [v_names; s_names; k_names; u_names(unique_u_names); x_names(unique_x_names); z_names; u_full_names; x_full_names];

%% Assign extracted components to m
m = this.m;

m.v_names = v_names;
m.vx_names = vx_names;
m.vu_names = vu_names;
m.vxInd = vxInd;
m.vuInd = vuInd;

m.xu_full_names = xu_full_names;

m.unique_x_names = unique_x_names;
m.unique_u_names = unique_u_names;
m.all_names = all_names;

this.m = m;

end

function assert_unique_name(names)
names = vec(names);
n_new = numel(names);

for i = 1:n_new
    % Check if item by this name already exists
    if ~isempty(names)
        name = names{i};
        if any(strcmp(name, names(1:i-1)))
            error('KroneckerBio:FinalizeModel:RepeatName', 'There are multiple components with the name %s', name)
        end
    end
end
end