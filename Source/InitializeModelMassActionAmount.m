function m = InitializeModelMassActionAmount(name)
%InitializeModelMassActionAmount Create an empty mass action Kronecker
%   model structure which can then have compartment, species, outputs,
%   parameters, and reactions added to it.
%
%   m = InitializeModel(name)
%
%   Inputs
%   name: [ string {''} ]
%       The name of the new Kronecker model
%
%   Outputs
%   m: [ model struct scalar ]
%       This structure is a Kronecker model. It has fields that organize
%       the properties of the model, summarize the constants of the model,
%       or provide handles to functions that describe the changing state of
%       the system.
%       .Type 'Model.MassActionKronecker'
%       .Name [ string ]
%       .Compartments [ struct vector ]
%           Each compartment is described by an element of this vector
%           .Name [ string ]
%               The name of the compartment
%           .Dimension [ 0 1 2 3 ]
%               The dimensionality of the compartment
%           .Size [ nonnegative scalar ]
%               The volume of the compartment
%       .Parameters [ struct vector ]
%           .Name [ string ]
%               The name of the kinetic parameter
%           .Value [ nonnegative scalar ]
%               The macroscopic value of the kinetic parameter
%       .Seeds [ struct vector ]
%           .Name
%               The name the seed parameter
%           .Value [ nonnegative scalar ]
%               The value of the seed parameter
%       .Inputs [ struct vector ]
%           .Name [ string ]
%               The name of the input species
%           .Compartment [ string ]
%               The name of the compartment in which this input resides
%           .DefaultValue [ nonnegative scalar ]
%               The value assigned when an experiment is constructed but no
%               value is given for the inputs
%       .States [ struct vector ]
%           .Name [ string ]
%               The name of the state species
%           .Compartment [ string ]
%               The name of the compartment in which this species resides
%           .InitialValue [ cell matrix ]
%               This matrix is ? by 2. Each row is a string followed by a
%               nonnegative scalar. The string represents one seed and the
%               scalar represents the amount that the seed contributes to
%               this state's initial value.
%       .Reactions [ struct vector ]
%           Each reaction is elementary, no reaction in this vector
%           corresponds to both the forward and reverse reverse of a
%           reversible reaction.
%           .Name [ string ]
%               The name of the reaction. This can be empty and is only
%               used by the model modification functions to alter or remove
%               reactions by name.
%           .Reactants [ cell vector 1 by 2 of strings ]
%               Each string corresponds to the full name a specific species
%               that is a reactant in this reaction. If a string is
%               empty, it means that no species takes part.
%           .Products [ cell vector 1 by 2 of strings ]
%               Each string corresponds to the full name a specific species
%               that is a product of this reaction. If a string is
%               empty, it means that no species is created.
%           .Parameter [ cell 1 by 2 ]
%               The first element is a string name of the kinetic parameter
%               associated with this reaction. The second element is a
%               scalar that is the modifier on this parameter for this
%               reaction.
%       .Outputs [ struct vector ]
%           .Name [ string ]
%               The name of the output
%           .Expressions [ cell vector of regexp strings | '' ]
%               The regular expressions match species of the model. Each
%               species whose full name contains the described substring
%               contributes to the output an amount equal to the
%               corresponding element in Values. An empty string means
%               that the value is constituitive (the output simply has
%               that value regardless of any species value).
%           .Values [ nonnegative vector ]
%               The value of the output is the sum of each value in
%               this vector times the sum of the amounts of each species
%               that is matched by Expressions.
%       .nv [ whole scalar ]
%           The number of compartments
%       .nk [ whole scalar ]
%           The number of kinetic parameters
%       .ns [ whole scalar ]
%           The number of seed parameters
%       .nu [ whole scalar ]
%           The number of input species
%       .nx [ whole scalar ]
%           The number of state species
%       .nr [ whole scalar ]
%           The number of elementary reactions
%       .ny [ whole scalar ]
%           The number of outputs
%       .dv [ whole vector nv ]
%           The dimensions of each compartment
%       .k [ positive vector nk ]
%           The values of each kinetic parameter
%       .s [ nonegative vector ns ]
%           The values of the the seed parameters
%       .u [ nonnegative vector nu ]
%           The vector of default values
%       .vxInd [ natural vector nx ]
%           The indexes to the compartments containing each state species
%       .vuInd [ natural vector nu ]
%           The indexes to the compartments containing each input species
%       .rOrder [ 0 1 2 vector nr ]
%           The order of each reaction
%       .krInd [ natural vector nr ]
%           The indexes to the parameters associated with each reaction
%       .x0 [ function handle returning nonnegative vector nx ]
%           The values of the initial conditions of each state species.
%       .dx0ds [ function of seeds returning nonnegative matrix nx by ns ]
%           First derivative of the initial conditions with respect to the
%           seeds. Constant valued with respect to s for mass action
%           models.
%       .d2x0ds2 [ function of seeds returning nonnegative matrix
%                  nx*ns-by-ns]
%           Second derivative of the initial conditions with respect to the
%           seeds. Zero matrix for mass action models.
%       .A1 [ real matrix nx by nx]
%           Map of how each state contributes to the mass action ODEs
%       .A2 [ real matrix nx by nx*nx ]
%           Map of how each product of states in x kron x/vx contributes to
%           the mass action ODEs
%       .A3 [ real matrix nx by nu*nx ]
%           Map of how each product in u kron x/vx  contributes to
%           the mass action ODEs
%       .A4 [ real matrix nx by nx*nu ]
%           Map of how each product in x kron u/vu contributes to
%           the mass action ODEs
%       .A5 [ real matrix nx by nu*nu ]
%           Map of how each product of inputs in u kron u/vu contributes to
%           the mass action ODEs
%       .A6 [ real matrix nx by nu ]
%           Map of how each input contributes to the mass action ODEs
%       .a [ nonnegative vector nx]
%           Contribution of constituitive synthesis to the mass action ODEs
%       .B1 [ nonnegative matrix nv by nx ]
%           Contribution of each state species to the volume of each
%           compartment
%       .B2 [ nonnegative matrix nv by nu ]
%           Contribution of each input species to the volume of each
%           compartment
%       .b [ nonnegative vector nv ]
%           Constituitive volume of each compartment
%       .dA1dk [ real matrix nx*nx by nk] (fx_k)
%           Derivative of A1 wrt k
%       .dA2dk [ real matrix nx*nx*nx by nk ] (fxx_k)
%           Derivative of A2 wrt k
%       .dA3dk [ real matrix nx*nu*nx by nk ] (fux_k)
%           Derivative of A3 wrt k
%       .dA4dk [ real matrix nx*nx*nu by nk ] (fxu_k)
%           Derivative of A4 wrt k
%       .dA5dk [ real matrix nx*nu*nu by nk ] (fuu_k)
%           Derivative of A5 wrt k
%       .dA6dk [ real matrix nx*nu by nk ] (fu_k)
%           Derivative of A6 wrt k
%       .dadk [ real matrix nx by nk ] (f_k)
%           Derivative of a wrt k
%       .dA1dk_fk_x [ real matrix nx*nk by nx ] (fk_x)
%           Reshape of dA1dk for faster computation of dfdk
%       .dA2dk_fk_xx [ real matrix nx*nk by nx*nx ] (fk_xx)
%           Reshape of dA2dk for faster computation of dfdk
%       .dA3dk_fk_ux [ real matrix nx*nk by nu*nx ] (fk_ux)
%           Reshape of dA3dk for faster computation of dfdk
%       .dA4dk_fk_xu [ real matrix nx*nk by nx*nu ] (fk_xu)
%           Reshape of dA4dk for faster computation of dfdk
%       .dA5dk_fk_uu [ real matrix nx*nk by nu*nu ] (fk_uu)
%           Reshape of dA5dk for faster computation of dfdk
%       .dA6dk_fk_u [ real matrix nx*nk by nu ] (fk_u)
%           Reshape of dA6dk for faster computation of dfdk
%       .f [ handle @(t,x,u) returns real vector nx ]
%           This is the ODE function of the system. The time t is scalar,
%           but not used since mass action systems are constant coefficient
%           differential equations. The state x is a nonnegative vector nx.
%           The input function u is the handle above.
%       .dfdx [ handle @(t,x,u) returns real matrix nx by nx ]
%           The partial derivative of f wrt x
%       .dfdu [ handle @(t,x,u) returns real matrix nx by nu ]
%           The partial derivative of f wrt u
%       .dfdk [ handle @(t,x,u) returns real matrix nx by nk ]
%           The partial derivative of f wrt k
%       .d2fdx2 [ handle @(t,x,u) returns real matrix nx*nx by nx ] (fx_x)
%           The partial derivative of dfdx wrt x
%       .d2fdu2 [ handle @(t,x,u) returns real matrix nx*nu by nu ] (fu_u)
%           The partial derivative of dfdu wrt u
%       .d2fdk2 [ handle @(t,x,u) returns real matrix nx*nk by nk ] (fk_k)
%           The partial derivative of dfdk wrt k
%       .d2fdudx [ handle @(t,x,u) returns real matrix nx*nx by nu ] (fx_u)
%           The partial derivative of dfdx wrt u. 
%       .d2fdxdu [ handle @(t,x,u) returns real matrix nx*nu by nx ] (fu_x)
%           The partial derivative of dfdu wrt x. 
%       .d2fdkdx [ handle @(t,x,u) returns real matrix nx*nx by nk ] (fx_k)
%           The partial derivative of dfdx wrt k. 
%       .d2fdxdk [ handle @(t,x,u) returns real matrix nx*nk by nx ] (fk_x)
%           The partial derivative of dfdk wrt x.
%       .d2fdkdu [ handle @(t,x,u) returns real matrix nx*nu by nk ] (fu_k)
%           The partial derivative of dfdu wrt k. 
%       .d2fdudk [ handle @(t,x,u) returns real matrix nx*nk by nu ] (fk_u)
%           The partial derivative of dfdk wrt u. 
%       .S [ integer matrix nx by nr ]
%           The stoichiometry matrix. Maps how each elementary reaction
%           changes the number of each state species.
%       .D1 [ real matrix nr by nx]
%           Map of how each state contributes to the reaction rates
%       .D2 [ real matrix nr by nx*nx ]
%           Map of how each product of states in x kron x/vx contributes to
%           the reaction rates
%       .D3 [ real matrix nr by nu*nx ]
%           Map of how each product in u kron x/vx  contributes to
%           the reaction rates
%       .D4 [ real matrix nr by nx*nu ]
%           Map of how each product in x kron u/vu contributes to
%           the reaction rates
%       .D5 [ real matrix nr by nu*nu ]
%           Map of how each product of inputs in u kron u/vu contributes to
%           the reaction rates
%       .D6 [ real matrix nr by nu ]
%           Map of how each input contributes to the reaction rates
%       .d [ nonnegative vector nr]
%           Contribution of constituitive synthesis to the reaction rates
%       .dD1dk [ real matrix nr*nx by nk] (rx_k)
%           Derivative of D1 wrt k
%       .dD2dk [ real matrix nr*nx*nx by nk ] (rxx_k)
%           Derivative of D2 wrt k
%       .dD3dk [ real matrix nr*nu*nx by nk ] (rux_k)
%           Derivative of D3 wrt k
%       .dD4dk [ real matrix nr*nx*nu by nk ] (rxu_k)
%           Derivative of D4 wrt k
%       .dD5dk [ real matrix nr*nu*nu by nk ] (ruu_k)
%           Derivative of D5 wrt k
%       .dD6dk [ real matrix nr*nu by nk ] (ru_k)
%           Derivative of D6 wrt k
%       .dddk [ real matrix nr by nk ] (r_k)
%           Derivative of d wrt k
%       .dD1dk_rk_x [ real matrix nr*nk by nx ] (rk_x)
%           Reshape of dD1dk for faster computation of drdk
%       .dD2dk_rk_xx [ real matrix nr*nk by nx*nx ] (rk_xx)
%           Reshape of dD2dk for faster computation of drdk
%       .dD3dk_rk_ux [ real matrix nr*nk by nu*nx ] (rk_ux)
%           Reshape of dD3dk for faster computation of drdk
%       .dD4dk_rk_xu [ real matrix nr*nk by nx*nu ] (rk_xu)
%           Reshape of dD4dk for faster computation of drdk
%       .dD5dk_rk_uu [ real matrix nr*nk by nu*nu ] (rk_uu)
%           Reshape of dD5dk for faster computation of drdk
%       .dD6dk_rk_u [ real matrix nr*nk by nu ] (rk_u)
%           Reshape of dD6dk for faster computation of drdk
%       .r [ handle @(t,x,u) returns real vector nr ]
%           The rate of each elementary reaction according to the curent
%           state of the system.
%       .drdx [ handle @(t,x,u) returns real matrix nr by nx ]
%           The partial derivative of r wrt x
%       .drdu [ handle @(t,x,u) returns real matrix nr by nu ]
%           The partial derivative of r wrt u
%       .drdk [ handle @(t,x,u) returns real matrix nr by nk ]
%           The partial derivative of r wrt k
%       .d2rdx2 [ handle @(t,x,u) returns real matrix nr*nx by nx ] (rx_x)
%           The partial derivative of drdx wrt x
%       .d2rdk2 [ handle @(t,x,u) returns real matrix nr*nk by nk ] (rk_k)
%           The partial derivative of drdk wrt k
%       .d2rdxdk [ handle @(t,x,u) returns real matrix nr*nk by nx ] (rk_x)
%           The partial derivative of drdk wrt x
%       .d2rdkdx [ handle @(t,x,u) returns real matrix nr*nx by nk ] (rx_k)
%           The partial derivative of drdx wrt k. This is a reshape of
%           d2rdxdk.
%       .v [ handle @(t,x,u) returns positive vector nv ]
%           The size of each compartment
%       .dvdx [ handle @(t,x,u) returns positive matrix nv by nx ]
%           Partial derivative of v wrt x
%       .dvdu [ handle @(t,x,u) returns positive matrix nv by nu ]
%           Partial derivative of v wrt u
%       .d2vdx2 [ handle @(t,x,u) returns positive matrix nv*nx by nx ]
%           Partial derivative of dvdx wrt x
%       .d2vdu2 [ handle @(t,x,u) returns positive matrix nv*nu by nu ]
%           Partial derivative of dvdu wrt u
%       .d2vdudx [ handle @(t,x,u) returns positive matrix nv*nx by nu ]
%           Partial derivative of dvdx wrt u
%       .d2vdxdu [ handle @(t,x,u) returns positive matrix nv*nu by nx ]
%           Partial derivative of dvdu wrt x
%       .Ready [ boolean scalar ]
%           Indicates whether all the mathematical properties of the model
%           (x0, f, S, ...) are up-to-date with its structural properties
%           (Compartments, Reactions, ...). Anytime the model is modified
%           by one of the Kronecker model-building functions, this field is
%           set to false. The model may contain inconsistent or out-of-date
%           information until FinalizeModel is called, which refreshes the
%           structure of the model.
%       .add [ struct scalar ]
%           Functions that add components to the model put them here until
%           FinalizeModel is called. Storing them here and processing them
%           all at once at the end is faster than lengthing each struct
%           vector one element a time.
%       .Update [ handle @(k,x0,q) return struct scalar ]
%           This function handle allows the parameter values of the model
%           to be changed without having to rebuild the model. Each vector
%           k, x0, and q must have the same size as their model
%           counterparts. The structure returned is a mass action Kronecker
%           model identical to the one to which Update is attached except
%           that the parameter values have been changed and the appropriate
%           matrices and function handles have been updated.

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if nargin < 1
    name = [];
end

m.Type = 'Model.MassActionAmount';
m.Name = fixModelName(name);

m.Compartments = emptyCompartments(0);
m.Parameters   = emptyParameters(0);
m.Seeds        = emptySeeds(0);
m.Inputs       = emptyInputs(0);
m.States       = emptyStates(0);
m.Reactions    = emptyReactions(0);
m.Outputs      = emptyOutputs(0);
m.Rules        = growRulesAnalytic;

m = initializeModelBase(m);

m.add.Compartments = growCompartments([], 0);
m.add.Parameters   = growParameters([], 0);
m.add.Seeds        = growSeeds([], 0);
m.add.Inputs       = growInputs([], 0);
m.add.States       = growStates([], 0);
m.add.Reactions    = growReactions([], 0);
m.add.Outputs      = growOutputs([], 0);
m.add.Rules        = growRulesAnalytic;
