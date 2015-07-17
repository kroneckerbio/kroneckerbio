function valid = isValidSymbolicModel(s, detailed)
% Verify that struct is a valid symbolic model. Contains the most up-to-date
% spec for Model.Symbolic.
% Inputs:
%   s [ struct ]
%       Test symbolic model
%   detailed [ true | {false} ]
%       Whether to perform detailed analysis on valid names/ids
% Outputs:
%   valid [ true | false ]
%       Whether s is a valid symbolic model
%
% Note: Throws error when model is invalid for now. Consider fixing to errors
%   are caught returning valid = false (w/ helpful message)
% TODO: implement detailed checking
% TODO: if this function remains fast, put in all model conversion code as
%   a sanity check

if nargin < 2
    detailed = false;
end

%% Correct type signature
assert(strcmp(s.Type, 'Model.Symbolic'),'isValidSymbolicModel:TypeError', 'Model has invalid type %s', s.Type)

%% Correct fields
expectedFields = {'Type', 'Name', ...
    'nv', 'vNames', 'vIDs', 'v', 'dv', ...
    'nx', 'xNames', 'xIDs', 'xvNames', 'x', ...
    'nu', 'uNames', 'uIDs', 'uvNames', 'u', ...
    'ns', 'sNames', 'sIDs', 's', ...
    'nk', 'kNames', 'kIDs', 'k', ...
    'nr', 'rNames', 'rIDs', 'r', ...
    'nz', 'zNames', 'zIDs', 'z'};
fields = fieldnames(s);
missingFields = setdiff(expectedFields, fields); % (arg1 not in arg2)
extraFields = setdiff(fields, expectedFields);
assert(isempty(missingFields), 'isValidSymbolicModel:MissingFields', 'Model struct missing fields: %s', cellstr2str(missingFields))
if ~isempty(extraFields)
    warning('isValidSymbolicModel:ExtraFields', 'Model struct has extra fields: %s', cellstr2str(extraFields))
end

%% Compartments
assert(isscalar(s.nv)   && isnumeric(s.nv) && s.nv >= 1 && mod(s.nv,1) == 0, 'isValidSymbolicModel:nvType', 'nv must be a scalar positive integer >= 1')
assert(iscell(s.vNames) && all(size(s.vNames) == [s.nv, 1]), 'isValidSymbolicModel:vNamesType', 'vNames must be a nv x 1 cell vector')
assert(iscell(s.vIDs)   && all(size(s.vIDs)   == [s.nv, 1]), 'isValidSymbolicModel:vIDsType', 'vIDs must be a nv x 1 cell vector')
assert(isnumeric(s.v)   && all(size(s.v)      == [s.nv, 1]), 'isValidSymbolicModel:vType', 'v must be a nv x 1 double vector')
assert(isnumeric(s.dv)  && all(size(s.dv)     == [s.nv, 1]), 'isValidSymbolicModel:dvType', 'dv must be a nv x 1 double vector')

assert(all(cellfun(@ischar, s.vNames)), 'isValidSymbolicModel:vNamesCellTypes', 'cells in vNames must be strings')
assert(all(cellfun(@ischar, s.vIDs)), 'isValidSymbolicModel:vIDsCellTypes', 'cells in vIDs must be strings')

%% States
assert(isscalar(s.nx)    && isnumeric(s.nx) && s.nx >= 0 && mod(s.nx,1) == 0, 'isValidSymbolicModel:nxType', 'nx must be a scalar positive integer')
assert(iscell(s.xNames)  && all(size(s.xNames)  == [s.nx, 1]), 'isValidSymbolicModel:xNamesType', 'xNames must be a nx x 1 cell vector')
assert(iscell(s.xIDs)    && all(size(s.xIDs)    == [s.nx, 1]), 'isValidSymbolicModel:xIDsType', 'xIDs must be a nx x 1 cell vector')
assert(iscell(s.xvNames) && all(size(s.xvNames) == [s.nx, 1]), 'isValidSymbolicModel:xvNamesType', 'xvNames must be a nx x 1 cell vector')
assert(iscell(s.x)       && all(size(s.x)       == [s.nx, 1]), 'isValidSymbolicModel:xType', 'x must be a nx x 1 cell vector')

assert(all(cellfun(@ischar, s.xNames)), 'isValidSymbolicModel:xNamesCellTypes', 'cells in xNames must be strings')
assert(all(cellfun(@ischar, s.xIDs)), 'isValidSymbolicModel:xIDsCellTypes', 'cells in xIDs must be strings')
assert(all(cellfun(@ischar, s.xvNames)), 'isValidSymbolicModel:xvNamesCellTypes', 'cells in xvNames must be strings')
assert(all(cellfun(@ischar, s.x)), 'isValidSymbolicModel:xCellTypes', 'cells in x must be strings')

%% Inputs
assert(isscalar(s.nu) && isnumeric(s.nu) && s.nu >= 0 && mod(s.nu,1) == 0, 'isValidSymbolicModel:nuType', 'nu must be a scalar positive integer')
assert(iscell(s.uNames)  && all(size(s.uNames)  == [s.nu, 1]), 'isValidSymbolicModel:uNamesType', 'uNames must be a nu x 1 cell vector')
assert(iscell(s.uIDs)    && all(size(s.uIDs)    == [s.nu, 1]), 'isValidSymbolicModel:uIDsType', 'uIDs must be a nu x 1 cell vector')
assert(iscell(s.uvNames) && all(size(s.uvNames) == [s.nu, 1]), 'isValidSymbolicModel:uvNamesType', 'uvNames must be a nu x 1 cell vector')
assert(isnumeric(s.u)    && all(size(s.u)       == [s.nu, 1]), 'isValidSymbolicModel:uType', 'u must be a nu x 1 double vector')

assert(all(cellfun(@ischar, s.uNames)), 'isValidSymbolicModel:uNamesCellTypes', 'cells in uNames must be strings')
assert(all(cellfun(@ischar, s.uIDs)), 'isValidSymbolicModel:uIDsCellTypes', 'cells in uIDs must be strings')
assert(all(cellfun(@ischar, s.uvNames)), 'isValidSymbolicModel:uvNamesCellTypes', 'cells in uvNames must be strings')

%% Seeds
assert(isscalar(s.ns) && isnumeric(s.ns) && s.ns >= 0 && mod(s.ns,1) == 0, 'isValidSymbolicModel:nsType', 'ns must be a scalar positive integer')
assert(iscell(s.sNames) && all(size(s.sNames) == [s.ns, 1]), 'isValidSymbolicModel:sNamesType', 'sNames must be a ns x 1 cell vector')
assert(iscell(s.sIDs)   && all(size(s.sIDs)   == [s.ns, 1]), 'isValidSymbolicModel:sIDsType', 'sIDs must be a ns x 1 cell vector')
assert(isnumeric(s.s)   && all(size(s.s)      == [s.ns, 1]), 'isValidSymbolicModel:sType', 's must be a ns x 1 double vector')

assert(all(cellfun(@ischar, s.sNames)), 'isValidSymbolicModel:sNamesCellTypes', 'cells in sNames must be strings')
assert(all(cellfun(@ischar, s.sIDs)), 'isValidSymbolicModel:sIDsCellTypes', 'cells in sIDs must be strings')

%% Parameters
assert(isscalar(s.nk)   && isnumeric(s.nk) && s.nk >= 0 && mod(s.nk,1) == 0, 'isValidSymbolicModel:nkType', 'nk must be a scalar positive integer')
assert(iscell(s.kNames) && all(size(s.kNames) == [s.nk, 1]), 'isValidSymbolicModel:kNamesType', 'kNames must be a nk x 1 cell vector')
assert(iscell(s.kIDs)   && all(size(s.kIDs)   == [s.nk, 1]), 'isValidSymbolicModel:kIDsType', 'kIDs must be a nk x 1 cell vector')
assert(isnumeric(s.k)   && all(size(s.k)      == [s.nk, 1]), 'isValidSymbolicModel:kType', 'k must be a nk x 1 double vector')

assert(all(cellfun(@ischar, s.kNames)), 'isValidSymbolicModel:kNamesCellTypes', 'cells in kNames must be strings')
assert(all(cellfun(@ischar, s.kIDs)), 'isValidSymbolicModel:kIDsCellTypes', 'cells in kIDs must be strings')

%% Reactions
assert(isscalar(s.nr) && isnumeric(s.nr) && s.nr >= 0 && mod(s.nr,1) == 0, 'isValidSymbolicModel:nrType', 'nr must be a scalar positive integer')
assert(iscell(s.rNames) && all(size(s.rNames) == [s.nr, 1]), 'isValidSymbolicModel:rNamesType', 'rNames must be a nr x 1 cell vector')
assert(iscell(s.rIDs)   && all(size(s.rNames) == [s.nr, 1]), 'isValidSymbolicModel:rIDsType', 'rIDs must be a nr x 1 cell vector')
assert(iscell(s.r)      && all(size(s.r)      == [s.nr, 3]), 'isValidSymbolicModel:rType', 'r must be a nr x 3 cell matrix')

assert(all(cellfun(@ischar, s.rNames)), 'isValidSymbolicModel:rNamesCellTypes', 'cells in rNames must be strings')
assert(all(cellfun(@ischar, s.rIDs)), 'isValidSymbolicModel:rIDsCellTypes', 'cells in rIDs must be strings')
% Reactants and products must be 1xn cell arrays, including a 1x0 empty cell array
assert(all(cellfun(@iscell, s.r(:,1))), 'isValidSymbolicModel:rReactantTypes', 'Reaction reactants must be in cell arrays')
assert(all(cellfun(@iscell, s.r(:,2))), 'isValidSymbolicModel:rProductTypes', 'Reaction products must be in cell arrays')
assert(all(cellfun(@ischar, s.r(:,3))), 'isValidSymbolicModel:rExprTypes', 'Reaction rate expressions must be strings')

%% Rules
assert(isscalar(s.nz)    && isnumeric(s.nz) && s.nz >= 0 && mod(s.nz,1) == 0, 'isValidSymbolicModel:nzType', 'nz must be a scalar positive integer')
assert(iscell(s.zNames)  && all(size(s.zNames) == [s.nz, 1]), 'isValidSymbolicModel:zNamesType', 'zNames must be a nz x 1 cell vector')
assert(iscell(s.zIDs)    && all(size(s.zIDs)   == [s.nz, 1]), 'isValidSymbolicModel:zIDsType', 'zIDs must be a nz x 1 cell vector')
assert(iscell(s.z)       && all(size(s.z)      == [s.nz, 3]), 'isValidSymbolicModel:zType', 'z must be a nz x 3 cell matrix')

assert(all(cellfun(@ischar, s.zNames)), 'isValidSymbolicModel:zNamesCellTypes', 'cells in zNames must be strings')
assert(all(cellfun(@ischar, s.zIDs)), 'isValidSymbolicModel:zIDsCellTypes', 'cells in zIDs must be strings')
assert(all(cellfun(@ischar, s.z(:,1))), 'isValidSymbolicModel:RuleTargetType', 'Rule targets must be strings')
assert(all(cellfun(@ischar, s.z(:,2))), 'isValidSymbolicModel:RuleExprType', 'Rule expressions must be strings')
assert(all(cellfun(@isValidRuleType, s.z(:,3))), 'isValidSymbolicModel:RuleType', 'Rules must be repeated assignment or initial assignment')

%% Detailed checks
% TODO: implement these
%   - Names aren't blank strings
%   - IDs are valid variables (+ possibly not substrings of each other - required for symbolic subs later)
%   - Reactants and products in reactions refer to existing things

%% Completed all checks
valid = true;

%% Helper functions
function valid = isValidRuleType(ruleType)
valid = true;
if ~ischar(ruleType)
    valid = false;
    return
end
switch ruleType
    case 'repeated assignment'
        % pass
    case 'initial assignment'
        % pass
    otherwise
        valid = false;
end