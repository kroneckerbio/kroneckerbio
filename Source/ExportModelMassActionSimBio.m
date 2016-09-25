function simbio = ExportModelMassActionSimBio(m, opts)
% See the documentation for `massaction2simbio`
if nargin < 2
    opts = [];
end

simbio = massaction2simbio(m, opts);