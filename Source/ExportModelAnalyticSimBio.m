function simbio = ExportModelAnalyticSimBio(m, opts)
% See the documentation for `analytic2simbio`
if nargin < 2
    opts = [];
end

simbio = analytic2simbio(m, opts);