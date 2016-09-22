function ExportModelAnalyticSBML(m, filename)
% Export kroneckerbio analytic model to SBML. Uses a SimBiology model as an
%   intermediate, relying on SimBiology's SBML exporter. See the
%   documentation for `analytic2simbio` for features and limitations.
% See https://www.mathworks.com/help/simbio/ref/sbmlexport.html

simbio = analytic2simbio(m);

sbmlexport(simbio, filename);