function ExportModelMassActionSBML(m, filename)
% Export kroneckerbio massaction model to SBML. Uses a SimBiology model as an
%   intermediate, relying on SimBiology's SBML exporter. See the
%   documentation for `massaction2simbio` for features and limitations.
% See https://www.mathworks.com/help/simbio/ref/sbmlexport.html
%
% Inputs:
%   m [ Model.MassActionAmount struct ]
%       Kroneckerbio analytic model
%   filename [ string ]
%       Name of output SBML file

simbio = massaction2simbio(m);

sbmlexport(simbio, filename);