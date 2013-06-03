function [odeFcn,thresholdNonNegative] = odenonnegative(ode,y0,threshold,idxNonNegative)  
%ODENONNEGATIVE  Helper function for handling nonnegative solution constraints
%   Modify the derivative function to prevent the solution from crossing zero.
%
%   See also ODE113, ODE15S, ODE23, ODE23T, ODE23TB, ODE45.

%   Jacek Kierzenka
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/05/31 16:31:08 $

neq = numel(y0);
thresholdNonNegative = [];
if any( (idxNonNegative < 1) | (idxNonNegative > neq) )
  error('MATLAB:odenonnegative:NonNegativeIndicesInvalid',...
        ['All indices specified in ODESET(''NonNegative'',idx) '...
         'must satisfy 1 <= idx <= neq']);
end
if any(y0(idxNonNegative) < 0)
  error('MATLAB:odenonnegative:NonNegativeViolatedAtT0',...
        'Non-negativity constraint violated at t0.');
end  
if length(threshold) == 1
  thresholdNonNegative = threshold(ones(size(idxNonNegative)));
else
  thresholdNonNegative = threshold(idxNonNegative);
end
thresholdNonNegative = thresholdNonNegative(:);
odeFcn = @local_odeFcn_nonnegative;   

% -----------------------------------------------------------
% Nested function: ODE with nonnegativity constraints imposed
%
  function yp = local_odeFcn_nonnegative(t,y,varargin)
    ndx = idxNonNegative(y(idxNonNegative) <= 0);
    y(ndx) = max(y(ndx),0);
    yp = feval(ode,t,y,varargin{:}); 
    yp(ndx) = max(yp(ndx),0);
  end  % local_odeFcn_nonnegative
% -----------------------------------------------------------

end  % odenonnegative

