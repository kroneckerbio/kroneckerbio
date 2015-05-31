function [Sxint,Spxint] = deval(sol,xint,idx)
%DEVAL  Evaluate the solution of a differential equation problem.
%   SXINT = DEVAL(SOL,XINT) evaluates the solution of a differential equation 
%   problem at all the entries of the vector XINT. SOL is a structure returned 
%   by an initial value problem solver (ODE45, ODE23, ODE113, ODE15S, ODE23S, 
%   ODE23T, ODE23TB, ODE15I), a boundary value problem solver (BVP4C, BVP5C), 
%   or a delay differential equations solver (DDE23, DDESD). The elements of 
%   XINT must be in the interval [SOL.x(1) SOL.x(end)]. For each I, SXINT(:,I) 
%   is the solution corresponding to XINT(I). 
%
%   SXINT = DEVAL(SOL,XINT,IDX) evaluates as above but returns only
%   the solution components with indices listed in IDX.   
%
%   SXINT = DEVAL(XINT,SOL) and SXINT = DEVAL(XINT,SOL,IDX) are also acceptable.
%
%   [SXINT,SPXINT] = DEVAL(...) evaluates as above but returns also the value 
%   of the first derivative of the polynomial interpolating the solution.
%
%   For multipoint boundary value problems or initial value problems extended 
%   using ODEXTEND, the solution might be discontinuous at interfaces. 
%   For an interface point XC, DEVAL returns the average of the limits from 
%   the left and right of XC. To get the limit values, set the XINT argument 
%   of DEVAL to be slightly smaller or slightly larger than XC.
%
%   Class support for inputs SOL and XINT:
%     float: double, single
%
%   See also ODE45, ODE23, ODE113, ODE15S, ODE23S, ODE23T, ODE23TB, ODE15I, 
%            DDE23, DDESD, DDENSD, BVP4C, BVP5C.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2012 The MathWorks, Inc.

if ~isa(sol,'struct')  
  % Try  DEVAL(XINT,SOL)
  temp = sol;
  sol  = xint;
  xint = temp;
end

try
  t = sol.x;    
  y = sol.y;    
catch
  error(message('MATLAB:deval:SolNotFromDiffEqSolver', inputname( 1 )));
end

if nargin < 3
  idx = 1:size(y,1);  % return all solution components
else 
  if any(idx < 0) || any(idx > size(y,1))
    error(message('MATLAB:deval:IDXInvalidSolComp', inputname( 3 )));
  end  
end  
idx = idx(:);
  
if isfield(sol,'solver')
  solver = sol.solver;
else
  if isfield(sol,'yp')
    warning(message('MATLAB:deval:MissingSolverField', inputname(1), inputname(1)));  
    solver = 'bvp4c';
  else
    error(message('MATLAB:deval:NoSolverInStruct',inputname(1)));
  end
end

% Select appropriate interpolating function.
switch solver
 case 'ode113'
  interpfcn = @ntrp113;
 case 'ode15i'
  interpfcn = @ntrp15i;
 case 'ode15s'
  interpfcn = @ntrp15s;
 case 'ode23'
  interpfcn = @ntrp23;
 case 'ode23s'
  interpfcn = @ntrp23s;
 case 'ode23t'
  interpfcn = @ntrp23t;
 case 'ode23tb'
  interpfcn = @ntrp23tb;
 case 'ode45'
  interpfcn = @ntrp45;
 case {'bvp4c','dde23','ddesd','ddensd'}
  interpfcn = @ntrp3h;
 case 'bvp5c'
  interpfcn = @ntrp4h;  
 otherwise
  error(message('MATLAB:deval:InvalidSolver', solver, inputname( 1 ))); 
end

% If necessary, convert sol.idata to MATLAB R14 format. 
if ~isfield(sol,'extdata') && ismember(solver,{'ode113','ode15s','ode23','ode45'})   
  sol.idata = convert_idata(solver,sol.idata);  
end

% Determine the dominant data type.
dataType = superiorfloat(sol.x,xint);

% Evaluate the first derivative?
Spxint_requested = (nargout > 1);   

% Non-negativity constraint
if isfield(sol, 'idata') && isfield(sol.idata, 'idxNonNegative')
    idxNonNegative = sol.idata.idxNonNegative;
else
    idxNonNegative = [];
end

% Allocate output arrays.
n = length(idx);
Nxint = length(xint);
Sxint = zeros(n,Nxint,dataType);
if Spxint_requested
  Spxint = zeros(n,Nxint,dataType);
end

% Make tint a row vector and if necessary, 
% sort it to match the order of t.
tint = xint(:).';  
tdir = sign(t(end) - t(1));
had2sort = any(tdir*diff(tint) < 0);
if had2sort
  [tint,tint_order] = sort(tdir*tint);
  tint = tdir*tint;
end  

% Using the sorted version of tint, test for illegal values.
if any(isnan(tint)) || ~isempty(tint) && ((tdir*(tint(1) - t(1)) < 0) || (tdir*(tint(end) - t(end)) > 0))
  error(message('MATLAB:deval:SolOutsideInterval',sprintf('%e',t(1)),sprintf('%e',t(end))));
end

evaluated = 0;
bottom = 1;
while evaluated < Nxint
  
  % Find right-open subinterval [t(bottom), t(bottom+1)) containing the next entry of tint. 
  % Unless t(bottom) == t(end), a proper interval is returned: t(bottom+1) ~= t(bottom).
  Index = find( tdir*(t(bottom:end) - tint(evaluated+1)) <= 0, 1, 'last'); % one-based Index
  bottom = bottom + (Index - 1);  % convert to zero-based Index

  % Is it [t(end), t(end)]?
  at_tend = (t(bottom) == t(end));

  % Return solution already available at t(bottom)
  index1 = find(tint(evaluated+1:end) == t(bottom));
    
  % Interpolate solution inside (t(bottom), t(bottom+1))
  if at_tend
    index2 = [];
  else
    index2 = find( (tdir*(tint(evaluated+1:end) - t(bottom)) > 0) & ...
                   (tdir*(tint(evaluated+1:end) - t(bottom+1)) < 0) );
  end
  
  % Return the (adjusted) solution at t(bottom)
  if ~isempty(index1)
      if at_tend          
          yint1 = y(:,end);
          if Spxint_requested
              % Extrapolate derivative from [t(bottom-1),t(bottom))
              interpdata = extract_idata(solver,sol,t,bottom-1,idxNonNegative);
              [~,ypint1] = interpfcn(t(bottom),t(bottom-1),y(:,bottom-1),...
                                   t(bottom),y(:,bottom),interpdata{:});    
          end
          
      elseif (bottom > 2) && (t(bottom) == t(bottom-1)) % Interface point
          % Average the solution (and its derivative) across the interface.
          yLeft  = y(:,bottom-1);       
          yRight = y(:,bottom);
          yint1 = (yLeft + yRight)/2;          
          if Spxint_requested
              % Get the 'left' derivative by extrapolating in [t(bottom-2), t(bottom-1)).        
              interpdata =  extract_idata(solver,sol,t,bottom-2,idxNonNegative); 
              [~,ypLeft] = interpfcn(t(bottom-1),t(bottom-2),y(:,bottom-2),...
                                    t(bottom-1),y(:,bottom-1),interpdata{:});
              % Get the 'right' derivative by interpolating in [t(bottom), t(bottom+1)).        
              interpdata =  extract_idata(solver,sol,t,bottom,idxNonNegative); 
              [~,ypRight] = interpfcn(t(bottom),t(bottom),y(:,bottom),...
                                    t(bottom+1),y(:,bottom+1),interpdata{:});
              ypint1 = (ypLeft + ypRight)/2;
          end
          warning(message('MATLAB:deval:NonuniqueSolution',sprintf('%g',t(bottom))));
      else  
          % 'Regular' mesh point
          yint1 = y(:,bottom); 
          if Spxint_requested
              % Interpolate derivative from [t(bottom),t(bottom+1))
              interpdata = extract_idata(solver,sol,t,bottom,idxNonNegative);
              [~,ypint1] = interpfcn(t(bottom),t(bottom),y(:,bottom),...
                                   t(bottom+1),y(:,bottom+1),interpdata{:});    
          end                    
      end
      
      % Accumulate the output.
      Sxint(:,evaluated+index1) = yint1(idx,ones(1,numel(index1)));  
      if Spxint_requested
          Spxint(:,evaluated+index1) = ypint1(idx,ones(1,numel(index1)));  
      end  
  end
    
  % Interpolate solution inside (t(bottom), t(bottom+1)).
  if ~isempty(index2) 
      % Get solver-dependent interpolation data for [t(bottom), t(bottom+1)).
      interpdata = extract_idata(solver,sol,t,bottom,idxNonNegative);
      
      % Evaluate the interpolant at all points from (t(bottom), t(bottom+1)).
      if Spxint_requested
          [yint2,ypint2] = interpfcn(tint(evaluated+index2),t(bottom),y(:,bottom),...
                                   t(bottom+1),y(:,bottom+1),interpdata{:});    
      else  
          yint2 = interpfcn(tint(evaluated+index2),t(bottom),y(:,bottom),...
                           t(bottom+1),y(:,bottom+1),interpdata{:});    
      end
  
      % Accumulate the output.
      Sxint(:,evaluated+index2) = yint2(idx,:);  
      if Spxint_requested
          Spxint(:,evaluated+index2) = ypint2(idx,:);  
      end  
  end
    
  evaluated = evaluated + length(index1) + length(index2);      
end

if had2sort     % Restore the order of tint in the output.
  Sxint(:,tint_order) = Sxint;  
  if Spxint_requested
    Spxint(:,tint_order) = Spxint;  
  end  
end

% --------------------------------------------------------------------------

function idataOut = convert_idata(solver,idataIn)
% Covert an old sol.idata to the MATLAB R14 format

idataOut = idataIn;
switch solver
 case 'ode113'
  idataOut.phi3d = shiftdim(idataIn.phi3d,1);
 case 'ode15s'
  idataOut.dif3d = shiftdim(idataIn.dif3d,1);
 case {'ode23','ode45'}
  idataOut.f3d = shiftdim(idataIn.f3d,1);
end

% --------------------------------------------------------------------------

function interpdata = extract_idata(solver,sol,t,tidx,idxNonNegative)
% Data for interpolation in [t(tidx), t(tidx+1))

switch solver
 case 'ode113'
  interpdata = { sol.idata.klastvec(tidx+1), ...
                 sol.idata.phi3d(:,:,tidx+1), ...
                 sol.idata.psi2d(:,tidx+1), ...
                 idxNonNegative };    
  
 case 'ode15i'
  k = sol.idata.kvec(tidx+1);
  interpdata = { sol.x(tidx:-1:tidx-k+1), ...        
                 sol.y(:,tidx:-1:tidx-k+1) };
  
 case 'ode15s'
  interpdata = { t(tidx+1)-t(tidx), ...
                 sol.idata.dif3d(:,:,tidx+1), ...
                 sol.idata.kvec(tidx+1), ...
                 idxNonNegative };
  
 case {'ode23','ode45'} 
  interpdata = { t(tidx+1)-t(tidx), ...
                 sol.idata.f3d(:,:,tidx+1), ...
                 idxNonNegative };
  
 case 'ode23s'          
  interpdata = { t(tidx+1)-t(tidx), ...
                 sol.idata.k1(:,tidx+1), ...
                 sol.idata.k2(:,tidx+1) };
  
 case 'ode23t'
  interpdata = { t(tidx+1)-t(tidx), ...
                 sol.idata.z(:,tidx+1), ...
                 sol.idata.znew(:,tidx+1), ...
                 idxNonNegative };
  
 case 'ode23tb'       
  interpdata = { sol.idata.t2(tidx+1) ...
                 sol.idata.y2(:,tidx+1), ...
                 idxNonNegative };
  
 case {'bvp4c','dde23','ddesd','ddensd'}  
  interpdata = { sol.yp(:,tidx), ...
                 sol.yp(:,tidx+1) };

 case 'bvp5c'  
  interpdata = { sol.idata.ymid(:,tidx), ...
                 sol.idata.yp(:,tidx), ...
                 sol.idata.yp(:,tidx+1) };    
  
end

