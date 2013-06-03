function [yinterp,ypinterp] = ntrp15s(tinterp,t,y,tnew,ynew,h,dif,k,idxNonNegative)
%NTRP15S  Interpolation helper function for ODE15S.
%   YINTERP = NTRP15S(TINTERP,T,Y,TNEW,YNEW,H,DIF,K,IDX) uses data computed in
%   ODE15S to approximate the solution at time TINTERP. TINTREP may be a
%   scalar or a row vector.   
%   [YINTERP,YPINTERP] = NTRP15S(TINTERP,T,Y,TNEW,YNEW,H,DIF,K,IDX) returns
%   also the derivative of the polynomial approximating the solution. 
%   
%   IDX has indices of solution components that must be non-negative. Negative 
%   YINTERP(IDX) are replaced with zeros and the derivative YPINTERP(IDX) is 
%   set to zero.
%   
%   See also ODE15S, DEVAL.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.12.4.6 $  $Date: 2005/04/18 22:11:49 $

s = (tinterp - tnew)/h;     

ypinterp = [];
if k == 1
  yinterp = ynew(:,ones(size(tinterp))) + dif(:,1) * s;
  if nargout > 1
    hdif = (1/h)*dif(:,1);
    ypinterp = hdif(:,ones(size(tinterp)));
  end  
else                    % cumprod collapses vectors
  K = (1:k)';
  kI = K(:,ones(size(tinterp)));  
  yinterp = ynew(:,ones(size(tinterp))) + ...
      dif(:,K) * cumprod((s(ones(k,1),:)+kI-1)./kI);
  
  if nargout > 1
    ypinterp = dif(:,ones(size(tinterp)));
    S  = ones(size(tinterp));
    dS = ones(size(tinterp)); 
    for i=2:k
      S = S .* (i-2+s)/i;
      dS = dS .* (i-1+s)/i + S;
      ypinterp = ypinterp + dif(:,i)*dS;    
    end
    ypinterp = ypinterp/h;
  end  
  
end

% Non-negative solution
if ~isempty(idxNonNegative)
  idx = find(yinterp(idxNonNegative,:)<0); % vectorized
  if ~isempty(idx)
    w = yinterp(idxNonNegative,:);
    w(idx) = 0;
    yinterp(idxNonNegative,:) = w;
    if nargout > 1   % the derivative
      w = ypinterp(idxNonNegative,:);
      w(idx) = 0;
      ypinterp(idxNonNegative,:) = w;
    end      
  end
end  
