function sol = integrateMfk(m, con, opts)

% Constants
nx = m.nx;
nu = m.nu; 
omega = m.Compartments(1).Values;
S = m.S;

A1=m.A1; A2=m.A2; A3=m.A3; A4=m.A4; A5=m.A5; A6=m.A6;
D1=m.D1; D2=m.D2; D3=m.D3; D4=m.D4; D5=m.D5; D6=m.D6;

T_1_nx  = sparse(Tmatrix([1 nx]));
T_nx_nx = sparse(Tmatrix([nx nx]));

aux1 = A3*(kron(T_1_nx,speye(nu)));
aux2 = A4*(kron(T_1_nx,speye(nu)));

Inx=speye(nx);
Inu=speye(nu);

% Construct system
[der, jac] = constructSystem();


% Initial conditions
if opts.UseModelICs
    ic = m.x0;
else
    ic = con.x0;
end


% Append the initial conditions of Covariance Matrix
ic = cat(1,ic,opts.V0(:));


% Input
if opts.UseModelInputs
    u = m.u;
else
    u = con.u;
end

% Integrate x over time
sol = accumulateOde(der, jac, 0, con.tF, ic, u, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol);
sol.u = u;
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating x and dxdT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac] = constructSystem()
        
%        Ip = speye(nT);

        VStart = nx+1;
        VEnd   = nx+nx*nx;
        
        f       = m.f;
        dfdx    = m.dfdx;
        dfdT    = @dfdTSub;
        d2fdx2  = m.d2fdx2;
        d2fdTdx = @d2fdTdxSub;
        
        der = @derivative;
        jac = @jacobian;
        
        % Derivative of [x; dxdT] with respect to time
        function val = derivative(t, joint, u)
            u = u(t);            
            x = joint(1:nx);
            V = reshape(joint(VStart:VEnd), nx,nx);
            
            % and devide by volume to turn into concentrations
            x=x/omega;
            u=u/omega;
            V=V/(omega^2);
            
            
            
            dXdt = f(t, x*omega, u)/omega + A2*V(:);
            
            M =   A1    +    A2 * ( kron(speye(nx), x) + kron(x,speye(nx))) ...
                +    A3*kron(T_1_nx,speye(nu))  * kron(speye(nx),u) ...
                +    A4*kron(speye(nu),T_1_nx) * kron(u,speye(nx)) ;
            
%           M =     A1  +      A2 *     ( kron(Inx, x) + kron(x,Inx))     +                 aux1*(kron(Inx,u))         +                     aux2*kron(u,Inu) ;          
            
            
            Lamb = diag( D1*x + D2*kron(x,x) + D3*kron(u,x) + D4*kron(x,u) + D5*kron(u,u) + D6*u ) ;
                       
            
            dVdt = M*V + V*M' + (1/omega)*(S*Lamb*S');
            
            % Combine and change back to amounts and not concetrations
            val = [dXdt*(omega); dVdt(:)*(omega^2)];
            
            
           
        end
        
        % Jacobian of [X; V] derivative
        function val = jacobian(t, joint, u)
            u = u(t);            
            x = joint(1:nx); % x_
            V = reshape(joint(VStart:VEnd), nx,nx);
            
            % and devide by volume to turn into concentrations
             x=x/omega;
             u=u/omega;
             V=V/(omega^2);
             
             
            M =  A1 + A2 * ( kron(speye(nx), x) + kron(x,speye(nx))) +  A3*(kron(T_1_nx,speye(nu))*(kron(speye(nx),u)))  + A4*(kron(T_1_nx,speye(nu))*kron(u,eye(nx))) ;
           

            dMdx_dA2_2D = kron(vec(speye(nx)), speye(nx)) + kron(T_nx_nx, speye(nx)) * (kron(speye(nx),vec(speye(nx))));
            
            
            dMdx_2D = sparse(0,0);
                for i=1:nx
                   dMdx_2D = cat(1,dMdx_2D,A2*dMdx_dA2_2D( ((i-1)*nx^2+1):((i)*nx^2) , :));
                end
                %check
                
            dMVdx_2D = sparse(0,0);
                for i=1:nx
                   dMVdx_2D = cat(1,dMVdx_2D,dMdx_2D( ((i-1)*nx+1):((i)*nx) , :)*V);
                end
                % inconc
                
            dVMtdx_2D =sparse(0,0);
                for i=1:nx
                   dVMtdx_2D = cat(1,dVMtdx_2D,V*dMdx_2D( ((i-1)*nx+1):((i)*nx) , :));
                end
                % inconc
                
            dvRdx = D1 + D2 * ( kron(speye(nx), x) + kron(x,speye(nx))) +  D3*(kron(T_1_nx,speye(nu))*kron(u,speye(nx)))   + D4*(kron(T_1_nx,speye(nu))*(kron(speye(nx),u))) ;
            
            
            dSRStdx_2D =sparse(0,0);
                for i=1:nx
                   dSRStdx_2D = cat(1,dSRStdx_2D,S*diag(dvRdx(:,i))*S');
                end
            % check
           
            ddvdtdx_2D =  dMVdx_2D + dVMtdx_2D + (1/omega)*dSRStdx_2D;
            ddvdtdx_2D =  reshape(ddvdtdx_2D', [nx*nx nx]);
            % check    
                
              
             % d(dVdt)/dv
             
             dVdv_dM_2D = [];
             for i=1:nx^2
                 aux=sparse(nx,nx);
                 aux(i) = 1;
                 dVdv_dM_2D = cat(1, dVdv_dM_2D, aux);
             end
             %check
             dMVdv_2D = sparse(0,0);
             for i=1:nx^2
                 aux=sparse(nx,nx);
                 aux(i) = 1;
                 
                 dMVdv_2D = cat(2, dMVdv_2D, M*dVdv_dM_2D( ((i-1)*nx+1):((i)*nx) , :)    );
             end
             %check
             
             dVMtdv_2D = sparse(0,0);
             for i=1:nx^2
                 aux=sparse(nx,nx);
                 aux(i) = 1;
                 
                 dVMtdv_2D = cat(2, dVMtdv_2D, dVdv_dM_2D(((i-1)*nx+1):((i)*nx) , :)      *M');
             end
             %check
             
            ddVdtdv_2D = dMVdv_2D + dVMtdv_2D;
            ddVdtdv_2D = reshape(ddVdtdv_2D, [nx*nx nx*nx]);
             
            % Combine
            val = [m.dfdx(t, x*omega, u),  A2*speye(nx^2)/omega;
                       ddvdtdx_2D*(omega), ddVdtdv_2D];
       
        end
        
        
        
    end

end