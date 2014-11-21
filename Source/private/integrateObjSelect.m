function sol = integrateObjSelect(m, con, obj, tGet, opts)

% Constants
nx = m.nx;
nObj = numel(obj);

% Construct system
[der, jac, del] = constructSystem();

% Initial conditions
if opts.UseModelSeeds
    s = m.s;
else
    s = con.s;
end

if ~con.SteadyState
    x0 = m.dx0ds * s + m.x0c;
    ic = [x0; 0];
else
    ic = [steadystateSys(m, con, opts); 0];
end

% Input
if opts.UseModelInputs
    u = m.u;
    q = m.q;
else
    u = con.u;
    q = con.q;
end

% Integrate
sol = accumulateOdeFwdSelect(der, jac, 0, con.tF, ic, u, con.Discontinuities, tGet, 1:nx, opts.RelTol, opts.AbsTol(1:nx+1)), del);
sol.u = u(tGet);
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;
sol.k = m.k;
sol.s = s;
sol.q = q;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating x and g %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
        f     = m.f;
        dfdx  = m.dfdx;
        d     = con.d;
        dx0ds = m.dx0ds;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [x; G] with respect to time
        function val = derivative(t, joint, u)
            u = u(t);            
            x = joint(1:nx);
            
            % Sum continuous objective functions
            g = 0;
            for i = 1:nObj
                g = g + opts.ObjWeights(i) * obj(i).g(t, x, u);
            end
            
            val = [f(t, x, u); g];
        end
        
        % Jacobian of [x; G] derivative
        function val = jacobian(t, joint, u)
            u = u(t);
            x = joint(1:nx);
            
            % Sum continuous objective gradients
            dgdx = zeros(1,nx);
            for i = 1:nObj
                dgdx = dgdx + opts.ObjWeights(i) * vec(obj(i).dgdx(t, x, u)).';
            end
            
            val = [dfdx(t, x, u), sparse(nx,1);
                            dgdx,            0];
        end

        % Dosing
        function val = delta(t, joint)
            val = [dx0ds * d(t); 0];
        end
    end
end