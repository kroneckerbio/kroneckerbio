function sol = integrateObj(m, con, obj, opts)

% Constants
nx = m.nx;
nObj = numel(obj);

% Construct system
[der, jac, del] = constructSystem();

if ~con.SteadyState
    order = 0;
    x0 = extractICs(m,con,opts,order);
    ic = [x0; 0];
else
    x0 = steadystateSys(m, con, opts);
    ic = [x0; 0];
end

% Integrate [f; g] over time
sol = accumulateOdeFwd(der, jac, 0, con.tF, ic, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx+1), del);
sol.u = con.u;
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;
sol.k = m.k;
sol.s = con.s;
sol.q = con.q;
sol.h = con.h;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating f and g %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [der, jac, del] = constructSystem()
        f     = m.f;
        dfdx  = m.dfdx;
        uf    = con.u;
        d     = con.d;
        x0    = m.x0;
        nd    = m.ns;
        
        der = @derivative;
        jac = @jacobian;
        del = @delta;
        
        % Derivative of [x; G] with respect to time
        function val = derivative(t, joint)
            u = uf(t);            
            x = joint(1:nx);
            
            % Sum continuous objective functions
            g = 0;
            for i = 1:nObj
                g = g + opts.ObjWeights(i) * obj(i).g(t, x, u);
            end
            
            val = [f(t, x, u); g];
        end
        
        % Jacobian of [x; G] derivative
        function val = jacobian(t, joint)
            u = uf(t);
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
            val = [x0(d(t)) - x0(zeros(nd,1)); 0];
        end
    end
end
