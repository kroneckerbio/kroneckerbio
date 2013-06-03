function sol = integrateObj(m, con, obj, opts)
% m is a scalar
% con is a scalar
% obj is a vector

% Constants
nx = m.nx;
nObj = numel(obj);

% Construct system
[der, jac, events] = constructSystem();

%TODO: events

% Initial conditions
if ~con.SteadyState
    if opts.UseModelICs
        ic = [m.x0; 0];
    else
        ic = [con.x0; 0];
    end
else
    ic = [steadystateSys(m, con, opts); 0];
end

% Input
if opts.UseModelInputs
    u = m.u;
else
    u = con.u;
end

% Integrate
sol = accumulateOde(der, jac, 0, con.tF, ic, u, con.Discontinuities, 1:nx, opts.RelTol, opts.AbsTol(1:nx+1));
sol.u = u;
sol.C1 = m.C1;
sol.C2 = m.C2;
sol.c  = m.c;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating x and g %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [gDer, gJac, gEvents] = constructSystem()
        
        f      = m.f;
        dfdx   = m.dfdx;
        
        gDer = @derivative;
        gJac = @jacobian;
        if false%~isempty(obj.Events)
            %TODO: events
            events  = obj.Events;
            gEvents = @eventFun;
        else
            gEvents = [];
        end
        
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
        
        % Event function for system
        function [val, dir, term] = eventFun(t, joint, u)
            u = u(t);
            x = joint(1:nx);
            [val, dir, term] = events(t, x, u);
        end
    end
end