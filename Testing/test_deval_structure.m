% Test minimal two-point solution with deval
ic = [1; 2; 3];
t0 = 0;
epsilon = 1e-10;
sol.solver = 'ode15s';
sol.extdata = [];
sol.x = [t0, t0 + epsilon];
sol.y = [ic, ic];
sol.stats = [];
sol.idata.kvec = [0, 0];
sol.idata.dif3d = zeros(size(ic,1), 5, 2);
sol.idata.idxNonNegative = [];

disp('Testing deval with our two-point solution:');
disp(['sol.x = ' mat2str(sol.x)]);
disp(['size(sol.y) = ' mat2str(size(sol.y))]);
disp(['size(sol.idata.dif3d) = ' mat2str(size(sol.idata.dif3d))]);
disp(['length(sol.idata.kvec) = ' num2str(length(sol.idata.kvec))]);

% Try to call deval on it
try
    result = deval(sol, 0);
    disp('Success with t=0');
    disp(result);
catch e
    disp(['Error with t=0: ' e.message]);
    for i = 1:min(3, length(e.stack))
        fprintf('  %s line %d\n', e.stack(i).name, e.stack(i).line);
    end
end

try
    result = deval(sol, epsilon/2);
    disp('Success with t=epsilon/2');
    disp(result);
catch e
    disp(['Error with t=epsilon/2: ' e.message]);
    for i = 1:min(3, length(e.stack))
        fprintf('  %s line %d\n', e.stack(i).name, e.stack(i).line);
    end
end

% Compare with real ode15s solution
disp(' ');
disp('For comparison, real ode15s solution:');
real_sol = ode15s(@(t,y) -y, [0 1e-10], ic);
disp(['real_sol.x = ' mat2str(real_sol.x)]);
disp(['size(real_sol.y) = ' mat2str(size(real_sol.y))]);
disp(['size(real_sol.idata.dif3d) = ' mat2str(size(real_sol.idata.dif3d))]);
disp(['length(real_sol.idata.kvec) = ' num2str(length(real_sol.idata.kvec))]);

try
    result = deval(real_sol, 0);
    disp('Real solution: Success with t=0');
catch e
    disp(['Real solution: Error with t=0: ' e.message]);
end
