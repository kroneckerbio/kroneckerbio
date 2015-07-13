function tests = UT14_AnalyticModels(useMEX)
% More extensive testing of the creation and simulation of analytic models
% from symbolic models
if nargin < 1
    useMEX = false;
end

% Get function handles
funnames = {'testBasicModel';'testNoStates';'testOneState';'testNoInputs';'testOneInput';'testNoOutputs';'testOneOutput';'testArbitraryOutput';'testModelUpdate'};
if useMEX
    funnames = strcat(funnames,'MEX');
end
testfuns = cellfun(@str2func,funnames,'UniformOutput',false);
tests = functiontests(testfuns);

end

%% Non-MEX functions

function testBasicModel(a)

useMEX = false;

opts = struct;

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testNoStates(a)

useMEX = false;

opts.Inputs = {'A';'B';'C';'D'};

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testOneState(a)

useMEX = false;

opts.Inputs = {'B';'C';'D'};

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end


end

function testNoInputs(a)

useMEX = false;

opts.Inputs = [];

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testOneInput(a)

useMEX = false;

opts.Inputs = {'A'};

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testNoOutputs(a)

useMEX = false;

opts.Outputs = [];

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testOneOutput(a)

useMEX = false;

opts.Outputs = {'A+B'};

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testArbitraryOutput(a)

useMEX = false;

opts.Outputs = {'koff^kon*(A^2+sqrt(B))/(C*DD)'};

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testModelUpdate(a)

useMEX = false;

opts = struct;

[m,expectedexprs] = getModel(opts,useMEX);

newk = m.k-[0.1;0.2;0.4;0.7];
m = m.Update(newk);

order = 0;
[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs,order);

% Get f and r indices. Don't test y changes because I set x and u to
% particular values, so parameter changes won't affect y.
vis = find(ismember(funnames,{'f';'r'}));

for vi = vis(:)'
    a.verifyGreaterThan(max(abs(funvals{vi} - expectedvals{vi})), 1e-9, [funnames{vi} ' did not change after model parameter update.']);
end
a.verifyEqual(m.k,newk,'Model parameters did not change following parameter update.')

end

%% MEX functions

function testBasicModelMEX(a)

useMEX = true;

opts = struct;

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testNoStatesMEX(a)

useMEX = true;

opts.Inputs = {'A';'B';'C';'D'};

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testOneStateMEX(a)

useMEX = true;

opts.Inputs = {'B';'C';'D'};

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end


end

function testNoInputsMEX(a)

useMEX = true;

opts.Inputs = [];

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testOneInputMEX(a)

useMEX = true;

opts.Inputs = {'A'};

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testNoOutputsMEX(a)

useMEX = true;

opts.Outputs = [];

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testOneOutputMEX(a)

useMEX = true;

opts.Outputs = {'A+B'};

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testArbitraryOutputMEX(a)

useMEX = true;

opts.Outputs = {'koff^kon*(A^2+sqrt(B))/(C*DD)'};

[m,expectedexprs] = getModel(opts,useMEX);

[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs);

for vi = 1:length(funvals)
    a.verifyEqual(full(funvals{vi}),expectedvals{vi},'AbsTol',1e-9,'RelTol',1e-6,[funnames{vi} ' differed from its expected value.']);  
end

end

function testModelUpdateMEX(a)

useMEX = true;

opts = struct;

[m,expectedexprs] = getModel(opts,useMEX);

newk = m.k-[0.1;0.2;0.4;0.7];
m = m.Update(newk);

order = 0;
[funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs,order);

% Get f and r indices. Don't test y changes because I set x and u to
% particular values, so parameter changes won't affect y.
vis = find(ismember(funnames,{'f';'r'}));

for vi = vis(:)'
    a.verifyGreaterThan(max(abs(funvals{vi} - expectedvals{vi})), 1e-9, [funnames{vi} ' did not change after model parameter update.']);
end
a.verifyEqual(m.k,newk,'Model parameters did not change following parameter update.')

end

%% Auxiliary functions

function [m, expectedexprs] = getModel(opts,useMEX)

if nargin < 2
    useMEX = false;
    if nargin < 1
        opts = struct;
    end
end

% Clear mex functions to allow overwriting of old ones, if necessary
clear mex

% Get symbolic model
[m, expectedexprs] = analytic_model_syms(opts);

% Set up mex directory
buildopts.UseMEX = useMEX;
kroneckerdir = cd(cd([fileparts(which('FinalizeModel.m')) filesep '..']));
buildopts.MEXDirectory = fullfile(kroneckerdir,'Testing','mexfuns');

% Compile MEX functions, if necessary
if useMEX
    warning('off','MATLAB:mex:GccVersion_link');
    compileMEXFunctions(buildopts.MEXDirectory,false)
    warning('on', 'MATLAB:mex:GccVersion_link');
end

end

function [funvals, expectedvals, funnames] = evaluateModel(m, expectedexprs, order)

if nargin < 3
    order = 2;
end

% Get model function names
if order == 2
    funnames = setdiff(fieldnames(expectedexprs),{'x','u','k','s','x0','dx0ds','d2x0ds2'});
    funnames = funnames(:);
    funnames_x0 = {'x0';'dx0ds';'d2x0ds2'};
elseif order == 1
    funnames = {'f';'r';'y';'dfdx';'dfdk';'dfdu';'drdx';'drdk';'drdu';'dydx';'dydu';'dydk'};
    funnames_x0 = {'x0';'dx0ds'};
elseif order == 0
    funnames = {'f';'r';'y'};
    funnames_x0 = {'x0'};
end

% Get model function handles
funlist = cellfun(@(fun)m.(fun),funnames,'UniformOutput',false);
funlist_x0 = cellfun(@(fun)m.(fun),funnames_x0,'UniformOutput',false);

% Get t, x, u, and s used in the tests
t = 10; % doesn't matter what t is as long as f or r don't have t in their expressions
x = expectedexprs.x;
u = expectedexprs.u;
s = expectedexprs.s;

% Evaluate the functions
evalfun = @(fun) fun(t,x,u);
funvals = cellfun(evalfun,funlist,'UniformOutput',false);
evalfun_x0 = @(fun) fun(s);
funvals_x0 = cellfun(evalfun_x0,funlist_x0,'UniformOutput',false);

% Concatenate f, r, y, and x0 functions
funvals = [funvals; funvals_x0];
funnames = [funnames; funnames_x0];

% Rearrange the expected values in a cell array the same size as funvals
expectedvals = cellfun(@(field)expectedexprs.(field),funnames,'UniformOutput',false);

end