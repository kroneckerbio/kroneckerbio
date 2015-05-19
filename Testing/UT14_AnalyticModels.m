function tests = UT14_AnalyticModels(useMEX)
% More extensive testing of the creation and simulation of analytic models
% from symbolic models
if nargin < 1
    useMEX = false;
end

% Get function handles
funnames = {'testBasicModel';'testNoStates';'testOneState';'testNoInputs';'testOneInput';'testNoOutputs';'testOneOutput';'testModelUpdate'};
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

function [m,expectedexprs] = getModel(opts,useMEX)

if nargin < 2
    useMEX = false;
    if nargin < 1
        opts = struct;
    end
end

% Clear mex functions to allow overwriting of old ones, if necessary
clear mex

% Get symbolic model
[m,expectedexprs] = symbolicmodel(opts);

% Set up mex directory
buildopts.UseMEX = useMEX;
kroneckerdir = cd(cd([fileparts(which('symbolic2PseudoKronecker.m')) filesep '..']));
buildopts.MEXDirectory = fullfile(kroneckerdir,'Testing','mexfuns');

% Remove old mex directory from the path, if it exists
% if exist(buildopts.MEXDirectory,'dir') == 7
%     id = 'MATLAB:rmpath:DirNotFound';
%     warning('off',id)
%     rmpath(buildopts.MEXDirectory);
%     warning('on',id)
%     %rmsuccess = rmdir(buildopts.MEXDirectory,'s');
% end

% Build analytical model
m = symbolic2PseudoKronecker(m,buildopts);

% Compile MEX functions, if necessary
if useMEX
    compileMEXFunctions(buildopts.MEXDirectory,false)
end

end

function [funvals,expectedvals,funnames] = evaluateModel(m,expectedexprs,order)

if nargin < 3
    order = 2;
end

% Get model function names
if order == 2
    funnames = setdiff(fieldnames(expectedexprs),{'x','u','k'});
elseif order == 1
    funnames = {'f';'r';'y';'dfdx';'dfdk';'dfdu';'drdx';'drdk';'drdu';'dydx';'dydu'};
elseif order == 0
    funnames = {'f';'r';'y'};
end

% Get model function handles
funlist = cellfun(@(fun)m.(fun),funnames,'UniformOutput',false);

% Get t, x, and u used in the tests
t = 10; % doesn't matter what t is as long as f or r don't have t in their expressions
x = expectedexprs.x;
u = expectedexprs.u;

% Evaluate the functions
evalfun = @(fun) fun(t,x,u);
funvals = cellfun(evalfun,funlist,'UniformOutput',false);

% Rearrange the expected values in a cell array the same size as funvals
expectedvals = cellfun(@(field)expectedexprs.(field),funnames,'UniformOutput',false);

end