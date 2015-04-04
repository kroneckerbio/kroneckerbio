function dummyfunlist = getDummyFunctions(mexdir)

% Get function names
filestruct = dir(mexdir);
names = {filestruct.name};
mexnames = names(~cellfun(@isempty,regexp(names,'\.mex')));
mexfuns = regexp(mexnames,'(.*)\.mex','tokens');
mexfuns = [mexfuns{:}];
mexfuns = [mexfuns{:}];

% Run functions and examine error messages for the phrase "dummy function"
% to determine whether functions are dummies
isdummy = false(size(mexfuns));
for fi = 1:length(mexfuns)
    try
        eval([mexfuns{fi} '(1,1,1,1)'])
    catch ME
        if isempty(regexp(ME.message,'dummy function','ONCE'))
            isdummy(fi) = false;
        else
            isdummy(fi) = true;
        end
    end
end

dummyfunlist = mexfuns(isdummy);