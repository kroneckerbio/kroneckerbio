function compileMEXFunctions(mexdir,useparallel,clibpath)
    % compileMEXFunctions(mexdir,useparallel,Clibrarypath) compiles the
    % Kronecker model MEX functions in the directory specified by the
    % string mexdir, then adds the directory to the path. Use this function
    % after building a model using the UseMEX option in
    % symbolic2PseudoKronecker.

    if nargin < 3
        clibpath = [];
        if nargin < 2
            useparallel = [];
            if nargin < 1
                mexdir = [];
            end
        end
    end    

    if isempty(clibpath)
        clibpath = 'default';
    end
    if isempty(useparallel)
        useparallel = false;
    end
    if isempty(mexdir)
        mexdir = cd;
    end

        % Make sure mexdir ends in a file separator
%         if ~strcmp(mexdir,filesep)
%             mexdir = [mexdir filesep];
%         end

        % Come up with a list of MEX functions to be compiled
        D = dir(fullfile(mexdir, '*fun.c'));
        filenames = {D.name};

        % Generate compiling commands
        switch clibpath
            case 'default'
                disp('Linking to most recent version of C libraries available')
                compilecommand = @(dfun)(['mex -largeArrayDims CFLAGS="\$CFLAGS -std=c99" -outdir "' mexdir '" "' dfun '"']);
            otherwise
                disp('Linking to user-provided C libraries')
                compilecommand = @(dfun)(['mex -L' clibpath ' -lc -largeArrayDims CFLAGS="\$CFLAGS -std=c99" -outdir "' mexdir '" "' dfun '"']);
        end
        compile = @(dfun) eval(compilecommand(dfun));

        % Do the compilations
        if useparallel
            parfor ci = 1:length(filenames)
                disp(['Compiling ' filenames{ci}])
                compile(fullfile(mexdir, filenames{ci}))
                disp(['Done compiling ' filenames{ci}])
            end
        else
            for ci = 1:length(filenames)
                disp(['Compiling ' filenames{ci}])
                compile(fullfile(mexdir, filenames{ci}))
                disp(['Done compiling ' filenames{ci}])
            end
        end
    
        % Add the directory to the path
        addpath(mexdir)
end
