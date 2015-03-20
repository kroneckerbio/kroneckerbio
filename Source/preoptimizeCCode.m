function preoptimizeCCode(codefile)

%error('Function is not yet complete!')

% Read code file into string
codestr = fileread(codefile);

% Use regexp to capture all the terms
fprintf('Capturing terms...')
terms = regexp(codestr,'T\[\d+\]\s*=\s*([^;])+;','tokens');
terms = [terms{:}];
Ts = regexp(codestr,'T\[(\d+)\]\s*=\s*','tokens');
Ts = [Ts{:}];
termfunname = 'termfun';
fprintf('Done.\n')
if isempty(Ts)
    newcodestr = codestr;
elseif any(~cellfun(@isempty,regexp(terms,termfunname)))
    fprintf('Code %20s already appears to be preoptimized. Writing back original code.\n',codefile)
    newcodestr = codestr;
else
    nTswithzeros = str2double(Ts{end})+1; % +1 accounts for the fact that T[0] is the first element
    ntermswithzeros = length(terms);
    assert(ntermswithzeros==nTswithzeros, 'Number of expressions found does not match number of expressions in c file.')

    % Filter out the zeros
    iszero = ~cellfun(@isempty,regexp(terms,'^0(.0)?$'));
    terms = terms(~iszero);
    Ts = Ts(~iszero);
    nterms = length(terms);
    nTs = length(Ts);

    %% Get patterns that define the inline functions to be written

    fprintf('Capturing variables in terms...')
    termsvars = regexp(terms,'[kxu]\[\d+\]','match');
    fprintf('Done.\n')
    fprintf('Removing indices from variables in terms...')
    terms_novarnums = regexprep(terms,'([kxu])\[\d+\]','$1\[\]');
    fprintf('Done.\n')

    % Find groups of matching structures to find which inline functions are
    % needed
    [~,Iterms,Iunique] = unique(terms_novarnums,'stable');
    uniqueterms = terms(Iterms);
    uniquetermsvars = termsvars(Iterms);

    % Construct each funoutput
    funoutputs = uniqueterms;
    fprintf('Constructing inline functions...')
    fprintf(repmat(' ',1,16))
    for ui = 1:length(uniqueterms)
        fprintf([repmat('\b',1,16) '%5.3g%% complete\n'],ui/length(uniqueterms)*100)
        for ti = 1:length(uniquetermsvars{ui})
            replaceterm = regexprep(uniquetermsvars{ui}{ti},'([kxu])\[\d+\]',['$1\[ii\[' num2str(ti-1) '\]\]']);
            funoutputs{ui} = regexprep(funoutputs{ui},regexptranslate('escape',uniquetermsvars{ui}{ti}),replaceterm,1); %'$1\[$1i\[${num2str(find(strcmp($0,uniquetermsvars{ui})))}\]\]');
        end
    end

    % Construct each list of index strings
    iistr = cell(nTs,1);
    fprintf('Constructing index strings...')
    fprintf(repmat(' ',1,17))
    for ii = 1:length(termsvars)
        fprintf([repmat('\b',1,16) '%5.3g%% complete\n'],ii/length(termsvars)*100)
        thisiistr = regexp(termsvars{ii},'\d+','match');
        thisiistr = [thisiistr{:}];
        thisiistr = [thisiistr; repmat({','},1,length(thisiistr))];
        thisiistr = thisiistr(:);
        thisiistr = [thisiistr{:}];
        thisiistr = thisiistr(1:end-1);
        iistr{ii} = thisiistr;
    end

    % Construct inline functions
    funnames = [repmat({termfunname},1,length(uniqueterms)); cellfun(@num2str,num2cell(1:length(uniqueterms)),'UniformOutput',false)];
    funnames = strcat(funnames(1,:),funnames(2,:));
    funstr = cell(size(uniqueterms));
    fprintf('Constructing inline functions...')
    fprintf(repmat(' ',1,16))
    for fi = 1:length(uniqueterms)
        fprintf([repmat('\b',1,16) '%5.3g%% complete\n'],fi/length(uniqueterms)*100)
        funstr{fi} = functiontemplate(funnames{fi},funoutputs{fi});
    end

    % Construct function calls
    initstr = cell(nTs,1);
    callstr = cell(nTs,1);
    fprintf('Constructing function calls...')
    fprintf(repmat(' ',1,16))
    for ti = 1:nTs
        fprintf([repmat('\b',1,16) '%5.3g%% complete\n'],ti/nTs*100)
        ui = Iunique(ti);
        initstr{ti} = inittemplate(ti,funnames{ui},iistr{ti});
        callstr{ti} = calltemplate(ti,funnames{ui},iistr{ti});
    end

    % Add inline functions before main function
    newcodestr = codestr;
    newcodestr = regexprep(newcodestr,'void (([rf])|(d(\d*\w)d(\w\d*)+))', [funstr{:} 'void $1']);

    % Replace T terms with functionalized forms
    newcodestr = regexprep(newcodestr,'(?<=/\*\sStart\smodified\sccode\s\*/).*(?=/\*\sEnd\smodified\sccode\s\*/)',[initstr{:} '\n' callstr{:}]);

    % Remove zero entries from ir and jc
    irentries = regexp(codestr,'ir\[\]\s*=\s*\{([\d,]*)\}','tokens');
    irentries = irentries{1}{1};
    ir = eval(['[' irentries ']']);
    jcentries = regexp(codestr,'jc\[\]\s*=\s*\{([\d,]*)\}','tokens');
    jcentries = jcentries{1}{1};
    jc = eval(['[' jcentries ']']);
    newir = ir(~iszero);
    newjc = zeros(size(jc));
    for jj = 1:length(jc)-1
        iwithzero = (jc(jj)+1):jc(jj+1);
        newjc(jj+1) = newjc(jj) + sum(~iszero(iwithzero));
    end
    newirentries = [cellfun(@num2str,num2cell(newir),'UniformOutput',false); repmat({','},1,length(newir))];
    newirentries = [newirentries{:}];
    newjcentries = [cellfun(@num2str,num2cell(newjc),'UniformOutput',false); repmat({','},1,length(newjc))];
    newjcentries = [newjcentries{:}];
    newcodestr = regexprep(newcodestr,'(ir\[\]\s*=\s*\{)[\d,]*(\})',['$1' newirentries '$2']);
    newcodestr = regexprep(newcodestr,'(jc\[\]\s*=\s*\{)[\d,]*(\})',['$1' newjcentries '$2']);

    % Adjust NZMAX to account for removed zeros
    newcodestr = regexprep(newcodestr,'(#define NZMAX )\d+','$1${num2str(length(newir))}');
end

% Rename old file, if necessary
notoptimizedsuffix = '_notoptimized.c';
if isempty(regexp(codefile,[regexptranslate('escape',notoptimizedsuffix) '$']))
    [pathstr,name,~] = fileparts(codefile);
    if exist([pathstr filesep 'NotOptimized'],'dir') ~= 7; mkdir([pathstr filesep 'NotOptimized']); end;
    movefile(codefile,[pathstr filesep 'NotOptimized' filesep name '_notoptimized.c'])
    fid = fopen(codefile,'w');
    fprintf(fid,newcodestr);
    fclose(fid);
else
    newfilename = [codefile(1:(end-length(notoptimizedsuffix))) '.c'];
    fid = fopen(newfilename,'w');
    fprintf(fid,newcodestr);
    fclose(fid);
end

    function out = inittemplate(ti,~,iistr)
        out = ['static const int ii' num2str(ti) '[] = {' iistr '};\n'];
    end

    function out = calltemplate(ti,funname,~)
        out = ['T[' num2str(ti-1) '] = ' funname '(t,x,u,k,ii' num2str(ti) ');\n'];
    end

    function out = functiontemplate(funname,funoutput)
        out = [ ...
            %'inline double ' funname '(double t[], double x[], double u[], double k[], const int ii[]);\n' ...
            'inline double ' funname '(double t[], double x[], double u[], double k[], const int ii[])\n' ...
            '{\n' ...
            '   return ' funoutput ';\n' ...
            '}\n\n' ...
            ];
    end

end
%% Old code
% % Write pow function
% pow = @(x,p) x^p;
% 
% % Convert square brackets into parentheses, and add one to the number in
% % between
% fprintf('Converting square brackets and indices...      \n')
% for ti = 1:nterms
%     fprintf('\b\b\b\b\b\b%6i',ti)
%     %terms{ti} = regexprep(terms{ti},'\[(\d+)\]','\(${num2str(str2num($1)+1)}\)');
%     terms{ti} = regexprep(terms{ti},'\[(\d+)\]','\(${sprintf(''%i'',sscanf($1,''%i'')+1)}\)');
% end
% 
% % Initialize symbolic vectors
% nk = regexp(codestr,'#define KROWS (\d+)','tokens');
% nk = [nk{:}];
% nk = str2double(nk{1});
% k = sym('k',[nk 1]);
% nx =  regexp(codestr,'#define XROWS (\d+)','tokens');
% nx = [nx{:}];
% nx = str2double(nx{1});
% x = sym('x',[nx 1]);
% nu =  regexp(codestr,'#define UROWS (\d+)','tokens');
% nu = [nu{:}];
% nu = str2double(nu{1});
% u = sym('u',[nu 1]);
% 
% % Evaluate terms to convert to symbolic
% termssymcell = cell(nterms,1);
% for ti = 1:length(terms)
%     fprintf('%5g\n',ti)
%     termssymcell{ti} = eval(terms{ti});
% end
% 
% disp();
% 
% % Split nonzero terms into numerators and denominators
% 
% % Simplify the numerators and denominators using the Horner representation,
% % which is good for calculation
% 
% % Reconstruct the numerators and denominators

% Construct each function call

% Check: same structure?
% fprintf('Finding terms with same structure...')
% matchingstructures = multistrcmp(terms_novarnums,terms_novarnums);
% matchingstructures = tril(matchingstructures,-1);
% [matchi,matchj] = find(matchingstructures);

% Check 2: same number pattern?
% termsvarspatterns = cellfun(@getpattern,termsvars);
% matchingstructuresandpatterns = false(nterms);
% matchingstructuresandpatterns(matchingstructures) = cellfun(@isequal,termsvarspatterns(matchi),termsvarspatterns(matchj));

%structurestofunctionalize = any(matchingstructures);

%     function pattern = getpattern(vars)
%         uvars = unique(vars,'stable');
%         pattern = multistrcmp(uvars,vars);
%     end
