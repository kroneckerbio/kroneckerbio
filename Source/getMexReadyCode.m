function getMexReadyCode(dsym,nzi,nzsizes,xSyms,uSyms,kSyms,dername,savedir)

m = size(dsym,1);
n = size(dsym,2);

% sub2ind requires a size input of length two
if length(nzsizes) == 1
    nzsizes = [nzsizes 1];
end

nzicell = mat2cell(nzi,size(nzi,1),ones(size(nzi,2),1)); % Separate columns into separate cells
nziindices = sub2ind(nzsizes,nzicell{:});

% Get code
if ~isempty(nziindices)
    cstr = ccode(vec(dsym(nziindices)));
else
    cstr = '';
end

% If only one nonzero element is in dsym, the evaluation is assigned to t0
% instead of T[index][0]. Replace t0 with T[0][0] for consistency.
if length(nziindices) == 1
    cstr = regexprep(cstr,'t0 =','T[0][0] =');
end

if ~isempty(cstr)
    % Remove extra [0] element from each term
    cstr = regexprep(cstr,'\]\[0\]','\]');

    % Get ir and jc indices
    [rs,cs] = ind2sub(size(dsym),nziindices);
    rs = rs-1; cs = cs-1; % Adjust for 0 indexing
    [indices,rctoirjc] = sortrows([rs,cs],2);
    ir = indices(:,1);
    jc = getjc(indices(:,2),n);
    irstr = strjoin(row(strtrim(cellstr(int2str(ir)))),',');
    jcstr = strjoin(row(strtrim(cellstr(int2str(jc)))),',');

    % Reorder terms to coincide with ir-jc indices
    cstr = strtrim(strsplit(cstr,'\n'));
    for ii = 1:length(cstr)
        cstr{rctoirjc(ii)} = ...
            strrep(...
                cstr{rctoirjc(ii)},...
                ['[' rctoirjc(ii)-1 ']'],...
                ['[' num2str(ii)-1 ']']...
            );
    end
    cstr = cstr(rctoirjc);

    % Replace standardized variable names with x, u, and k vector elements
    nx = length(xSyms);
    nu = length(uSyms);
    nk = length(kSyms);
    getNewVars = @(var,nvar) strsplit(strtrim(sprintf([var '[%d]\n'],0:nvar-1)),'\n')';
    xnew = getNewVars('x',nx);
    unew = getNewVars('u',nu);
    knew = getNewVars('k',nk);
    if nx == 0
        xnew = {};
    end
    if nu == 0
        unew = {};
    end
    if nk == 0
        knew = {};
    end
    xold = fastchar(xSyms);
    uold = fastchar(uSyms);
    kold = fastchar(kSyms);
    cstr = regexprep(cstr,[xold;uold;kold],[xnew;unew;knew]);
    
    % Concatenate strings in cells
    cstr = strjoin(row(cstr),'\n');
    
else
    
    nx = length(xSyms);
    nu = length(uSyms);
    nk = length(kSyms);
    
    cstr = '';
    ir = [];
    irstr = '';
    jcstr = '';
    
end

% Read template file into string
codestr = fileread('dsymtemplate.c');
% Set up string replacement tags
nzmaxstr = num2str(length(ir));
rowsstr = num2str(m);
colsstr = num2str(n);
xrowsstr = num2str(nx);
urowsstr = num2str(nu);
krowsstr = num2str(nk);
tags = struct('tag',{'%NZMAX%';'%ROWS%';'%COLS%';'%XROWS%';'%UROWS%';'%KROWS%';'%CSTR%';'%IRSTR%';'%JCSTR%';'%FUN%'}, ...
    'stringname',   {nzmaxstr ; rowsstr; colsstr; xrowsstr; urowsstr; krowsstr; cstr;   irstr;      jcstr;  dername});
ntags = length(tags);

% Perform string replacements
for ti = 1:ntags
    codestr = regexprep(codestr, tags(ti).tag, tags(ti).stringname);
end

% Output string to file
fid = fopen([savedir filesep dername 'fun.c'],'w+');
fprintf(fid,codestr);
fclose(fid);

    function jc = getjc(jindices,n)
        jc = zeros(n+1,1);
        for ci = 2:(n+1)
            jctemp = find(jindices==ci-2,1,'last');
            if isempty(jctemp)
                jc(ci) = jc(ci-1);
            else
                jc(ci) = jctemp;
            end
        end
    end

    function fstr = fastchar(f)
        filename = which('fastchar.mu');
        read(symengine, filename);
        fstr = feval(symengine,'fastchar',f);
        fstr = eval(fstr);
    end

end
