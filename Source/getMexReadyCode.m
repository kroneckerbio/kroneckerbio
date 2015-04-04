function getMexReadyCode(dsym,nzi,nzsizes,xSyms,uSyms,kSyms,dername,savedir)

% Set whether to pre-optimize code or not
preoptimize = true;

if preoptimize
    fprintf(['Writing preoptimized MEX code for ' dername '...\n'])
else
    fprintf(['Writing MEX code for ' dername '...\n'])
end

dsymstr = 'T'; % when drawing one column at a time from dsym, ccode uses the variable 'T'
m = size(dsym,1);
n = size(dsym,2);

fprintf('Getting C code...')
cstrcell = cell(size(nzi,1),1);
% for coli = 1:n
%     cstrtemp = ccode(dsym(:,coli));
%     if ~isempty(cstrtemp)
%         cstrtemp = [cstrtemp '\n'];
%     end
%     cstrcell{coli} = regexprep(cstrtemp,[dsymstr '\[(\d+)\]\[\d+\]'],[dsymstr '\[$1\]\[' num2str(coli-1) '\]']);
%     fprintf('%3g percent complete\n',100*coli/n)
% end

% sub2ind requires a size input of length two
if length(nzsizes) == 1
    nzsizes = [nzsizes 1];
end

fprintf('%25s', repmat(' ',1,25))
for nzii = 1:size(nzi,1)
    % Convert generalized indices to matrix indices
    nzicell = num2cell(nzi(nzii,:));
    thisnzi = sub2ind(nzsizes,nzicell{:});
    [r,c] = ind2sub(size(dsym),thisnzi);
    cstrtemp = ccode(dsym(r,c));
    if ~isempty(cstrtemp)
        cstrtemp = [cstrtemp '\n'];
    end
    cstrcell{nzii} = regexprep(cstrtemp,'t0',[dsymstr '\[' num2str(r-1) '\]\[' num2str(c-1) '\]']);
    fprintf([repmat('\b',1,25) '%25s'], sprintf('%5.4g%% complete',100*nzii/size(nzi,1)))
    %fprintf('%3g percent complete\n',100*nzii/size(nzi,1))
end
cstr = [cstrcell{:}];
fprintf('\n')

if ~isempty(cstr)
    fprintf('Getting indices and expressions...')
    ijsstruct = regexp(cstr,[dsymstr '\[(?<i>\d+)\]\[(?<j>\d+)\]\s*=\s*(?<s>[^;]+);'],'names');
    % str2numplus1 = @(str)(str2double(str));
    % for ii = 1:length(ijsstruct)
    %     ijsstruct(ii).i = str2numplus1(ijsstruct(ii).i);
    %     ijsstruct(ii).j = str2numplus1(ijsstruct(ii).j);
    % end
    fprintf('Done.\n')

    fprintf('Sorting elements by column...')
    indices = [cellfun(@str2double,{ijsstruct(:).i}') cellfun(@str2double,{ijsstruct(:).j}')];
    [indices,si] = sortrows(indices,2);
    ir = indices(:,1);
    jc = getjc(indices(:,2),n);
    fprintf('Done.\n')
    %pr = {ijsstruct(si).s}';

    fprintf('Writing ir and jc strings...')
    irstr = indexstr(ir);
    jcstr = indexstr(jc);
    fprintf('Done.\n')

    fprintf('Replacing output indices...')
    getoldindexstr = @(ii)([dsymstr '\[' ijsstruct(si(ii)).i '\]\[' ijsstruct(si(ii)).j '\]']);
    getnewindexstr = @(ii)([dsymstr '[' num2str(ii-1) ']']);
    for ii = 1:length(ijsstruct)
        cstr = regexprep(cstr,getoldindexstr(ii),getnewindexstr(ii));
    end
    fprintf('Done.\n')

    fprintf('Replacing variable names with indexed vector elements...')
    % Get new variable names
    varRepinit = @(varSyms)(num2cell((1:length(varSyms))'));
    kRep = varRepinit(kSyms);
    xRep = varRepinit(xSyms);
    uRep = varRepinit(uSyms);
    varnum2str = @(k,varname)([varname '[' num2str(k-1) ']']);
    knum2str = @(k)varnum2str(k,'k');
    xnum2str = @(x)varnum2str(x,'x');
    unum2str = @(u)varnum2str(u,'u');
    kRep = cellfun(knum2str,kRep,'UniformOutput',false);
    xRep = cellfun(xnum2str,xRep,'UniformOutput',false);
    uRep = cellfun(unum2str,uRep,'UniformOutput',false);
    % Convert old variable names to strings
    varsym2str = @(varSyms)(arrayfun(@char,varSyms,'UniformOutput',false));
    kStr = varsym2str(kSyms);
    xStr = varsym2str(xSyms);
    uStr = varsym2str(uSyms);
    % Replace variable strings
    cstr = replacevars(kStr,kRep,cstr);
    cstr = replacevars(xStr,xRep,cstr);
    cstr = replacevars(uStr,uRep,cstr);
    fprintf('Done.\n')
else
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
xrowsstr = num2str(length(xSyms));
urowsstr = num2str(length(uSyms));
krowsstr = num2str(length(kSyms));
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

% Preoptimize, if necessary
if preoptimize
    preoptimizeCCode([savedir filesep dername 'fun.c'])
end

    function strout = indexstr(ir)
        addcommas = @(x)([num2str(x) ',']);
        strout = arrayfun(addcommas,ir,'UniformOutput',false);
        strout = [strout{1:end}];
        strout = strout(1:end-1); % leave off last comma
    end

    
    function cstr = replacevars(varStr,varRep,cstr)
        for vi = 1:length(varStr)
            cstr = regexprep(cstr,varStr{vi},varRep{vi});
        end
    end

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

end
