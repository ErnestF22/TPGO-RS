function vargout=saveAscii(fileName,M)
[~,~,ext]=fileparts(fileName);
if isempty(ext)
    fileName=[fileName '.txt'];
elseif length(ext)==1
    fileName=[fileName '.'];
end

fid=fopen(fileName,'wt');
if fid<0
    error('Error opening file %s', fileName);
else
    fprintf(fid, '%+.17e ', M(:));
end

if nargout>0
    vargout=fileName;
end
