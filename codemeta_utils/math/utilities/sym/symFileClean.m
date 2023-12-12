function symFileClean(fileNameIn,fileNameOut,argumentStr)

[~, functionNameIn]=fileparts(fileNameIn);
[~, functionNameOut]=fileparts(fileNameOut);

fidIn=fopen(fileNameIn,'r');
fidOut=fopen(fileNameOut,'w');

while ~feof(fidIn)
    s=fgetl(fidIn);
    
    flagHeader=~isempty(regexp(s,['function \w* = ' functionNameIn  '\([\w,]*\)$'], 'once'));
    if ~flagHeader
        s=regexprep(s,'(\d+)_(\d+)','($1,$2)');
    else
        s=regexprep(s,'\([\w,]*\)',['(' argumentStr ')']);
        s=strrep(s,functionNameIn,functionNameOut);
    end
    fprintf(fidOut,'%s\n',s);
    if flagHeader
        fprintf(fidOut,'%%    Processed by symFileClean\n');
    end
end

fclose(fidIn);
fclose(fidOut);