%function l=fgetNonEmptyLine(fid)
%Same as fget but skips empty lines
function l=fgetNonEmptyLine(fid)
l=[];
while isempty(l)
    l=fgetl(fid);
end
