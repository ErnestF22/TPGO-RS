function fputIntegerLine(fid,v)
s=repmat('%d ',1,length(v));
s(end)=[];
s=[s '\n'];
fprintf(fid,s,v);
