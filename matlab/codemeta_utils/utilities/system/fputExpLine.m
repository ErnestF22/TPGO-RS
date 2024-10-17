function fputExpLine(fid,v)
s=repmat('%.10e ',1,length(v));
s(end)=[];
s=[s '\n'];
fprintf(fid,s,v);
