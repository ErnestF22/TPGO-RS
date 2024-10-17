%function ho_title(fid,title,level)
%Print a title in an HTML file
%Arguments:
%   fid     file handle to write to
%   title   string of text to insert
%   level   level of the title. Default: 1

%%AUTORIGHTS%%

function ho_title(fid,title,level)

if(exist('level','var')==0)
    level=1;
end
if(level<0 || fix(level)~=level)
    error('Level must be integer positive, it was %d', level);
end
fprintf(fid,'<h%d>%s</h%d>\n',level,title,level);
