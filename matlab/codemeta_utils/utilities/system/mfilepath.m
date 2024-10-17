%Return path of current mfile (usage similar to mfilename('fullpath'))
function p=mfilepath()
%get calling stack
s=dbstack('-completenames');
%s(1) is the current function
if length(s)==1
    %called from the command line
    p=pwd();
else
    %return only path of the calling file
    p=fileparts(s(2).file);
end


