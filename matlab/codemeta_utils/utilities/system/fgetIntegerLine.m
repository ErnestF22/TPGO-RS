%function varargout=fgetIntegerLine(fid)
%Read a line of intergers using the file identifier fid.
%Forms:
%   [a,b,c,..]=fgetIntegerLine(fid)
%       Each variable (a,b,c,...) gets one integer. A warning is issued if
%       there line in the file contains more variables than the number of
%       output arguments.
%   a=fgetIntegerLine(fid)
%       The output argument is an array containing all the numbers read
%The function returns empty arrays if an End Of File is encountered.
function varargout=fgetIntegerLine(fid)
l=fgetNonEmptyLine(fid);
if ~ischar(l)
    %EOF, return empty arrays
    varargout=cell(1,nargout);
else
    varargout{1}=sscanf(l,'%d');
    if nargout>1
        if length(varargout{1})>nargout
            warning('Number of elements in the line is greater then the number of output arguments');
        end
        varargout=num2cell(varargout{1});
    end
end