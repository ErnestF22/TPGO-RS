%Execute commands, and gives error if command fails
function systemOrError(cmd,varargin)
flagEscape=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'interpretescape'
            flagEscape=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
if flagEscape
    cmd=sprintf(cmd);
end

[b, output]=system(cmd);
if b>0
    disp(output)
    error('Execution of "%s" failed',cmd)
end

