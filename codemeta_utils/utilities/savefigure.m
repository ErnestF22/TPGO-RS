%function savefigure(name,ext,dim,flagSaveFile)
%
%Save the current figure in a file. On unix systems, if extension is eps or
%epsc, and epstopdf is available, then generate also a pdf version
%Arguments:
%   name            name of the file
%   ext             extension (determines format, see drivers for 'print'
%                   command) Default: 'png'
%   dim             2x1 array with dimension in pixels. Use 'savefigureResolution'
%                   command for approximate conversion to length units.
%                   Default: current dimensions
%   flagSaveFile    if 0, run without actually writing the file
%                   if 1, run and save the file (default)
%                   if 2, same as 1, but also print name
%
%   See also savefigureResolution, print
%

%%AUTORIGHTS%%

function savefigure(name,ext,dim,flagSaveFile)
if(exist('ext','var')==0)
    ext='png';
end
rect=get(gcf,'Position');
if(exist('dim','var')==0)
    dim=rect(3:4);
end
if(exist('flagSaveFile','var')==0)
    flagSaveFile=1;
end
rect(3)=dim(1);
rect(4)=dim(2);
set(gcf,'Units','inches');
set(gcf,'Position',rect);
set(gcf,'PaperPositionMode','auto')
set(gcf,'PaperUnits','inches')
set(gcf,'PaperSize',rect(3:4))
if(flagSaveFile)
    eval(['print -d' ext ' ' name]);
    if(flagSaveFile==2)
        disp([name '.' ext])
    end
    if(strcmp(ext,'eps') || strcmp(ext,'epsc'))
        for cmd={'epspdf','epstopdf'}
            % Note: 'which' will return 0 if the command was found
            if ispc()
                [flagToolNotAvailable,o]=system(['where ' cmd{1}]);
            else
                [flagToolNotAvailable,o]=system(['which ' cmd{1}]);
            end
            if ~flagToolNotAvailable
                epspdfCommand=cmd{1};
                break
            end
        end
        if ~flagToolNotAvailable
            [flagConversionFail]=system([epspdfCommand ' ' name '.eps']);
            if flagConversionFail
                warning('Could not convert to pdf')
            end
            if(flagSaveFile==2)
                disp([name '.pdf'])
            end
        else
            warning('eps to pdf command line utility not available')
        end
    end
end
