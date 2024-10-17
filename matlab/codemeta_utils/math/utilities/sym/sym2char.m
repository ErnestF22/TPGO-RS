%Convert an array of symbolic vars to a char array
%function ca=sym2char(a,varargin)
%Optional arguments
%   'l',l           initial number of columns of ca (final number could be
%                   larger, depending on the entries of a) 
%   'progressBar'   show progress bar during conversion
function ca=sym2char(a,varargin)
l=1;
flagProgressBar=false;

ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'progressbar'
            flagProgressBar=true;
        case 'l'
            ivarargin=ivarargin+1;
            l=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

Na=numel(a);
ca=repmat(' ',Na,l);
if flagProgressBar
    w=getTextWaitBar(Na);
    w(0)
end
for ia=1:Na
    s=char(a(ia));
    ca(ia,1:length(s))=s;
    if flagProgressBar
        w(ia)
    end
end

ca(ca==0)=' ';
