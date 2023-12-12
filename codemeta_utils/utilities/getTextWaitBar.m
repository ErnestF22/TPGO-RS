%function bar=getTextWaitBar(maxx,dotInc,numInc)
%Returns the handle to a function BAR which is used to display a text
%waitbar such as
%
%   0%....50%....100%
%
%MAXX represents the number which represents 100% on the bar
%DOTINC represents how much each dot is worth (default 0.02)
%NUMINC represents at what ratio increments to show a number (default 0.1)
%By calling BAR(newx), the bar is up to the percentage newx/maxx
%Only increments are possible. If you display something between calls to
%BAR(), things will be messed up
%
%Example
%
%   bar=getTextWaitBar(100);
%   bar(0)
%   for it=1:100
%       pause(0.1)
%       bar(it)
%   end

%%AUTORIGHTS%%

function bar=getTextWaitBar(maxx,dotInc,numInc)
if(~exist('maxx','var'))
    maxx=1;
end
if(~exist('dotInc','var'))
    dotInc=0.02;
end
if(~exist('numInc','var'))
    numInc=0.1;
end

nextDot=1;
nextNum=1;
x=0;                %internal state between zero and one
finished=false;
firstCall=true;     %becomes false if increment is called at lest one

bar=@increment;

    function increment(newx)
        x=newx/maxx;
        update;
    end

    function update
        if firstCall
            fprintf('0%%')
            firstCall=false;
        end
        while(nextDot*dotInc<=x || nextNum*numInc <=x)
            while(nextDot*dotInc<=nextNum*numInc && nextDot*dotInc<=x)
                fprintf('.')
                nextDot=nextDot+1;
            end
            if(nextNum*numInc<=x)
                fprintf('%u%%',round(nextNum*numInc*100))
                nextNum=nextNum+1;
            end
        end
        if(x>=1 && ~finished)
            fprintf('\n')
            finished=true;
        end
    end
end
