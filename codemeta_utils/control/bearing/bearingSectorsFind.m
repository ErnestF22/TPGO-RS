%Find the closest unit vector that is not inside a sector
%function [uProj,flagHit,dProj]=bearingSectorsFind(g,y,a)
function [uProj,flagHit,dProj]=bearingSectorsFind(g,y,a)
[uProj,flagHit,dProj]=projectSectorsUpDownIterate(g,y,a,'up','down');
if size(uProj,2)==2
    [dProj,idxMinDProj]=min(dProj);
    uProj=uProj(:,idxMinDProj);
end
    


function [uProj,flagHit,dProj]=projectSectorsUpDownIterate(g,y,a,varargin)
flagUp=false;
flagDown=false;
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'up'
            flagUp=true;
        case 'down'
            flagDown=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

[~,flagHit]=projectSectorsUpDown(g,y,a,'up');

uProj=[];
dProj=[];
if flagUp
    gProjUp=g;
    dProjUpCum=0;
    cnt=0;
    while 1
        [gProjUp,flagHitUp,dProjUp]=projectSectorsUpDown(gProjUp,y,a,'up');
        dProjUpCum=dProjUpCum+dProjUp;
        if ~flagHitUp || dProjUpCum>=pi
            break
        end
        cnt=cnt+1;
        if cnt>1000
            keyboard
        end
    end
    if dProjUpCum<pi
        uProj=[uProj gProjUp];
        dProj=[dProj dProjUpCum];
    end
end

if flagDown
    gProjDown=g;
    dProjDownCum=0;
    while 1
        [gProjDown,flagHitDown,dProjDown]=projectSectorsUpDown(gProjDown,y,a,'down');
        dProjDownCum=dProjDownCum+dProjDown;
        if ~flagHitDown || dProjDownCum>=pi
            break
        end
    end
    if dProjDownCum<pi
        uProj=[uProj gProjDown];
        dProj=[dProj dProjDownCum];
    end
end



function [uProj,flagHit,dProj]=projectSectorsUpDown(g,y,a,varargin)
flagUp=false;
flagDown=false;
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'up'
            flagUp=true;
        case 'down'
            flagDown=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

ay=atan2(y(2,:),y(1,:));
ayMax=modAngle(ay+a);
ayMin=modAngle(ay-a);

ag=atan2(g(2),g(1));
idxHit=find(and(angleDiff(ag,ayMin)>2*eps,angleDiff(ayMax,ag)>2*eps));
flagHit=~isempty(idxHit);

if ~flagHit
    uProj=cnormalize(g);
    dProj=0;
else
    uProj=[];
    dProj=[];
    if flagUp
        [dHitMax,idxIdxHitMax]=max(angleDiff(ayMax(idxHit),ag));
        idxHitMax=idxHit(idxIdxHitMax);
        uProj=[uProj angle2Vector(ayMax(idxHitMax))];
        dProj=[dProj dHitMax];
    end
    if flagDown
        [dHitMin,idxIdxHitMin]=max(angleDiff(ag,ayMin(idxHit)));
        idxHitMin=idxHit(idxIdxHitMin);
        uProj=[uProj angle2Vector(ayMin(idxHitMin))];
        dProj=[dProj dHitMin];
    end
end

function y=angle2Vector(a)
y=[cos(a); sin(a)];
