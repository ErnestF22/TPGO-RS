function POCSJC
N=4;
p=1:N;
d=-ones(1,N);
d(1)=0;

disp([num2str(p,' %02d'); num2str(d,' %+2d')])
disp(' ')
for n=1:factorial(N)-1
    idxNotZero=find(d~=0);
    [~,idxIdxMax]=max(p(idxNotZero));
    selIdx=idxNotZero(idxIdxMax);
    selDir=d(selIdx);
    selIdxNext=selIdx+selDir;
    selVal=p(selIdx);
    
    [p(selIdx),p(selIdxNext)]=swap(p(selIdx),p(selIdxNext));
    [d(selIdx),d(selIdxNext)]=swap(d(selIdx),d(selIdxNext));
    if selIdxNext==1 || selIdxNext==N || p(selIdxNext+selDir)>p(selIdxNext)
        d(selIdxNext)=0;
    end
    
    idxLeft=1:selIdxNext;
    idxRight=selIdxNext:N;
    d(idxLeft(p(idxLeft)>selVal))=+1;
    d(idxRight(p(idxRight)>selVal))=-1;
    
    disp([num2str(p,' %02d'); num2str(d,' %+2d')])
    disp(' ')
end
end

function [b,a]=swap(a,b)
end

  
function POCSJC2
  N=3;
  %Init
  dir=-ones(N);
  p=1:N;
  pi=1:N;
  %call
  for a=1:6
      [p,pi,dir]=perm(p,pi,dir,2);
  end
end

function [p,pi]=move(p,pi,x,d)
  z = p(pi(x)+d);
  p(pi(x)) = z;
  p(pi(x)+d) = x;
  pi(z) = pi(x);
  pi(x) = pi(x)+d;  
end

function [p,pi,dir]=perm(p,pi,dir,n)
    if n > length(p)
       disp(p)
    else
        [p,pi,dir]=perm(p,pi,dir,n+1);
        for i=1:n-1
           [p,pi]=move(p,pi,n,dir(n));
           [p,pi,dir]=perm(p,pi,dir,n+1); 
        end
        dir(n) = -dir(n);
    end
end
