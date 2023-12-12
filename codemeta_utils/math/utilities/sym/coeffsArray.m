%Array version of the coeffs function in the Symbolic Toolbox
%function [c,t]=coeffsArray(p,x)
function [c,t]=coeffsArray(p,x,varargin)
flagProgressBar=false;

ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'progressbar'
            flagProgressBar=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

sz=size(p);
NEl=numel(p);
cCell=cell(1,NEl);
tCell=cell(1,NEl);
tCellStr=cell(1,NEl);

if flagProgressBar
    w=getTextWaitBar(NEl);
    w(0);
end
for iEl=1:NEl
    [cEl,tEl]=coeffs(p(iEl),x);
    cCell{iEl}=cEl;
    tCell{iEl}=tEl;
    tCellStr{iEl}=sym2char(tEl);
end

t=tCell{1};
tStr=tCellStr{1};
c=cell(1,NEl);
c{1}=cCell{1};
for iEl=2:NEl
    %coefficients/embedding to compare
    tNew=tCell{iEl};
    tNewStr=tCellStr{iEl};
    cNew=cCell{iEl};
    
    %find additional embedding terms
    [tAddStr,idxAddNew]=setdiff(tNewStr,tStr,'rows');
    NAdd=length(idxAddNew);
    
    %find mapping between new and current embedding terms
    [tIntersectStr,idxIntersectNew,idxIntersectCurrent]=intersect(tNewStr,tStr,'rows');
    
    %prepare coeffs vector to add
    cNext=[sym(zeros(size(t))) cNew(idxAddNew)];
    cNext(idxIntersectCurrent)=cNew(idxIntersectNew);
    
    %add coeffs vectors and additional embedding terms
    c{iEl}=cNext;
    t=[t tNew(idxAddNew)];
    if NAdd>0
        %append new string version of embedding
        tStr=char(tStr,tNewStr(idxAddNew,:));
        %pad previous coeffs vectors 
        z=zeros(1,NAdd);
        for jEl=1:iEl-1
            c{jEl}=[c{jEl} z];
        end
    end
    if flagProgressBar
        w(iEl)
    end
end
c=reshape(c,sz);

