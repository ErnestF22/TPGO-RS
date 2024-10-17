%function h=getSeqPerm(v)
%Returns a handle to a function h() which, when called, sequentially
%generates all permutations of v
%
%Example
%
%  generator=getSeqPerm(1:5);
%
%Then
%
%  for(it=1:factorial(5)); disp(generator()); end
%
%displays all the permutations of the numbers between 1 and 5

function h=getSeqPerm(v,varargin)
method='sjc';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'method'
            ivarargin=ivarargin+1;
            method=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

N=length(v);
s=1:N;
switch lower(method)
    case 'wiki'
        h=@wikiCounter;
    case 'sjc'
        sgn=1;
        d=-ones(1,N);
        d(1)=0;
        h=@sjcCounter;
    otherwise
        error('Method for generating permutations not reconized')
end

    function out=wikiCounter()
        %Algorithm for in-place generation of permutations from Wikipedia
        out=v(s);
        %find highest index i such that s(i)<s(i+1)
        ii=find(s(1:end-1)<s(2:end),1,'last');
        %find highest index j>i such that s(j)>s(i)
        jj=find(s(ii+1:end)>s(ii),1,'last')+ii;
        %swap s(i) and s(j)
        s([ii jj])=s([jj ii]);
        %reverse the order of all elements after index i
        s(ii+1:end)=fliplr(s(ii+1:end));
    end

    function [out,p]=sjcCounter()
        out=v(s);
        p=sgn;
        %Steinhaus?Johnson?Trotter algorithm with Even's speedup
        sgn=sgn*-1;
        idxNotZero=find(d~=0);
        if ~isempty(idxNotZero)
            [~,idxIdxMax]=max(s(idxNotZero));
            selIdx=idxNotZero(idxIdxMax);
            selDir=d(selIdx);
            selIdxNext=selIdx+selDir;
            selVal=s(selIdx);

            [s(selIdx),s(selIdxNext)]=swap(s(selIdx),s(selIdxNext));
            [d(selIdx),d(selIdxNext)]=swap(d(selIdx),d(selIdxNext));
            if selIdxNext==1 || selIdxNext==N || s(selIdxNext+selDir)>s(selIdxNext)
                d(selIdxNext)=0;
            end

            idxLeft=1:selIdxNext;
            idxRight=selIdxNext:N;
            d(idxLeft(s(idxLeft)>selVal))=+1;
            d(idxRight(s(idxRight)>selVal))=-1;
        else
            s=1:N;
        end
    end
end

function [b,a]=swap(a,b)
end
