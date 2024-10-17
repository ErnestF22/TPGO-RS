%function sCat=catStructRecursive(dim,s1,s2,...)
%Obtain structure s12 by applying cat(dim,...) on each field pair. It is
%assumed that s1 and s2 have the same fields.
%If a field is itself a struct, the function calls itself recursively
%instead of concatenating.
%Example
%     >> a=struct('f',1,'g',2,'h',struct('ff','a'));
%     >> b=struct('f',3,'g',4,'h',struct('ff','b'));
%     >> c=struct('f',5,'g',6,'h',struct('ff','c'));
%     >> abc=catStructRecursive(2,a,b,c)
%     abc = 
%         f: [1 3 5]
%         g: [2 4 6]
%         h: [1x1 struct]
%     >> abc.h
%     ans = 
%         ff: 'abc'

function sCat=catStructRecursive(dim,varargin)
s=varargin;
Ns=length(s);
sCat=s{1};
fieldNames=fieldnames(sCat);
NFields=length(fieldNames);
for is=2:Ns
    s2=s{is};
    for iField=1:NFields
        fn=fieldNames{iField};
        s2field=s2.(fn);
        if isstruct(s2field)
            sCat.(fn)=catStructRecursive(dim,sCat.(fn),s2field);
        else
            sCat.(fn)=cat(dim,sCat.(fn),s2field);
        end
    end
end
