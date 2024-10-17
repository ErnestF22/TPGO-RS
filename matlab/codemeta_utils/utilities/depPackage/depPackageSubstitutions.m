%function depPackageSubstitutions(nameList,pairs,varargin)
%PAIRS is a Nx2 list of pairs of strings. The function performs all the
%substitutions on all the strings in nameList.
%WARNING: substitutions are executed one after the other. Hence earlier
%substitutions might affect the ones that follow.
%
%Input arguments
%   pairs   [NPairs x 2] cell array with pairs of strings to substitute
%Optional arguments
%   'removeRoot'    find and remove the common root of all the strings.
%                   This substitution is made BEFORE all the others.
%   'ensureFileSep' Add, if not present, a file separator at the end of
%                   each NON-EMPTY string in PAIRS
function nameList=depPackageSubstitutions(nameList,pairs,varargin)
flagRemoveRoot=false;
flagEnsureFileSep=false;

if size(pairs,2)~=2
    error('The argument pairs must be a cell array with two columns')
end

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'removeroot'
            flagRemoveRoot=true;
        case 'ensurefilesep'
            flagEnsureFileSep=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if flagRemoveRoot
    r=depPackageFindCommonRoot(nameList);
    nameList=strrep(nameList,r,'');
end

for iPair=1:size(pairs,1)
    s1=pairs{iPair,1};
    s2=pairs{iPair,2};
    if flagEnsureFileSep
        s1=ensureFileSep(s1);
        s2=ensureFileSep(s2);
    end
    nameList=strrep(nameList,s1,s2);
end

function s=ensureFileSep(s)
if ~isempty(s) && s(end)~=filesep
    s(end+1)=filesep;
end