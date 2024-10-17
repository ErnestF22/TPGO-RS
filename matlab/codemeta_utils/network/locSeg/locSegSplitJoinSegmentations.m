%Merge segmentations
%function [members,attributes]=locSegSplitJoinSegmentations(members,newMembers,attributes,newAttributes)
%If a new segmentation is the same as an old one, add the corresponding new
%attributes to the existing ones. If the old or new segmentations are a
%subset of the other, split the largest one, and assign the attributes
%correspondingly. If segmentations are not a subset of each other but
%overlap, throw an error.
%Inputs
%   members, newMembers         [NSegments x NNodes] boolean matrix with
%                               segmentation memberships
%   attributes, newAttributes   [NSegments x 1] cell array with attributes
%                               for each segment
function [members,attributes]=locSegSplitJoinSegmentations(members,newMembers,attributes,newAttributes)
flagUseAttributes=true;
if nargin<=2
    flagUseAttributes=false;
end

[NNewMembers,NNodes]=size(newMembers);
for iMember=1:NNewMembers
    n=newMembers(iMember,:);
    idxSame=findSameSegmentation(members,n);
    if length(idxSame)>1
        error('Detected repeated segment in members')
    elseif length(idxSame)==1
        %n is the same as an existing group
        if flagUseAttributes
            attributes{idxSame}=[attributes{idxSame} newAttributes{iMember}];
        end
    else
        %n is not the same as an existing group
        idxSupersets=findSupersetSegmentation(members,n);
        NIdxSupersets=length(idxSupersets);
        if NIdxSupersets~=0
            %n is a subset of one or more existing groups
            splitMembers=zeros(2*NIdxSupersets,NNodes);
            if flagUseAttributes
                splitAttributes=cell(1,2*NIdxSupersets);
            end
            for iIdxSupersets=1:NIdxSupersets
                idx=idxSupersets(iIdxSupersets);
                m=members(idx,:);
                m1=splitSegmentation(m,n);
                splitIdx=(iIdxSupersets-1)*2+1;
                splitMembers(splitIdx,:)=n;
                splitMembers(splitIdx+1,:)=m1;
                if flagUseAttributes
                    splitAttributes{splitIdx}=[attributes{idx} newAttributes{iMember}];
                    splitAttributes{splitIdx+1}=attributes{idx};
                end
            end
            idxNotSupersets=setdiff(1:size(members,1),idxSupersets);
            members=[members(idxNotSupersets,:); splitMembers];
            if flagUseAttributes
                attributes=[attributes(idxNotSupersets) splitAttributes];
            end
        else
            %n is not a subset of an existing group
            idxSubsets=findSubsetSegmentation(members,n);
            NIdxSubsets=length(idxSubsets);
            if NIdxSubsets~=0
                %n is a subset of one or more existing groups
                for iIdxSubsets=1:NIdxSubsets
                    idx=idxSubsets(iIdxSubsets);
                    m=members(idx,:);
                    n=splitSegmentation(n,m);
                    if flagUseAttributes
                        attributes{idx}=[attributes{idx} newAttributes{iMember}];
                    end
                end
                if any(n)
                    error('One of the rows in newMembers is incompatible with members')
                end
            else
                error('One of the rows in newMembers is incompatible with members')
            end
        end
    end
end

%return index(es) of rows of members which the same as n
function idx=findSameSegmentation(members,n)
NMembers=size(members,1);
idx=find(all(members==ones(NMembers,1)*n,2));

%return index(es) of rows of members which are supersets of n
function idx=findSupersetSegmentation(members,n)
NMembers=size(members,1);
idx=find(all(or(members==1, and(members==0,ones(NMembers,1)*n==0)),2));

%return index(es) of rows of members which are subsets of n
function idx=findSubsetSegmentation(members,n)
NMembers=size(members,1);
nRep=ones(NMembers,1)*n;
idx=find(all(or(nRep==1, and(members==0,nRep==0)),2));

%Given n which represent a subset of m, give m1, which contains
%elements only present in m
function m1=splitSegmentation(m,n)
m1=and(m==1,n~=1);
