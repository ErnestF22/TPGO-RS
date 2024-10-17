%Compute relative rotation from two absolute rotations
%function R12=relativeRotations(R1,R2,varargin)
%The rotations R12 transforms the data in the reference frame of R2 to
%the reference frame of R1
%Optional inputs
%   '21'    give the transformation from the frame of R2 to the one of R1
%           (default)
%   '12'    give the transformation from the frame of R1 to the one of R2
%   'edges' Interpret R2 as a [NEdges x 2] list of edges. The rotations R12
%       will contain all the relative rotations obtained from the pairs of
%       slices of R1 indicated in this edge list
function R12=relativeRotations_codemeta(R1,R2,varargin)
methodAbsolutePoses='reference';
flagInvertR=false;
flagG2IsEdges=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case '12'
            flagInvertR=true;
        case '21'
            flagInvertR=false;
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        case 'edges'
            flagG2IsEdges=true;
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

N1=size(R1,3);
N2=size(R2,3);

if flagG2IsEdges
    E=R2;
    R2=R1;
else
    E=[vec(repmat(1:N1,N2,1)) vec(repmat((1:N2)',1,N1))];
end
NEdges=size(E,1);
R12=zeros(3,3,NEdges);
for iEdge=1:NEdges
    i1=E(iEdge,1);
    i2=E(iEdge,2);
    switch methodAbsolutePoses
        case 'reference'
            R12(:,:,iEdge)=invR(R1(:,:,i1))*R2(:,:,i2);
        case 'pose'
            R12(:,:,iEdge)=R1(:,:,i1)*invR(R2(:,:,i2));
        otherwise
            error('Absolute pose method specification not valid')
    end
end
R12=squeeze(R12);
if flagInvertR
    R12=invR(R12);
end
