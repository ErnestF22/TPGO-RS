%Construct essential matrix from rotations and translations
%function epipolarBuildEFromRT(R1,T1,R2,T2,varargin)
%Input
%   R1,T1,R2,T2     Rotation and translation of the cameras
%                   If R2,T2 are omitted or empty, interpret R1,T1 as
%                   relative rotation and translation
%Optional input
%   'references'    Interpret R1,T1,R2,T2 as references (default).*
%   'poses'         Interpret R1,T1,R2,T2 as poses.*
%   'methodAbsolutePoses', method   Method for interpreting R1,T1,R2,T2.*
%                                   Can be 'reference' or 'pose'.
%   * These additional inputs have no effect if R2,T2 are omitted, i.e., if
%     R1,T1 represent relative poses
function E=epipolarBuildEFromRT(R1,T1,R2,T2,varargin)
methodAbsolutePoses='reference';
flagRelativePoses=false;

if ~exist('R2','var') || isempty(R2)
    flagRelativePoses=true;
end

NPoses=size(R1,3);
if NPoses>1
    E=zeros(3,3,NPoses);
    for iPoses=1:NPoses
        if flagRelativePoses
            E(:,:,iPoses)=epipolarBuildEFromRT(R1(:,:,iPoses),T1(:,iPoses),[],[],varargin{:});
        else
            E(:,:,iPoses)=epipolarBuildEFromRT(R1(:,:,iPoses),T1(:,iPoses),R2(:,:,iPoses),T2(:,iPoses),varargin{:});
        end
    end
else

    %optional parameters
    ivarargin=1;
    while(ivarargin<=length(varargin))
        switch(lower(varargin{ivarargin}))
            case 'references'
                methodAbsolutePoses='reference';
            case 'poses'
                methodAbsolutePoses='pose';
            case 'methodabsoluteposes'
                ivarargin=ivarargin+1;
                methodAbsolutePoses=lower(varargin{ivarargin});
            otherwise
                    error(['Argument ' varargin{ivarargin} ' not valid!'])
        end
        ivarargin=ivarargin+1;
    end

    if ~flagRelativePoses
        if strcmpi(methodAbsolutePoses,'pose')
            %if we are in the pose interpretation, pass to the reference
            %interpretation

            [R1,T1]=invRT(R1,T1);
            [R2,T2]=invRT(R2,T2);
        end
        E=-R1'*hat3(cnormalize(T2-T1))*R2;
    else
        E=-hat3(cnormalize(T1))*R1;
    end
end
