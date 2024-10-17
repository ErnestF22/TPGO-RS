%Differential of the mapping from individual rotations and translations to the QREM
%function D=essential_fromRT_Diff(R1,T1,R2,T2)
%Computes the matrix D such that
%[va;vb]=D*[vR1;vT1;vR2;vT2],
%where va,vb are the tangent vectors for the two components of the QREM in
%vector coordinates, vR1, vR2 are the tangent vectors in vector coordinates
%for R1 and R2, and vT1,vT2 are the tangent vectors for T1 and T2.
function D=essential_fromRT_Diff(R1,T1,R2,T2,varargin)
flagHorizontalProjection=true;
flagQProvided=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flaghorizontalprojection'
            ivarargin=ivarargin+1;
            flagHorizontalProjection=varargin{ivarargin};
        case 'q'
            ivarargin=ivarargin+1;
            Q=varargin{ivarargin};
            flagQProvided=true;
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

I=eye(3);
Z=zeros(3);

T=T2-T1;

DR0= householderRotation_DiffMat(T,3);
A1= R1'*DR0;
A2= R2'*DR0;

D= [I -A1 Z A1;
    Z -A2 I A2];

if flagHorizontalProjection
    if ~flagQProvided
        Q=essential_fromRT(R1,T1,R2,T2);
    end
    vVert=[Q(3,:)';Q(6,:)']/sqrt(2);
    D=D-vVert*(vVert'*D);
end
