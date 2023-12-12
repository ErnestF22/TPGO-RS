%Compute the matrix representation of differential of sphere_log for the 2-D sphere
%function D=sphere3_logDiffMat(y1,y2)
%Optional arguments
%   'project'   Add projection matrix on tangent space of y1
function D=sphere3_logDiffMat(y1,y2,varargin)
flagProject=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'project'
            flagProject=true;
    end
    ivarargin=ivarargin+1;
end

%rotation of pi around y2
Ry2pi=-eye(3)+2*(y2*y2');
%householder rotation aligning y1 to y2
H=householderRotation(y1,y2);
DH=householderRotation_DiffMat(y1,y2);
%rotation aligning y1 to y2 and closest to the origin
Rmin=H*Ry2pi;
DRmin=Ry2pi'*DH-2*[zeros(3) hat(y2)];
%logarithm of Rmin
l=logrot(Rmin);
Dl=rot3_logDiff(eye(3),l,'tangentVec');

%put everything together in the differential
%Note: sphere_log(y1,y2)=hat(y1)*l(t)
%D=[-hat(l) zeros(3)]+hat(y1)*Dl*DRmin;
D=[hat(l) zeros(3)]+hat(y1)*Dl*DRmin;

if flagProject
    D=orthComplementProjector(y1)*D;
end
