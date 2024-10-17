%Build the Jacobian matrix for different kinds of constraints
%function Jr=locSegBuildConstraintJacobian(constraintType,xi,xj,varargin)
% Typical inputs:
%   'fixedBearings',xi,xj
%   'relativeBearings',xi,xj,Ri,tij
%   'distances',xi,xj
%   'fixedRotations',xi,~
%   'relativeRotations',xi,xj,Ri,Rj,Rij
%   'fixedTranslations',xi,xj
%   'relativeTranslations',xi,xj,Ri,Tij
%   'fixedVector',xi,~,Ri,g
%   'fixedPlane',~,~,n
%   'fixedLine',~,~,n
%   'fixedPosition',~,~
%Returns a [K x 2D] matrix with the Jacobian of the constraint evaluated at
%the given state, where K is the dimension of the constraint and D is the
%dimension of the coordinates at each node (D=3 for 2-D, D=6 for 3-D)
function Jr=locSegJacobianConstraintsBuild(constraintType,xi,xj,varargin)
switch constraintType
    case 'fixedBearings'
        Jr=jacobianFixedBearingsConstraints(xi,xj);
    case 'distances'
        Jr=jacobianDistanceConstraints(xi,xj);
    case 'relativeBearings'
        Ri=varargin{1};
        tij=varargin{2};
        Jr=jacobianRelativeBearingsConstraints(xi,xj,Ri,tij);
    case 'fixedRotations'
        Jr=jacobianFixedRotations(xi);
    case 'relativeRotations'
        Ri=varargin{1};
        Rj=varargin{2};
        Rij=varargin{3};
        Jr=jacobianRelativeRotationsConstraints(Ri,Rj,Rij);
    case 'fixedTranslations'
        dimData=size(xi,1);
        Jr=jacobianFixedTranslationsConstraints(dimData);
    case 'relativeTranslations'
        Ri=varargin{1};
        Tij=varargin{2};
        Jr=jacobianRelativeTranslationsConstraints(xi,xj,Ri,Tij);
    case 'fixedVector'
        Ri=varargin{1};
        g=varargin{2};
        Jr=jacobianFixedVector(xi,Ri,g);
    case 'fixedPlane'
        n=varargin{1};
        Jr=jacobianFixedPlane(xi,n);
    case 'fixedLine'
        n=varargin{1};
        Jr=jacobianFixedLine(xi,n);
    case 'fixedPosition'
        dimData=size(xi,1);
        Jr=jacobianFixedPosition(dimData);
    otherwise
        error('typeConstraints not valid');
end

function Jr=jacobianFixedPosition(dimData)
JRotation=singleJacobianNullRotation(dimData,dimData);
JTranslation=eye(dimData);
J=[JTranslation JRotation];
Jr=[J zeros(size(J))];

function Jr=jacobianFixedLine(n)
dimData=size(Ri,1);
JRotation=singleJacobianNullRotation(dimData-1,dimData);
JTranslation=computeOrthogonalComplement(n)';
J=[JTranslation JRotation];
Jr=[J zeros(size(J))];

function Jr=jacobianFixedPlane(n)
dimData=size(Ri,1);
JRotation=singleJacobianNullRotation(1,dimData);
JTranslation=n';
J=[JTranslation JRotation];
Jr=[J zeros(size(J))];

function Jr=jacobianFixedVector(Ri,g)
dimData=size(Ri,1);
S=[0 -1; 1 0];
switch dimData
    case 2
        JRotation=Ri*S*g;
    case 3
        JRotation=-Ri*hat(g);
end
JTranslation=singleJacobianNullTranslationFromRotation(dimData);
J=[JTranslation JRotation];
Jr=[J zeros(size(J))];

function Jr=jacobianFixedTranslationsConstraints(dimData)
[~,dimRotations]=locSegDimData2DimCoordinates(dimData);
Jr=[eye(dimData) zeros(dimData,dimRotations) -eye(dimData) zeros(dimData,dimRotations)];

function Jr=jacobianRelativeRotationsConstraints(Ri,Rj,Rij)
S=[0 -1; 1 0];
dimData=size(Ri,1);
switch dimData
    case 2
        Jr=[zeros(4,2) reshape(Ri*S,[],1) zeros(4,2) -reshape(Rij*Rj*S,[],1)];
    case 3
        Jr=[zeros(9,3) vcthat(Ri) zeros(9,3) -vcthat(Rij*Rj)];
end

function Jr=jacobianFixedRotations(xi)
dimData=length(xi);
JRotation=singleJacobianFixedRotation(dimData);
JTranslationFromRotation=singleJacobianNullTranslationFromRotation(dimData);
J=[JTranslationFromRotation JRotation];
Jr=[J zeros(size(J))];

function Jr=jacobianFixedBearingsConstraints(xi,xj)
dimData=length(xi);
JTranslation=singleJacobianBearingsConstraintsTranslation(xi,xj);
JRotationFromTranslation=singleJacobianNullRotation(dimData,dimData);
Jr=[JTranslation JRotationFromTranslation -JTranslation -JRotationFromTranslation];

function Jr=jacobianDistanceConstraints(xi,xj)
dimData=length(xi);
JTranslation=(xi-xj)'/norm(xi-xj);
JRotation=singleJacobianNullRotation(1,dimData);
Jr=[JTranslation JRotation -JTranslation -JRotation];

function Jr=jacobianRelativeBearingsConstraints(xi,xj,Ri,tij)
dimData=length(xi);
JTranslation=singleJacobianBearingsConstraintsTranslation(xi,xj);
switch dimData
    case 2
        JRotation1=Ri*[0 1; -1 0]*tij;
    case 3
        JRotation1=Ri*hat(tij);
end
JRotation2=singleJacobianNullRotation(dimData,dimData);
Jr=[JTranslation JRotation1 -JTranslation -JRotation2];

function Jr=jacobianRelativeTranslationsConstraints(xi,xj,Ri,Tij)
dimData=length(xi);
JTranslation=eye(dimData);
switch dimData
    case 2
        JRotation1=Ri*[0 1; -1 0]*Tij;
    case 3
        JRotation1=Ri*hat(Tij);
end
JRotation2=singleJacobianNullRotation(dimData,dimData);
Jr=[JTranslation JRotation1 -JTranslation -JRotation2];

function J=singleJacobianBearingsConstraintsTranslation(xi,xj)
dimData=length(xi);
n=norm(xi-xj);
J=(eye(dimData)-(xi-xj)*(xi-xj)'/n^2)/n;

function J=singleJacobianNullRotation(dimConstraint,dimData)
switch dimData
    case 2
        J=zeros(dimConstraint,1);
    case 3
        J=zeros(dimConstraint,3);
end

function J=singleJacobianNullTranslationFromRotation(dimData)
switch dimData
    case 2
        J=zeros(1,dimData);
    case 3
        J=zeros(3,dimData);
end

function JRotation=singleJacobianFixedRotation(dimData)
switch dimData
    case 2
        JRotation=1;
    case 3
        JRotation=eye(3);
end
