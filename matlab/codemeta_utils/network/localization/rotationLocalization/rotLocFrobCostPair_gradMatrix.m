%Matrix representation of the gradient mapping from Ri,Rj to gradc
%function gradMat=rotLocFrobCostPair_gradMatrix(Ri,Rj,Rij,E)
%Output
%   gradMat     [2d x 2d] matrix with a matrix such that 
%               [gradc(:,:,1)',gradc(:,:,2)']==gradMat*[Ri';Rj']
function gradMat=rotLocFrobCostPair_gradMatrix(Ri,Rj,Rij,E)
if ~exist('E','var')
    E=Rij-Ri'*Rj;
end

d=size(Ri,1);

gradMat=zeros(2*d);
gradMat(1:d,d+(1:d))=-E;
gradMat(d+(1:d),1:d)=-E';
