%Create transformation matrix to create data with given covariance matrix
%function SigmaHalf=sigmaToTransformationMatrix(sigmaNoise)
%Inputs
%   sigmaNoise  If a scalar, the standard deviation of the distribution
%               (i.e., the squared root of the variance). The resulting
%               transformation will produce isotropically distributed data.
%               If a matrix, the resulting transformation is the squared
%               root of the pseudo inverse.
%   D           Dimension of the final data. If omitted, D is given by the
%               dimension of sigmaNoise 
function SigmaHalf=sigmaToTransformationMatrix(sigmaNoise,D)
if ~exist('D','var')
    D=size(sigmaNoise,1);
end
N=size(sigmaNoise,3);
SigmaHalf=zeros(D,D,N);
for iN=1:N
    if size(sigmaNoise,1)==1
        SigmaHalf(:,:,iN)=sigmaNoise(:,:,iN)*eye(D);
    else
        SigmaHalf(:,:,iN)=invSqrtPSDMatrix(sigmaNoise(:,:,iN));
    end
end