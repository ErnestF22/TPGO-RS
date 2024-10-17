%Add calibration matrices to data
%function data=sfm_addCalibration(data,K)
%Adds field calibration to data structure. If K is omitted or empty, use K
%computed from image coordinates.
%If K is given, its first two rows are multiplied by data.resizeFactor.
function data=sfm_addCalibration(data,K)
flagGivenK=exist('K','var') && ~isempty(K);

if ~flagGivenK
    NImg=size(data.imgSize,2);
    K=zeros(3,3,NImg);
    K(3,3,:)=1;
    for iImg=1:NImg
        halfImgSize=data.imgSize(:,iImg)/2;
        for k=1:2
            K(k,k,iImg)=halfImgSize(k);
            K(k,3,iImg)=halfImgSize(k);
        end
    end
else
    K(1:2,:,:)=data.resizeFactor*K(1:2,:,:);
end
data.calibration=K;
