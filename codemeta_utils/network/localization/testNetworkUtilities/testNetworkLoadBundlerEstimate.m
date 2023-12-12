function [G,X,br]=testNetworkLoadBundlerEstimate(dirName)

%load bundler results
br=bundlerReadFile(dirName);
cm=br.cameras;

%prepare ground truth poses
G=[cm.R permute(shiftdim(cm.T,-1),[2 1 3]); zeros(1,3,br.NCameras) ones(1,1,br.NCameras)];

%get also 3-D points
X=br.points.position;
