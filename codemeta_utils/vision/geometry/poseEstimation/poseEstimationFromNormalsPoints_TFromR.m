function T=poseEstimationFromNormalsPoints_TFromR(R,nx,M,d)
T=-M\(nx'*R(4:9)'-d');