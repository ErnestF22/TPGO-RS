% ##########################################################################
%% G L O B A L L Y   O P T I M A L   C O N S E N S U S   M A X I M I S A T I O N
%% This package contains the source code which implements optimal Consensus 
% Maximisation proposed in
% T.J. Chin, P. Purkait, A. Eriksson and D. Suter
% Efficient Globally Optimal Consensus Maximisation with Tree Search, 
% In Proceedings of the IEEE Conference on Computer Vision and Pattern 
% Recognition (CVPR), June 2015, Boston
% 
% Copyright (c) 2015 Pulak Purkait (pulak.purkait@adelaide.edu.au.)
% School of Computer Science, The University of Adelaide, Australia
% The Australian Center for Visual Technologies
% http://www.cs.adelaide.edu.au/directory/pulak.purkait
%% Please acknowledge the authors by citing the above paper in any academic 
%  publications that have made use of this package or part of it.
% ##########################################################################

clear; 
close all; 

imargb = im2double(imread('keble_a.jpg'));
imcrgb = im2double(imread('keble_b.jpg'));

matches.im1 = imargb; 
matches.im2 = imcrgb; 

[fa,da] = vl_sift(im2single(rgb2gray(imargb))) ; 
[fc,dc] = vl_sift(im2single(rgb2gray(imcrgb))) ;

xa = fa(1, :); 
ya = fa(2, :); 
xc = fc(1, :); 
yc = fc(2, :); 

% show detected points
% figure(1); clf;
% imagesc(imargb); hold on;
% plot(xa,ya,'+y');
% set(gca, 'XTick', [], 'YTick', []);
% hold off; 

% show all points
% figure(2); clf;
% imagesc(imcrgb); axis image; hold on;
% plot(xc,yc,'+y');
% set(gca, 'XTick', [], 'YTick', []);
% hold off; 

% %
% Compute tentative matches between image 1 (a) and 2 (b) by matching local features
% %
 
[match, scores] = vl_ubcmatch(da, dc, 2.222); 
xat       = xa(match(1, :));
yat       = ya(match(1, :));
xct       = xc(match(2, :));
yct       = yc(match(2, :));

% show all tentative matches
% figure(2); clf;
% imagesc(imargb); hold on;
% plot(xat,yat,'+g');
% hl = line([xat; xct],[yat; yct],'color','y');
% title('Tentative correspondences');
% axis off;

% Pad data with homogeneous scale factor of 1
matches.X1 = [[xat; yat]; ones(1,numel(xat))];
matches.X2 = [[xct; yct]; ones(1,numel(xct))];        

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robustly fit homography
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify the inlier threshold (in noramlized image co-ordinates)

X = unique([matches.X1', matches.X2'], 'rows');
matches.X1 = X(:, 1:3)'; 
matches.X2 = X(:, 4:6)'; 
 
nbpoints = size(matches.X1, 2); 

X = [ matches.X1; matches.X2 ]; 

figure(3), clf;
plot_match(matches, X, [1:nbpoints], 1);
title('Tentative correspondences');
hold off; 

[x1, T1] = normalise2dpts(matches.X1);
[x2, T2] = normalise2dpts(matches.X2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%
% Robustly fit homography
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify the inlier threshold (in noramlized image co-ordinates)
th        = 0.1; 
itrn = 1; 
seed = 0; 
rng(seed); 

trialavg = 0; 
vavrg = 0; 
tic; 
for i=1:itrn
    [H1, inliers1, trial] = ransacfithomography(x1, x2, th);
    vavrg = vavrg + numel(inliers1); 
    trialavg = trialavg + trial; 
end
t1 = toc/itrn; 
vavrg =nbpoints-vavrg/itrn; 
trialavg = trialavg/itrn; 
v1 = 1:nbpoints; 
v1(inliers1) = []; 

figure(4), clf;  
plot_match(matches, X, inliers1, 1); 
title('RANSAC Inliers');
hold off; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%
% ASTAR fit homography
%%%%%%%%%%%%%%%%%%%%%%%%%%%

u1 = [x1; zeros(6, nbpoints); x1];
u2 = x2(1:2, :); 
u2 = repmat(u2(:), 1, 3).*reshape([x1;x1], 3, 2*nbpoints)'; 
A = [reshape(u1, [6, 2*nbpoints])', -u2]; 
b = zeros(2, nbpoints); 
c = [zeros(6, nbpoints); x1]; 
d = zeros(1, nbpoints); 
 

tic
[P5, Iopt5, v5, nnum5, xnum5] = maxconASTAR(A, b, c, d, th); 
t5 = toc; 
H5 = reshape(P5(:),3,3)';

inliers5 = 1:nbpoints; 
inliers5(v5) = []; 
figure(5), clf; 
plot_match(matches,X,inliers5, 1); 
title('ASTAR Inliers');
hold off; 

fprintf(' RANSAC time %f, ASTAR time %f\n',t1, t5);
disp('Average Inliers:'); 
disp(['RANSAC = ' num2str(nbpoints - vavrg)]); disp(['ASTAR  = ' num2str(nbpoints - numel(v5))]);  


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warp and composite images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HR = T2\H1*T1;
HA = T2\H5*T1;

% figure(5); clf;
bbox=[-300 750 -100 600];    % image space for mosaic
% warp image b to mosaic image using an identity homogrpahy
% Image b is chosen as the reference frame
iwc = vgg_warp_H(matches.im2, eye(3), 'linear', bbox);
% imshow(iwc); axis image;
% warp image 1 to the reference mosaic frame (image 2) 
figure(6); clf;
iwa = vgg_warp_H(matches.im1, HR, 'linear', bbox);  % warp image a to the mosaic image
imshow(iwa); axis image;
imagesc(double(max(iwc,iwa))); % combine images into a common mosaic (take maximum value of the two images)
set(gca, 'XTick', [], 'YTick', []);
title('RANSAC fit Homography');

% warp image 1 to the reference mosaic frame (image 2) 
figure(7); clf;
iwa = vgg_warp_H(matches.im1, HA, 'linear', bbox);  % warp image a to the mosaic image
imshow(iwa); axis image;
imagesc(double(max(iwc,iwa))); % combine images into a common mosaic (take maximum value of the two images)
set(gca, 'XTick', [], 'YTick', []);
title('ASTAR fit Homography');

% 18746 sec 48092
