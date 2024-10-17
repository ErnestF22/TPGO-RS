clear all
close all
clc

addpath('subject15/')
addpath('util/')


%% frame selection from sequence

t =   1:1:150;

%% generate noisy data  

sigma = 0.25; % std of noise (isotropic zero mean normal)

[X,Xn,E] =  CreateData(sigma) ;

Xt = X(:,:,t);
Xnt = Xn(:,:,t);

nE = size(E,1); % number of links
nV = nE + 1; % number of joints,  (tree assumption)

 
%% estimate the lengths from groundtruth
lens = zeros(nE,1);
for ie  = 1:nE , 
    lens(ie)   = norm(Xt(:,E(ie,1),1)-Xt(:,E(ie,2),1)); 
end

sigma_lens = 10^(-8)*eye(nE); % uncertainty of lengths, lengths are in the state vector

%% measurements


y = reshape( Xn(:,:,t),[45 numel(t)]); % projections (calibrated coordinates)
sigma_y = sigma^2*eye(3*nV); % projection error, assumed isotropic, covariance





%% initialization 


% initialization
sigma_X0 = sigma^2*eye(3*nV); %
X0 =  Xn(:,:,t(1));


%% estimation here
% Xe: updates, Xp: predictions
[Xe,Xp] = IcraEstimation(y,sigma_y,E,lens,sigma_lens,X0,sigma_X0,t);


err = zeros(numel(t),1);
errinit = zeros(numel(t),1);
for i=1:numel(t)
    err(i) = norm(Xe(:,:,i)-Xt(:,:,i),'fro');
    err_init(i) = norm(Xnt(:,:,i)-Xt(:,:,i),'fro');
end
figure(20),
plot(err,'r'); 
hold on;
plot(err_init,'b');
hold off;
legend('REKF','Initial')
pause(1)






%% display results
figure(3),
for f = 1 :numel(t)
    
    plot3(squeeze(Xt(1,:,f)),squeeze(Xt(2,:,f)),squeeze(Xt(3,:,f)),'bo')
    title(['frame #' sprintf('%d',f)])
    axis equal
  

    hold on;
    
    for ie = 1:nE
    plot3(squeeze(Xt(1,E(ie,:),f)),squeeze(Xt(2,E(ie,:),f)),squeeze(Xt(3,E(ie,:),f)),'b-')
    end
    
    
    plot3(squeeze(Xe(1,:,f)),squeeze(Xe(2,:,f)),squeeze(Xe(3,:,f)),'r+')
    
     for ie = 1:nE
    plot3(squeeze(Xe(1,E(ie,:),f)),squeeze(Xe(2,E(ie,:),f)),squeeze(Xe(3,E(ie,:),f)),'r-')
     end
    
    pause(0.01)
     
    hold off;
   

end