%DISCRETE_MOTION	
%	Computes the displacement from image correspondences x1 and x2
%	using Ma-Kosecka-Sastry method
%
%	[R,p] = discrete_motion(x1,x2) computes the rotation R and translation p
%	from image correspondences x1 and x2 (both are 3 by n vectors, where n is the 
%	number of correspondences

function [R,p,F] = discrete_motion(x1,x2,robust)

if nargin < 3
    robust = 0;
end

% Estimate essential vector

A = [];
for i=1:size(x1,1)
    A = [A; x1.*repmat(x2(i,:),size(x1,1),1)];
end
A = A';

if rank(A) < 8
   disp('WARNING: Matrix A for estimating the essential matrix E (Ae=0) is') 
   disp('rank deficient. You might want to try using different data. You')
   disp('will possibly have to increase the Field of View or depth variability')
end

if robust == 1
    e = robustnull(A);
else
    [U,S,V] = svd(A);
    e=V(:,9);
end

% Projection onto the normalized essential space 

F=[e(1) e(2) e(3); e(4) e(5) e(6); e(7) e(8) e(9)];

%diag(x2'*F*x1)

[U,S,V] = svd(F);
F = U*diag([S(1,1) S(2,2) 0 ])*V';
S=diag([1 1 0]);
E=U*S*V';

% Displacement recovery from the projected essential matrix
% 4 displacements are recovered, and the one with positive depth is chosen

[U,S,V] = svd_essential(E);
R = V*rot([0 0 1],pi/2)*U';
p = vee3(V*rot([0 0 1],pi/2)*S*V');

%mean(sign((p'*cross(x1,cross(x1,R*x1))))) 
%mean(sign((p'*cross(x1,cross(p*ones(1,size(x2,2)),R*x2))))) 
%if p'*hat(x1(:,1))^2*R*x2(:,1)>0 & -x1(:,1)'*hat(p)^2*R*x2(:,1)>0
if mean(sign((p'*cross(x1,cross(x1,R*x2))))) > 0 &  ...
   mean(sign((p'*cross(x1,cross(p*ones(1,size(x2,2)),R*x2))))) >0
   return
end

R = V*rot([0 0 1],-pi/2)*U';
p = vee3(V*rot([0 0 1],-pi/2)*S*V');

%if p'*hat(x1(:,1))^2*R*x2(:,1)>0 & -x1(:,1)'*hat(p)^2*R*x2(:,1)>0
if mean(sign((p'*cross(x1,cross(x1,R*x2))))) > 0 &  ...
   mean(sign((p'*cross(x1,cross(p*ones(1,size(x2,2)),R*x2))))) >0
   return
end

E=-E;
[U,S,V] = svd_essential(E);

R = V*rot([0 0 1],pi/2)*U';
p = vee3(V*rot([0 0 1],pi/2)*S*V');

%if p'*hat(x1(:,1))^2*R*x2(:,1)>0 & -x1(:,1)'*hat(p)^2*R*x2(:,1)>0
if mean(sign((p'*cross(x1,cross(x1,R*x2))))) > 0 &  ...
   mean(sign((p'*cross(x1,cross(p*ones(1,size(x2,2)),R*x2))))) >0
   return
end

R = V*rot([0 0 1],-pi/2)*U';
p = vee3(V*rot([0 0 1],-pi/2)*S*V');
