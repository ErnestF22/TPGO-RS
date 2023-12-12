%an arbitrary graph
A=full(bucky); 
A=A(1:30,1:30); %this is to make things a little more interesting

%an arbitrary number of clusters
k=3;

%node degrees
d=sum(A,2);
%various versions of the normalization matrix
D=diag(d);      
Dsqrt=diag(sqrt(d));
DsqrtInv=diag(1./sqrt(d));

%the Laplacian and its normalized version
L=D-A;
Lnorm=DsqrtInv*L*DsqrtInv;

%find the solution with the normalized Laplacian (problem **2 in the notes)
[Vnorm Enorm]=eig(Lnorm);
Y=Vnorm(:,1:k);

%get the approximate solution with the original Laplacian (problem *2 in
%the notes)
Z=DsqrtInv*Y;

%check the constraints on Z (problem *1 in the notes)
disp('Check constraints on continuous solution Z')
disp('Z''*D*Z=')
disp(Z'*D*Z)
disp('Z*inv(Z''*Z)*Z''*ones(30,1)=')
disp(Z*inv(Z'*Z)*Z'*ones(30,1))
