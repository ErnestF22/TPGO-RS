%SVD_ESSENTIAL 
%	Computes the singular value decomposition of an essential matrix E.
%	[U,S,V] = svd_essential(E) returns U and V which are rotation matrices 
%	that satisfy det(U)=det(V)=1, and S = diag([lambda lambda 0]).
%
%	See also svd

function [U,S,V] = svd_essential(E)

[U,S,V] = svd(E);
du=det(U);
dv=det(V);

if du<0 & dv<0
   U=-U; V=-V;
   return
elseif du<0 & dv>0
   S1=rot([0 0 1],pi/2)*S;
   U=-U*expm(S1*pi)*rot([0 0 1],pi/2);
   V=V*rot([0 0 1],pi/2);
	return   
elseif du>0 & dv<0
   S1=rot([0 0 1],pi/2)*S;
   U=U*expm(S1*pi)*rot([0 0 1],pi/2);
   V=-V*rot([0 0 1],pi/2);
   return
end