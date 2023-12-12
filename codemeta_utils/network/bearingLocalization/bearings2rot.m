%Compute relative rotation matrix given bearings w.r.t. third point
%function Rij=bearings2rot(tij,tji,tik,tjk)
%Given two reference frames B_i,B_j and a third point x_k, and the bearings
%of B_j seen in B_i (tij), B_i seen in B_j (tji), x_k seen in B_i (tik) and
%x_k seen in B_j (tjk), compute the rotation to transform coordinates from
%Bj to Bi (Rij)
function Rij=bearings2rot(tij,tji,tik,tjk)

%normals
ni=cnormalize(cross(tij,tik));
nj=cnormalize(cross(-tji,tjk));

%bi-normals
bi=cross(ni,tij);
bj=cross(nj,-tji);

%bearing-aligned axes
Ai=[-tij ni bi];
Aj=[tji nj bj];

%estimated relative rotation
Rij=Ai*Aj';
