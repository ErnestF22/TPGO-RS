function rot3_plotvf_test

exp=2;
switch exp
    case 1
        clear variables
        S=eye(3);
        vf=@(R) -rot_log(R,S);
        R=@(t) rot_exp(S,t*hat([0;0;1]));
        rot3_plotvf(S,'5x3x7',vf)
        %myR=cat(3,R(0),R(-pi/2),R(pi/2));
        %myV=cat(3,vf(myR(:,:,1)),vf(myR(:,:,2)),vf(myR(:,:,3)));
        %rot3_plotvf(S,myR,myV)
        axis equal
        set(gcf,'renderer','opengl')
    case 2
        grid='7x3x7';
        figure(1)
        vf=@(R) -rot_log(R,eye(3));
        R=@(t) rot_exp(eye(3),t*hat([0;0;1]));
        rot3_plotvf(eye(3),grid,vf,'axes','points')
        axis equal
        title('Distance')

        figure(2)
        A=hat([1;0;0]);
        vf=@(R) rot_eucl2RiemGrad(R,A);
        rot3_plotvf(eye(3),grid,vf,'axes','points')
        axis equal
        title('trace(A^TR)')
        
        figure(3)
        epsilon=0.8;
        v=hat([1;0;0]);
        vf=@(R) gradMinGeodDist(R,eye(3),v);
        rot3_plotvf(eye(3),grid,vf,'points','axes','intervals',linspace(-pi+epsilon,pi-epsilon,7),linspace(0,pi,3),linspace(-pi+epsilon,pi-epsilon,7))

end

function g=gradMinGeodDist(S,R0,v)
[m,g]=minGeodDist(S,R0,v,'grad');


        
