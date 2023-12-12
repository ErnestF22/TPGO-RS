function rot3_tangentPlotFun_test

exp=1;
switch exp
    case 1
        E=eye(3);
        vRId=hat([0;0;1]);%lie_randTangentNormVector(rot_funs,eye(3));

        %S=rot_exp(eye(3),lie_randTangentNormVector(rot_funs,eye(3)));
        S=eye(3);
        %S=rot_exp(eye(3),0.99*pi*v0);
        for iv0=1:2
            v0=hat(E(:,iv0));
            for s=linspace(-pi+1e-6,pi-1e-6,7)
                R0=rot_exp(eye(3),s*v0);
                v=rot_parallel(eye(3),s*v0,vRId);
                Rt=@(t) rot_exp(R0,t*v);
                vt=@(t) rot_parallel(R0,t*v,v);

                rot3_tangentPlotFun(S,Rt,'angle','axes','style','b-','vf',vt,'vfstep',30)
                hold on
            end
        end
        hold off
        axis equal
        set(gcf,'renderer','opengl')
    case 2
        S=eye(3);
        N=[7,7,11];
        R=genRotGrid(N);
        
        r=rot_vee(S,rot_log(S,R));
        plot3(r(1,:),r(2,:),r(3,:),'.')
        hold off
        axis equal
        set(gcf,'renderer','opengl')
end

function R=genRotGrid(N)
epsilon=1e-6;
%N(1), alpha - distance base rotation from origin rotation
%N(2), beta - angle in the plane perpendicular to direction to the origin
%N(3), gamma - number of points on the geodesic

R=[];
v0=hat([1;0;0]);
for alpha=linspace(-pi+epsilon,pi-epsilon,N(1))
    R0=rot_exp(eye(3),alpha*v0);
    
    for beta=linspace(0,pi,N(2))
        vRId=cos(beta)*hat([0;1;0])+sin(beta)*hat([0;0;1]);
        
        v=rot_parallel(eye(3),alpha*v0,vRId);        
        Rt=@(t) rot_exp(R0,t*v);
        for gamma=linspace(-pi,pi,N(3))
            R=cat(3,R,Rt(gamma));
        end
    end
end

