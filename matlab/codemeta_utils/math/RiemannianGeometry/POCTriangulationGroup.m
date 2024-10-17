function POCTriangulationGroup
switch 2
    case 1
        e3=[0;0;1];
        E=hat(e3);
        X=[randn(2,2); 5+rand(1,2)];
        x1=homogeneous(X,3);
        x2=homogeneous(X-[0;0;2]*[1 1],3);
        xp1=x1(:,1).*x1(:,2);
        xp2=x2(:,1).*x2(:,2);

        disp([x1(:,1)'*E*x2(:,1)  x1(:,2)'*E*x2(:,2) xp1'*E*xp2])

    case 2
        resetRands()
        e=triang_eye();
        v=triang_randTangentNormVector(e);
        x=@(t) triang_exp(e,t*v);

        figure(1)
        fE=@(t) triang_epipolarConstraint(x(t));
        plotfun(fE)
        title('Check that x(t) stays in the manifold')
        
        dx=@(t) x(t).*v;
        
        figure(2)
        check_der(x,dx)
        title('Check that dx(t) is the derivative of x(t)')
        
        figure(3)
        fV=@(t) triang_tangentProjVerticalScalar(x(t),dx(t));
        plotfun(fV)
        title('Check that vertical projection of tangent is zero')

    case 3
        e=triang_eye();
        v=randn(2);
        x=@(t) e+t*v;
        %fE=@(t) triang_epipolarConstraint(x(t));
        %fEx=@(x) [x(1);x(2)].'*[0 -1; 1 0]*[x(3);x(4)];
        fEx=@(x) x(3)*x(2)-x(1)*x(4);
        fE=@(t) fEx(x(t));
        %dfExv=@(x,v) -x(1)*v(4)+x(2)*v(3)+x(3)*v(2)-x(4)*v(1);
        %dfExv=@(x,v) [-x(4); x(3); x(2); -x(1)]'*v(:);
        dfExv=@(x,v) triang_tangentBasisOrth(x,'flagNormalize',false)'*v(:);
        
        dfE=@(t) dfExv(x(t),v);
        check_der(fE,dfE)
        
        
        
end

