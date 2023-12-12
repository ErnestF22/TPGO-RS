function POCComputeCurvatureBasisTSO3
bracket=@(R,A,B) R*(R'*A*R'*B-R'*B*R'*A);
R=rot_randn(eye(3));
I=eye(3);
for d1=1:3
    Ed1=R*hat(I(:,d1));
    for d2=1:3
        Ed2=R*hat(I(:,d2));
        for d3=1:3
            Ed3=R*hat(I(:,d3));
            r=trace(bracket(R,bracket(R,Ed2,Ed1),Ed3)'*Ed1);
            if r<1e-15
                r=0;
            end
            if abs(r-2)<1e-15
                r=2;
            end
            fprintf('(%d,%d,%d): %d\n',d1,d2,d3,r)
        end
    end
end
