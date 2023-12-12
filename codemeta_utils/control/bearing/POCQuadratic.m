function POCQuadratic
switch 2
    case 1
        N=4;
        A=randn(N,N);
        c=randn(N,1);
        a=0;
        for i=1:N
            for j=i+1:N
                cij=[c(i); c(j)];
                Aij=[A(i,i)/(N-1) A(i,j); A(j,i) A(j,j)/(N-1)];
                a=a+cij'*Aij*cij;
            end
        end

        [a c'*A*c]
    case 2
        syms b1 b2 g N1
        M=[b1/N1 b2*g; b1*g b2/N1];
        l=eig(M)
        subs(l,[b1 b2 N1 g], [rand rand 2 1])
        
        
end