function POCsimilarityTransformEstimation
resetRands()
x1=rand(2,1);
x2=rand(2,1);
cvx_begin
    variable A(2,2)
    variable k
    e1=A*x1-[1;0];
    e2=A*x2-k*x2;
    minimize (e1'*e1+e2'*e2)
cvx_end

M=[
    x1(1) 0 x1(2) 0 0
    0 x1(1) 0 x1(2) 0
    x2(1) 0 x2(2) 0 -x2(1)
    0 x2(1) 0 x2(2) -x2(2)
];
    
m=[1;0;0;0];

%M*[A(:);k]-m
ak=M\m;
Ab=reshape(ak(1:4),2,2);

keyboard