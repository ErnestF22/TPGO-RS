function POCReviewProbabilityUpdate
N=20;
M=10;
a=0.1;
p=zeros(M,N);
p(:,1)=linspace(0,1,M)';
for n=2:N
    p(:,n)=0.3*(p(:,n-1)-0.5)+0.5;
end

plot(p')

