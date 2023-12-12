function l=POCReviewInequalitySkew()
z=zeros(3);
w=cnormalize(randn(3,1));
hw=hat3(w); 
l=eig([z -hw; hw z]+eye(6));
