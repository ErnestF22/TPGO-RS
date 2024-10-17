function POCReviewR_test
%%
syms w1 w2 z
%w1=0.5*rand;
%w2=0.5*rand;
%z=0;
Rw=POCReviewRwz(w1,w2,z);
[w1b,w2b,a,b,c]=POCReviewRzTow(Rw);
%Rwb=POCReviewRwz(w1b,w2b,0);
%simplify(Rw.'*Rw)
simplify([w1b-w1 w2b-w2])
%%
%keyboard
