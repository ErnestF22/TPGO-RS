function [Rw,Rz,R]=POCReviewRwz(w1,w2,z)
Rz=[cos(z) sin(z) 0; -sin(z) cos(z) 0; 0 0 1];
Rw=1/(1+w1^2+w2^2)*...
    [1+w1^2-w2^2 2*w1*w2 -2*w2;
    2*w1*w2 1-w1^2+w2^2 2*w1;
    2*w2 -2*w1 1-w1^2-w2^2];
R=Rw*Rz;
