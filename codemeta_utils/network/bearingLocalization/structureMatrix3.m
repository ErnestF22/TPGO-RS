%function C=structureMatrix3()
%Return a matrix representing the structure of the cross product in R^3.
%C(i,j)=k such that cross(e_i,e_j)=sign(k) e_abs(k) (with the convention
%that the result is zero if k=0).
function C=structureMatrix3()
C=-hat3([1;2;3]);

