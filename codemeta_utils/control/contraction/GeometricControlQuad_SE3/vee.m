% This function returns a vector that comes from its skew symmetric matrix
function vec = vee(skewSymmetricMat)
    vec = [skewSymmetricMat(3,2); skewSymmetricMat(1,3); skewSymmetricMat(2,1)] ;
end