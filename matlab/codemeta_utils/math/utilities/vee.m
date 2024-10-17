%p = vee(phat) computes the vee operator (inverse of wedge), that is
%p = vee(hat(p))
%
%See also hat

function p=vee(phat)
warning('This function is deprecated, use vee3 instead')
p = [phat(3,2)-phat(2,3);phat(1,3)-phat(3,1);phat(2,1)-phat(1,2)]/2;
