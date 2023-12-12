function [Arf] = find_tableau_after_flip(z_b,bit,P,As,Ar,Ab_set)
% z_b is z before flipping, bit is the flipped bit, As is initial tableau for
% z_b, Ar is the final tableau for z_b
% Arf returns the final tableau after flipping
z_b(bit,:)=~z_b(bit,:);
[bcolumn,d] = getbcolumn(z_b,Ab_set); % get the b column
d = 2*d;
m1 = (bcolumn-As(:,end));
if z_b(bit,:)==1
    v = m1./(-bcolumn(bit,end));
else
    v = m1./As(bit,end);
end
% if z_b(bit,:)==1
%     v = m1./As(bit,end);
% else
%     v = m1./(-bcolumn(bit,end));
% end
v(end) = 0;
Arf = Ar;
Arf(:,bit+d)=Arf(:,bit+d)+P*v;
end