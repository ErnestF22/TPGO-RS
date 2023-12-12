function [Arf, Pnew] = find_tableau_after_flip_P(z_before,bit,P,Ar,Ab_set)
% z_before is z before flipping, bit is the flipped bit
% Ar is the final tableau for z_before
% Arf returns the final tableau after flipping
% Pnew returns the updated P matrix

z_new = z_before;
z_new(bit) = ~z_new(bit); % z after flipping
b1 = getbcolumn(z_new,Ab_set); % b column after flipping
b = getbcolumn(z_before,Ab_set); % b column before flipping
d = size(Ar,2)-size(Ar,1);
if z_before(bit,:)==1
    m2 = (b1-b)./(-b1(bit)); % bit column of M2
else
    m2 = (b1-b)./b(bit,end);
end
Arf = Ar;
Arf(:,bit+d)=Arf(:,bit+d)+P*m2;

M = eye(size(P,2)); % I
m1 =(b-b1)/b1(bit); % bit column of M1
m1(1:bit-1,:) = 0;
M(:,bit) = M(:,bit)+m1; % = I+M1
Pnew = P*M; % = P(I+M1)
Pnew = round(Pnew,15);
% [A1,~] = getdualAmatrix(z_new,Ab_set); % initial tableau after flipping
% [A3,~] = getdualAmatrix(z_before,Ab_set); % initial tableau before flipping
% A1 = M*A1;
% A1(:,bit+d) = [];
% A3(:,bit+d) = [];
% if abs(sum((A3-A1),'all'))>=1e-3
%     disp('conversion fail!'); % all the columns are same except for column bit
% end
% 
% [A1,~] = getdualAmatrix(z_new,Ab_set); % initial tableau after flipping
% sumall = abs(sum((Pnew*A1-Arf),'all'));
% if sumall>=1e-3
%     disp('flip error!');
% end
end