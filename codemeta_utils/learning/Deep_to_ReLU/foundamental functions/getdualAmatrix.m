function [As,basic] = getdualAmatrix(z,Ab_set)
Ad = [];
bd = [];
j=1;
for i=1:length(Ab_set)
    if i==1
        Ai = Ab_set(1).A;
        bi = Ab_set(1).b;
        z_size = size(Ai,1);
        zi = z(j:j+z_size-1,:);
        A_cum = zi.*Ai;
        b_cum = zi.*bi;
        j = j+z_size;
    else
        Ai = Ab_set(i).A*A_cum;
        bi = Ab_set(i).A*b_cum+Ab_set(i).b;
        z_size = size(Ai,1);
        zi = z(j:j+z_size-1,:);
        A_cum = zi.*Ai;
        b_cum = zi.*bi;
        j = j+z_size;
    end
    Ad = [Ad; Ai];
    bd = [bd; bi];
end
z(z==1)=-1; % Ad*x+bd>=0 *(-1)
z(z==0)=1; % Ad*x+bd<=0 *(1)
Ad = z.*Ad;
bd = z.*bd;

d = size(Ad,2);
s = size(bd,1); %slack
var = d*2+s;
I = eye(s); % slack variables
As = [Ad -Ad I]; % As(x^+ - x^-)+s=-bs
basic = (var-s+1):var; % basic = slack variables
As(:,end+1) = -bd; % As(x^+ - x^-)+s=-bs
As(end+1,:) = 0;