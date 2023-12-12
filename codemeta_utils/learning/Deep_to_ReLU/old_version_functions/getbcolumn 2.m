function [bcolumn,d] = getbcolumn(z,Ab_set)
bd = [];
j=1;
for i=1:length(Ab_set)
    if i==1
        bi = Ab_set(1).b;
        z_size = size(bi,1);
        zi = z(j:j+z_size-1,:);
        b_cum = zi.*bi;
        j = j+z_size;
    else
        bi = Ab_set(i).A*b_cum+Ab_set(i).b;
        z_size = size(bi,1);
        zi = z(j:j+z_size-1,:);
        b_cum = zi.*bi;
        j = j+z_size;
    end
    bd = [bd; bi];
end
z(z==1)=-1; % Ad*x+bd>=0 *(-1)
z(z==0)=1; % Ad*x+bd<=0 *(1)
bd = z.*bd;
bd(end+1) = 0;
bcolumn = -bd;
d = size(Ab_set(1).A,2);