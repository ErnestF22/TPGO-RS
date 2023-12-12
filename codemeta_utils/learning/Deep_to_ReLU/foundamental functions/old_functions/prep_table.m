function [As,basic] = prepTabla(A,d)
for i = 1:size(A,1)
    Aa(i,:) = A(i,1:d);
    ba(i,1) = A(i,d+1);
end
As = [Aa -Aa eye(size(Aa,1)) ba];
As(end+1,:) = 0;
basic = [2*d+1:size(As,2)-1];
end