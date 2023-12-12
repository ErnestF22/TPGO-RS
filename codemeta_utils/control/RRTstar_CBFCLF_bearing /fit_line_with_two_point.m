function [A,b] = fit_line_with_two_point(from_node,to_node)
n = size(to_node,2);
A=[];
b=[];
for i=1:n
    x = to_node(:,i);
    coeff = polyfit([from_node(1) x(1)],[from_node(2) x(2)],1);
    a = [coeff(1) 1];
    A = [A;a];
    b =[b;coeff(2)];
end
end