clear 
close all

%finding neighbors of z1=[1;1] z2=[1;1] z3=1:
 

z = [1 1 1 0 1];
Z = generate_bit_neighbors(z);


x = -10:0.05:10;

nbColors=8;
colors=parula(nbColors);

A1 = [1 0;0 1];
b1 = -2*[1;1];

A2 = [4 1;3 1];
b2 = -5*[1;1];


A3 = [1 1];
b3 = 1;



for i=1:size(Z,1)
    z1 = Z(i,1:2)';
    z2 = Z(i,3:4)';
    z3 = Z(i,5);
 
    cnt = i;
    
    [A,b]= cumulativeCoeff(A1,b1,A2,b2,A3,b3,z1,z2,z3);
    R2(i,:) = [A b];
    
    Y = funImage(x,x,@(x) network(A1,b1,A2,b2,A3,b3,z1,z2,z3,x),...
        'method','surf','methodOpts',{'FaceColor',...
        colors(mod(cnt-1,nbColors)+1,:)});
    if ~all(all(isnan(Y)))
        legendInfo{cnt} = ['Z = ' num2str([z1',z2',z3])];
%         R2(c,:) = [A b];
%         c=c+1;
    else
        legendInfo{cnt} = [''];
    end
    hold on
end
legend(legendInfo)
xlabel('x')
ylabel('y')
A=[];
A = [R2(:,1:2) -ones(size(R2,1),1) eye(size(R2,1)) -R2(:,3)];
A = [A; zeros(1,size(A,2))];

basic = 4:3+size(R2,1);
[basic,result,A2]=dual_simplex(A,basic);
V = find_vertices(A,basic);

for i=1:size(V,1)
    plot3(V(i,1),V(i,2),V(i,3),'r*')
    hold on
end

% figure(3)
% for i=1:size(R2,1)
% [x,y] = meshgrid(-10:0.2:10); % Generate x and y data
% z = R2(i,1)*x+R2(i,2)*y+R2(i,3);
% surf(x,y,z) %Plot the surface
% hold on
% end
