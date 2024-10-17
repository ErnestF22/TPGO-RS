clear all
close all
clc

%% load data

load joints15
load 15_02

Sl = S(:,mocapIndices);

E = [];
for i=1:numel(skel.tree)
    for j=1:numel(skel.tree(i).children)
       E = [E;i skel.tree(i).children(j) ] ;
    end  
end



%% change depth to make it approximately constant

F= size(Sl,1)/3;
P = size(Sl,2);


Zdes = 10; % desired depth

AnglePerFrame = 5;

Angle = (1:F)*(AnglePerFrame/180*pi);

for i=1:F
    
    meanX = mean(Sl(3*i-2,:));
    meanY = mean(Sl(3*i-1,:));
    meanZ = mean(Sl(3*i,:));
    
    Sl(3*i,:) =  Sl(3*i,:) -meanZ;
     Sl(3*i-2,:) =  Sl(3*i-2,:) -meanX;
    Sl(3*i-1,:) =  Sl(3*i-1,:) -meanY;
    
    
    % rotate properly here
    
    rotaxis = [0 1 0]';
    Ry = expm(hatop(rotaxis)*Angle(i));
    Sl(3*i-2:3*i,:) = Ry*Sl(3*i-2:3*i,:);
    
    
    % translate to the desired depth
      Sl(3*i,:) =  Sl(3*i,:) + Zdes;
end








%% plot

a=[0 0 ];
figure,plot3(Sl(1,:),Sl(2,:),Sl(3,:),'bo');view(a)
return;

if 0
        for i=1200:1500

        plot3(Sl(3*i-2,:),Sl(3*i,:),Sl(3*i-1,:),'bo');view(a)
        hold on;
        plot3(Sl(3*i-2,[2 3 4]),Sl(3*i,[2 3 4]),Sl(3*i-1,[2 3 4]),'r+');
        plot3(Sl(3*i-2,[5 6 7]),Sl(3*i,[5 6 7]),Sl(3*i-1,[5 6 7]),'ro');
        plot3(Sl(3*i-2,[10 11 12]),Sl(3*i,[10 11 12]),Sl(3*i-1,[10 11 12]),'k+');
        plot3(Sl(3*i-2,[13 14 15]),Sl(3*i,[13 14 15]),Sl(3*i-1,[13 14 15]),'g+');

        for j=1:size(E,1)
            plot3( [Sl(3*i-2,E(j,1)) Sl(3*i-2,E(j,2))],...
                [Sl(3*i-0,E(j,1)) Sl(3*i-0,E(j,2))],...
                [Sl(3*i-1,E(j,1)) Sl(3*i-1,E(j,2))])
        end
        axis equal
        hold off

        pause(0.0001)
        end

end


X = zeros([3 F P]);
x = zeros([2 F P]);
for p=1:P, 
    X(:,:,p) = reshape(Sl(:,p),[3 F]); 
    x(:,:,p) =  [ X(1,:,p)./X(3,:,p);...
                  X(2,:,p)./X(3,:,p)];
end


return;
figure,plot(squeeze(x(1,1,:)),squeeze(x(2,1,:)),'bo')
for i=1200:1500
    plot(squeeze(x(1,i,:)),squeeze(x(2,i,:)),'bo')
    hold on;
    for j=1:size(E,1)
    plot( [x(1,i,E(j,1)) x(1,i,E(j,2))],[x(2,i,E(j,1)) x(2,i,E(j,2))])
    end
    axis equal
        hold off
    pause(0.00001)
end






P13 = [Sl(1:3:end,13) Sl(2:3:end,13) Sl(3:3:end,13)];
P14 = [Sl(1:3:end,14) Sl(2:3:end,14) Sl(3:3:end,14)];
   
%save('data5','X','x')







dx = P13-P14;    
dis = sqrt(sum(dx.^2,2));
figure,plot(dis)   
dx = dx/norm(dx(1,:));
figure,plot3(dx(1,1),dx(1,2),dx(1,3),'b.');
axis([-2 2 -2 2 -2 2]);
hold on;
[x,y,z] = sphere;
surf(x,y,z)

for i=2:5
    plot3(dx(i,1),dx(i,2),dx(i,3),'b.')
    pause(0.001);
end