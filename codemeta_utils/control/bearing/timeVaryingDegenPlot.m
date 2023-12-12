switch 5
    case 1 %degen 1
        xs = [1.5 1; 1 2];
        plot(xs(1,:),xs(2,:),'kx');
        hold on
        plot([0 2],[4 0],'k:');
        x=linspace(0,2.5);
        y=-0.2*cos(2*x)+1.5;
        %ydot=0.4*sin(2*x)
        plot(x,y,'r');
        point=[1.1791, 1.1791-.1, 1.1791+.1];
        pointy=-2.*point+4;
        plot(point,pointy,'b*');
        ydot=0.4*sin(2*point);
        %quiver(point,pointy,1,ydot,'b');
        axis([0.5 1.75 0.75 2.25])
        point2=0.75;
        point2y=-0.2*cos(2*point2)+1.5;
        plot([xs(1,1),point2],[xs(2,1),point2y],'k:',[xs(1,2),point2],[xs(2,2),point2y],'k:');
        plot(point2,point2y,'c*')
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        hold off
        clear
    case 2 %degen 2
        xs = [1 1.2;1 1.7];
        plot(xs(1,:),xs(2,:),'kx');
        hold on
        x=linspace(0,2.5);
        y=-0.2*cos(4*x)+1.5;
        plot(x,y,'r');
        point = 0.7854;
        pointy = -0.2*cos(4*point)+1.5;%1.7 , slope -3.2618, 0.30657
        plot([xs(1,1),point],[xs(2,1),pointy],'k:',[xs(1,2),point],[xs(2,2),pointy],'k:');
        plot(point,pointy,'b*');
        q=quiver(point, pointy, 0.2, 0, 0,'b');
        q.MaxHeadSize = 0.5;
        q2=quiver(xs(1,1),xs(2,1),-0.14,-0.14*0.30657,0,'k');
        q2.MaxHeadSize = 0.5;
        hold off
        axis([0.4 1.6 0.8 2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        clear
    case 3 %prop 1
        xs = [0.5 1; 1.5, -2];
        plot(xs(1,:),xs(2,:),'kx');
        hold on
        point=[-0.5,0];
        plot([-4,1],[-3,4.5],':k',[-5,1],[4,-4],'k:');%slope 1.5, -1.333; -2/3, 0.75
        q=quiver(point(1),point(2),1.5,0,0,'r');
        q.MaxHeadSize = 0.5;
        q2=quiver(point(1),point(2),-1.5,0,0,'--r');
        q2.MaxHeadSize = 0.5;
        q3=quiver(xs(1,1),xs(2,1),-0.5385-.5,2.192-1.5,0,'k');
        q3.MaxHeadSize = 0.5;
        q4=quiver(xs(1,2),xs(2,2),0.04-1,-2.72+2,0,'k');
        q4.MaxHeadSize = 0.5;
        plot(point(1),point(2),'r*');
        hold off
        axis([-3.5 3.5 -3.5 3.5]);
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        clear
    case 4 %prop 2
        xs = [1.5 -.5; 2 -2];
        plot(xs(1,:),xs(2,:),'kx');
        hold on
        point = [-1;1.5];
        plot(0,0,'b*');
        plot([0,xs(1,1)],[0,xs(2,1)],'b',[0,xs(1,2)],[0,xs(2,2)],'b');
        plot([-3.5,4],[1,2.5],'k:',[-1.5,0],[5,-5.5],'k:');
        plot(point(1),point(2),'r*');
        hold off
        axis([-3 3 -3 3]);
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        clear
    case 5 %2d planar solution
        xs = [1 1 1.1875; 2.25, -2, -2.25];
        plot(xs(1,1:2),xs(2,1:2),'kx');
        hold on
        %plot(xs(1,3),xs(2,3),'rx');
        point=[-0.5,0];
        plot([-4,1],[-3,4.5],':k',[-5,1],[4,-4],'k:');%slope 1.5, -1.333; -2/3, 0.75
        %plot([-2,4],[4.2502,-3.7498],'r:');
        %q=quiver(point(1),point(2),0,2.25,0,'--r');
        %q=quiver(point(1),point(2), 0, 2.25, 0,'--r');
        %q.MaxHeadSize = 0.5;
        q2=quiver(point(1),point(2),-1.5,0,0,'--b');
        q2.MaxHeadSize = 0.5;
        q3=quiver(xs(1,1),xs(2,1),-0.5385-.5,2.192-1.5,0,'k');
        q3.MaxHeadSize = 0.5;
        q4=quiver(xs(1,2),xs(2,2),0.04-1,-2.72+2,0,'k');
        q4.MaxHeadSize = 0.5;
        %q5=quiver(xs(1,3),xs(2,3),2.2676-1.1875,-1.44+2.25,0,'r');
        %q5.MaxHeadSize = 0.5;
        plot(point(1),point(2),'b*');
        hold off
        axis([-3.5 3.5 -3.5 3.5]);
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        clear
    case 6 %2d parallel incorrect
        xs = [1 1 1.1875; 2.25, -2, -2.25];
        plot(xs(1,1:2),xs(2,1:2),'kx');
        hold on
        plot(xs(1,3),xs(2,3),'rx');
        point=[-0.5,0];
        plot([-4,1],[-3,4.5],':k');%slope 1.5, -1.333; -2/3, 0.75
        plot([-2,4],[4.2502,-3.7498],'r:');
        q=quiver(point(1),point(2), 0, 2.25, 0,'--r');
        q.MaxHeadSize = 0.5;
        q3=quiver(xs(1,1),xs(2,1),-0.5385-.5,2.192-1.5,0,'k');
        q3.MaxHeadSize = 0.5;
        q5=quiver(xs(1,3),xs(2,3),2.2676-1.1875,-1.44+2.25,0,'r');
        q5.MaxHeadSize = 0.5;
        plot(point(1),point(2),'k*');
        hold off
        axis([-3.5 3.5 -3.5 3.5]);
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        clear
end