switch 1
    case 1 %top
        %plot(simout.data(:,1),simout.data(:,2),'b','LineWidth',1)
        plot(data(:,1),data(:,2),'k','LineWidth',1.2)
        hold on
        plot(data(:,3),data(:,4),data(:,5),data(:,6),data(:,7),data(:,8),data(:,9),...
            data(:,10),data(:,11),data(:,12),data(:,13),data(:,14))
        %plot(simout.data(:,3),simout.data(:,4),'r')
        plot(simout1.data(1,:,1),simout1.data(2,:,1),'kx')
        axis([-10 8 -8 10])
        %axis([-1.5 1.5 -1.5 1.5])
        hold off
    case 2 %position vs time
        plot(simout.time(:),simout.data(:,1))
        hold on
        grid on
        plot(simout.time(:),simout.data(:,2))
        plot(simout.time(:),simout.data(:,3))
        plot(simout.time(:),simout.data(:,4))
        %axis([0 26 -2 1.5])
        axis([0 50 -9 5])
        hold off
end