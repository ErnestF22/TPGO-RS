function icra13_plots
close all
flagSaveFigure=2;
figureDir='../../../papers/control/ICRA13-BearingControl/figures';
figSize=180*[1,0.8]*0.85;
fs=setFigFontSize(7);
%setFigFont('CMU Sans Serif')

optsCostDisplay={'optsX',{'Marker','o','MarkerFaceColor','b','MarkerSize',5}};
optsStartLocation={'ro','MarkerFaceColor','r','MarkerSize',5};
for iPlot=[1]
    switch iPlot
        case 1
            load('icra13_singleTrajectoryComparison_data')
            baseFileName=fullfile(figureDir,'trajectories');
            for iFun=1:NFuns
                fName=funsName{iFun};
                figure
                funs=bearingCostFunctions(fName);
                bearingCostDisplay(sceneData,funs,optsCostDisplay{:})

                for ix0=1:LPoints
                    for jx0=1:LPoints
                        x0=result(iFun).x0{ix0,jx0};
                        x=result(iFun).x{ix0,jx0};
                        hold on
                        plot(x0(1),x0(2),optsStartLocation{:})
                        plot(x(1,:),x(2,:),'k')
                        hold off
                    end
                end
                axis([-10 10 -8 10])
                set(gca,'XTick',[])
                set(gca,'YTick',[])
                savefigure([baseFileName 'Integrator' fName '_' num2str(LPoints)],'epsc',figSize,flagSaveFigure);
            end
        case 2
            load('icra13_singleTrajectoryComparisonUnicycle_data')
            baseFileName=fullfile(figureDir,'trajectories');
            style={{'b'},{'r'},{'g'},{'m'}};

            NFuns=length(funsName);
            for iFun=1:NFuns
                fName=funsName{iFun};
                figure
                funs=bearingCostFunctions(fName);
                bearingCostDisplay(sceneData,funs,'flagContour',false,optsCostDisplay{:})
                pause(0.01)

                for ix0=1:LPoints
                    for jx0=1:LPoints
                        for kx0=1:4
                            x=result(iFun).x{ix0,jx0,kx0};
                            x0=result(iFun).x0{ix0,jx0,kx0};
                            hold on
                            plot(x0(1),x0(2),optsStartLocation{:})
                            plot(x(1,:),x(2,:),style{kx0}{:})
                            hold off
                        end
                    end
                end
                axis([-10 10 -8 10])
                set(gca,'XTick',[])
                set(gca,'YTick',[])
                savefigure([baseFileName 'Unicycle' fName '_'  num2str(LPoints)],'epsc',figSize,flagSaveFigure);
            end
        case 3
            %figSize=180*[1,1]*0.79;
            for nLandmarkSet=1:4
                load(['icra13_varianceBearingLocalization_' num2str(nLandmarkSet) '_data'])
                vContour=linspace(0,13,200);
                baseFileName=fullfile(figureDir,'varianceMap');
                for iD=1:ND
                    figure(iD)
                    z1=log(z{iD});
                    %contourf(xx,yy,z1,vContour,'LineColor','none');
                    imagesc(x,x,z1,[0 13])
                    colormap(flipud(gray))
                    disp(min(z1(:)))
                    disp(max(z1(:)))
                    hold on
                    plot(XLandmarks{iD}(1,:),XLandmarks{iD}(2,:),'bo','MarkerSize',5,'MarkerFaceColor','b');
                    hold off
                    axis equal
                    axis tight
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    savefigure([baseFileName '_S' num2str(nLandmarkSet) '_D' num2str(iD)]...
                        ,'epsc',figSize,flagSaveFigure);
                end
            end            
    end
end
setFigFontSize(fs);