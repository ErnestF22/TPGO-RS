function ReluPlexPlotNetwork(layer,L)
nbColors=32;
colors=parula(nbColors);
c=0;
gs1=getSeqCombSet([0 1],2);
x=linspace(-15,15,200);


for is1=1:4
    s1=gs1()';
    gs2=getSeqCombSet([0 1],2);
    for is2=1:4
        s2=gs2()';
        gs3=getSeqCombSet([0 1],1);
        for is3=1:2
            s3=gs3()';
%             gs4=getSeqCombSet([0 1],1);
%             for is4=1:2
%                 s4=gs4()';
                %             disp([s1; s2; s3]);
                S = [s1;s2;s3];
                figure(1)
                Y = funImage(x,x,@(x) linearRestrict(x,layer,L,S)...
                    ,'method','surf','methodOpts',{'FaceColor',colors(mod(c-1,nbColors)+1,:)});
                hold on
                c=c+1;
                if ~all(all(isnan(Y)))
                    legendInfo{c} = ['Z = ' num2str([s1',s2',s3])];
                else
                    legendInfo{c} = [''];
                end
                hold on
%             end
        end
    end
end
% legendInfo = legendInfo(~cellfun('isempty',legendInfo));
legend(legendInfo)
xlabel('x')
ylabel('y')

end