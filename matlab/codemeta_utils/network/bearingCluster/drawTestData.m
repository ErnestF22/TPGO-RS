function drawTestData(fileName)

pictureFileName=['testData/' fileName '.png'];
im=imread(pictureFileName);

mode='vertices';
endpoint=1;

x=zeros(2,0);
E=[];
figure(1)
while true
    imshow(im)
    hold on
    if ~isempty(x)
        plotPoints(x,'go','MarkerSize',15)
        if ~isempty(E)
            plot([x(1,E(:,1)); x(1,E(:,2))],[x(2,E(:,1)); x(2,E(:,2))],'g')
        end
    end
    hold off
    text(-50,size(im,1)/2,'Exit')
    text(size(im,2)/2,-10,mode)
    
    c=ginput(1);
    if c(1)<0
        break
    end
    if c(2)<0
        switch mode
            case 'vertices'
                mode='edges';
            case 'edges'
                mode='vertices';
        end
    else
        switch mode
            case 'vertices'
                x=[x c'];
            case 'edges'
                dx=sqrt(euclideanDistMatrix(x,c'));
                [mindx,idxMindx]=min(dx);
                if mindx<5
                    fprintf('Selected vertex %d\n',idxMindx)
                    switch endpoint
                        case 1
                            eprev=idxMindx;
                            endpoint=2;
                        case 2
                            ENew=[eprev idxMindx];
                            E=[E;ENew];
                            fprintf('Added edge [%d %d]\n',ENew)
                            endpoint=1;
                    end
                else
                    disp('No selection')
                end
        end
    end
end
close(1)

fileNameSave=['testData/' fileName];
save(fileNameSave)
disp(['Data saved to ' fileNameSave])


