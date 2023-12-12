sclear

desc = [];
pts = [];
imageLabel = [];

% load ../dataset/FG3DCar/selected-sedan
% datapath = '../dataset/multilane/location1';
datapath = '../../../data/graffiti/wall';
filelist = dir([datapath '/*.png']);
imgList = {filelist.name};
nImg = length(imgList);

for i = 1:nImg
    load(sprintf('%s/%s.view.mat',datapath,imgList{i}));
    im{i} = viewInfo.img;
    desc = [desc,viewInfo.desc];
    pts = [pts,single(viewInfo.frame(1:2,:))];
    imageLabel = [imageLabel,single(i*ones(1,size(viewInfo.desc,2)))];
end

D = pdist(desc');
D = squareform(D);

%%
[membership,info]=quickshift_matching(D,imageLabel,'threshold',0.5,'ratioInterCluster',0.1);

clear cluster
for i = 1:max(membership)
    ptsID = find(membership==i);
    cluster(i).ptsID = ptsID;
    cluster(i).density = mean(info.density(ptsID));
    cluster(i).counts = length(ptsID);
end 

figure('position',[100,100,800,500]);
ha = tight_subplot(2,3,0.01,0.01,0.01);
[~,ind] = sort([cluster.density],'descend');
for j = 1:nImg
    axes(ha(j));
    imshow(im{j});
    hold on
end

N = 10;
cmap = hsv(N);
i = 0;
for clusterID = ind
    cluster(clusterID).density
    cluster(clusterID).counts
    if cluster(clusterID).counts == nImg
        i = i+1;
        for j = cluster(clusterID).ptsID
            axes(ha(imageLabel(j)));
            plot(pts(1,j),pts(2,j),'+','MarkerEdgeColor',cmap(i,:),'markersize',14,'linewidth',3);
        end
%         pause
    end
    if i >= N
        break
    end
end
