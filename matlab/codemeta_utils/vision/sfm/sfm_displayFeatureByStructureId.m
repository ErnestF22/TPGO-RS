function sfm_displayFeatureByStructureId(data,iStructure,varargin)
structureMemberName='structure';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'member'
            ivarargin=ivarargin+1;
            structureMemberName=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

[x,iImageList]=sfm_getFeatureLocationByStructureId(data,iStructure,'member',structureMemberName);
NImages=length(iImageList);
NCols=ceil(sqrt(NImages));
NRows=ceil(NImages/NCols);
for iiImage=1:NImages
    iImage=iImageList(iiImage);
    img=sfm_getImageById(data,iImage);
    subplot(NCols,NRows,iiImage)
    sfm_rawDisplayFeature(img,x(:,iiImage));
end
