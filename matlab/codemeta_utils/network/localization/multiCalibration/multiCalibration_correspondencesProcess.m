%function c=multiCalibration_correspondencesProcess(c)
%Process raw data in the correspondences loaded by
%multiCalbration_datasetToCorrespondencesFile() and add fields with the
%data which can be used by multiCalibration
%In detail, it uses the field 'type' of each element in the cell array:
%   'plane-plane'   normalize the raw data and produce fields n1,d1,n2,d2
%                   with the two plane information
%   'plane-point'   normalize the raw data of the plane and produce fields
%                   n1,d1, then reshape raw data of points to produce field
%                   x22D
function c=multiCalibration_correspondencesProcess(c)
NCorrespondences=length(c);
for iCorrespondence=1:NCorrespondences
    cCurrent=c{iCorrespondence};
    switch cCurrent.type
        case 'plane-plane'
            [n1,d1]=cnormalize(cCurrent.rawData{1}');
            d1=-d1;
            [n2,d2]=cnormalize(cCurrent.rawData{2}');
            d2=-d2;
            cCurrent.n1=n1;
            cCurrent.d1=d1;
            cCurrent.n2=n2;
            cCurrent.d2=d2;
        case 'plane-point'
            [n1,d1]=cnormalize(cCurrent.rawData{1}');
            d1=-d1;
            n1=reshape(repmat(n1,[2 1]),3,[]);
            d1=reshape(repmat(d1,[2 1]),1,[]);
            %use num. of dim. of matrix to determine version of the data
            %and use the right extraction code
            sz=size(cCurrent.rawData{2});
            if length(sz)==2
                x22D=reshape(cCurrent.rawData{2}(:,[2 3 5 6])',2,[]);
            else
                x22D=reshape(permute(cCurrent.rawData{2}(:,:,2:3),[3 2 1]),2,[]);
            end
            cCurrent.n1=n1;
            cCurrent.d1=d1;
            cCurrent.x22D=x22D;
    end
    c{iCorrespondence}=cCurrent;
end
