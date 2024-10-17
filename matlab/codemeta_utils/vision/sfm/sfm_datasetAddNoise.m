%Add random Gaussian noise to data.feature.locationNormalized
function data=sfm_datasetAddNoise(data,sigmaNoise)

NCameras=length(data.feature);
for iCamera=1:NCameras
    data.feature(iCamera).locationNormalized=....
        data.feature(iCamera).locationNormalized...
        +sigmaNoise*randn(size(data.feature(iCamera).locationNormalized));
end
