function [results] = rotRef_contractionTest_ProcessResults(varargin)
% Load all *.mat files in "rotRef_contractionResults" directory or the
% provided file path
% INPUTS:
%   (optional) folderpath := a string with full file path to where all results of
%       rotRef_contractionTest.m is stored
% OUTPUTS:
%   results := a struct array with elements: kd, kv, kp, beta, mag_R,
%       mag_RRef, mag_W, 'filename.mat'

% Load optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'folderpath'
            % Initial gains must be sent as a 3x1 vector
            ivarargin=ivarargin+1;
            folderpath = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

if ~exist('folderpath','var')
    % Assume results are stored in the following subfolder
    folderpath = 'rotRef_contractionResults';
end
% Get all files in folder
resultFiles = dir(folderpath);
results = struct('kd',{},'kv',{},'kp',{},'beta',{},'mag_R',{},...
    'mag_RRef',{},'mag_W',{},'filename',{});

% Load data and
iCount = 1;
for ii = 1:size(resultFiles,1)
    if ~resultFiles(ii).isdir
        filename = [resultFiles(ii).folder '/' resultFiles(ii).name];
        % Run in try-catch since some files may be corrupted
        try
            load(filename,'kd','kv','kp','beta','mag_R','mag_RRef','mag_W');
            results(iCount).kd = kd;
            results(iCount).kv = kv;
            if exist('kp','var')
                results(iCount).kp = kp;
            end
            results(iCount).beta = beta;
            results(iCount).mag_R = mag_R;
            if exist('mag_RRef','var')
                results(iCount).mag_RRef = mag_RRef;
            end
            results(iCount).mag_W = mag_W;
            results(iCount).filename = filename;
            iCount = iCount + 1;
            clear kd kv kp beta mag_R mag_RRef mag_W filename
        catch
        end
    end
end
end

