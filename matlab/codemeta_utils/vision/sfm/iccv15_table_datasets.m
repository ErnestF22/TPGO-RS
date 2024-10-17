function iccv15_table_datasets
datasets={'fountain','herzjesu','herzjesularge','castle','castlelarge','castleentry'};
NDatsets=length(datasets);

fid=fopen('../../../papers/vision/ICCV15-RotationAveragingSurvey/tableDatasets.tex','w');
fprintfTee(fid,'Name \t& \\# poses\t& \\# edges\t& \\%% edges\\\\\n')
fprintfTee(fid,'\\midrule\n')
for iDataset=1:NDatsets
    thisDataset=datasets{iDataset};
    s=load(['sfmdata_' thisDataset]);
    NPoses=length(s.data.imgFileName);
    NEdges=length(s.data.matchFiltered);
    percEdges=NEdges/length(s.data.match)*100;
    fprintfTee(fid,'%s \t& %d \t& %d \t& %.2f \\\\\n',...
        thisDataset, NPoses, NEdges, percEdges);
end
fclose(fid);

function fprintfTee(fid,varargin)
fprintf(varargin{:})
fprintf(fid,varargin{:});
