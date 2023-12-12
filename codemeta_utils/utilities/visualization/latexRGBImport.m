function latexRGBImport

filename = 'latexRGB.txt';
formatSpec = '%3s%4s%4s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec,...
    'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
flagComment=cellfun(@(x) x(1)=='#', dataArray{1});

processColumn=@(idx) cell2mat(cellfun(@(x) str2double(x),dataArray{idx}(~flagComment),'UniformOutput',false));

rgbValues=[processColumn(1) processColumn(2) processColumn(3)]/255;
rgbNames=cellfun(@(x) strrep(x,' ',''),dataArray{4}(~flagComment),'UniformOutput',false);

[~,idxUnique]=unique(rgbNames);
rgbValues=rgbValues(idxUnique,:);
rgbNames=rgbNames(idxUnique);
structParams=vec([rgbNames num2cell(rgbValues,2)]');
rgb=struct(structParams{:});

save('latexRGB','rgb')
