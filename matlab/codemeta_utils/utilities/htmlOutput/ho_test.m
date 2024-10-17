%function ho_test
%Example of use of the group of functions ho_* for writing HTML pages

function ho_test

fileName='test.html';

%open file for writing
fid=fopen(fileName,'w');

%insert HTML header
ho_start(fid,'Title of the test page')

%insert one level 1 and one level 2 titles
ho_title(fid,'Big title: test page');
ho_title(fid,'Subtitle: this is really a test page',2);

%insert a paragraph of text
ho_text(fid,'This is an interesting paragraph of text.')

%start a table
ho_tablestart(fid,char('Titles','Data1','Data2'));
%insert rows of numbers, strings, formatted numbers
ho_tablerow(fid,'Numbers',[1 2])
ho_tablerow(fid,'Percentages',[1 2],'format','%.3f %%')
ho_tablerow(fid,'Strings',char('ciao','mondo'))
ho_tablerow(fid,'Numbers formatted<br />as strings',char(num2str(1,'%0.3f %%'),num2str([1 2 3])))
%finish the table
ho_tableend(fid)

%draw, save, and insert a figure
hf=figure();
plot(sin(linspace(0,3*pi)));
ho_figure(fid,'testImage','currentFigure')
close(hf);

%insert the content of a matrix
A=magic(8);
ho_matrix(fid,A,'format','%.3f');

%insert HTML footer
ho_end(fid);

fclose(fid);

disp(['HTML written in file ' fileName])