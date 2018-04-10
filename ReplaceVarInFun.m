function   ReplaceVarInFun(filename,replaceLine,newText)

%this function open a file as text and replaces a specific line with
%another one
%intended to be used for post processing in witlis when something needs to
%be changed iteratively in the processingconfig file

%get number of lines
fid = fopen(filename);
res={};
while ~feof(fid)
  res{end+1,1} =fgetl(fid); %#ok<AGROW>
end
fclose(fid);
numLines=numel(res);


fileID = fopen(filename,'r');
mydata = cell(1, numLines);
for k = 1:numLines
   mydata{k} = fgetl(fileID);
end
fclose(fileID);

mydata{replaceLine} = newText;

fileID = fopen(filename,'w');
fprintf(fileID,'%s\n',mydata{:});
fclose(fileID);


end