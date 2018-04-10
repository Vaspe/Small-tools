function MakeInvisibleFigVisible(Directory)

% When batch of plots are saved in an invisible mode their 'visible' property is 
%set to off and cannot be seen when opened. This functions patches the
%issue by opening opening and altering the fig properties and re-saving
%them. Only the absolute path of the folder containing the figures is
%needed. 

curdir=pwd;
cd (Directory)
fileN_int=dir('*.fig');
for kk=1:length(fileN_int)
    ADfiles{kk,1}=[fileN_int(kk).name];
    aa=open(ADfiles{kk,1});
    aa.Visible='on';
    saveas(gcf, ADfiles{kk,1})
    close(aa)
    clear aa
end
cd (curdir)