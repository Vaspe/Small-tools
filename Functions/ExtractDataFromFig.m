%GET DATA FROM FIGURE
function [x,y] = ExtractDataFromFig(figname,OneXaxis)

openfig(figname,'invisible');
fh=gcf;
h = gca; 
D = get(h,'children');

y  = get(D,'Ydata');
x1 = get(D,'Xdata');
if OneXaxis==1
    x = x1(1);
else
    x = x1;
end
% axesObjs = get(h, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
% 
% ydata = get(dataObjs, 'YData');
close (fh)