function generate_gif_heat2d

% Number of frames
numFrames=400;
% Time step between 2 frames
step = 0.1;
animated(1,1,1,numFrames) = 0;

% Main loop
for i=1:numFrames
figure(1);
% Set sizex and sizey from 'param' file
fileParam = fopen('param','r');
data=textscan(fileParam,'%s');
sizex=str2num(data{1}{2});
sizey=str2num(data{1}{4});
fclose(fileParam);
% Load data and plot with surfc
[X,Y]=meshgrid(0:sizex+1,0:sizey+1);
Z=load(strcat('outputPar',num2str(i),'.dat'));
surfc(X,Y,Z);
shading interp;
view([0,0,1]);
hc=colorbar;
set(hc,'position',[0.932 0.3 0.02 0.6]);
caxis([-10 10]);
xlabel('x domain');
ylabel('y domain');
zlabel('temperature');
xlim([0 sizex+1]);
ylim([0 sizey+1]);
% Set pause 
pause(step);    
frame = getframe(figure(1));    
if i == 1
  [animated, cmap] = rgb2ind(frame.cdata, 256, 'nodither');
else
  animated(:,:,1,i) = rgb2ind(frame.cdata, cmap, 'nodither');
end  
end
% Write final animated gif
imwrite(animated,cmap,'Heat_2D.gif','DelayTime',step,'LoopCount',inf); 
end
