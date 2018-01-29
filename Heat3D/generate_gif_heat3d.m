function generate_gif_heat3d

% Load dimensions 
file=textread('param','%s','delimiter','\n');
x_dim=str2double(file(2));
y_dim=str2double(file(4));
z_dim=str2double(file(6));
x0=zeros(y_dim+2,x_dim+2,z_dim+2);
% Number of frames : must be equal to number of output files
numFrames=400;
% Time step between 2 frames
step=0.1;
animated(1,1,1,numFrames)=0;
char_f='%f';
for m=1:x_dim+1
  char_f=strcat(char_f,' %f');
end    
% Main loop on number of frames
for l=1:numFrames
  % Open current output file
  fid=fopen(strcat('outputPar',num2str(l),'.dat'),'r');        
  % Read all values
  for k=1:z_dim+2
    x=fscanf(fid,char_f,[x_dim+2 y_dim+2]);
    x=x';
    x0(1:y_dim+2,1:x_dim+2,k)=x;
    fgetl(fid);
  end
  % Close current output file
  fclose(fid);
  % Create meshgrid
  [x1 y1 z1]=meshgrid(1:x_dim+2,1:y_dim+2,1:z_dim+2);
  % Initialize plot
  hFig=figure(1);
  set(hFig,'Position',[400 400 750 600]);
  slice(x1,y1,z1,x0,[(x_dim+2)/2 x_dim+2], [(y_dim+2)/2 y_dim+2], [1 (z_dim+2)/2+1]);
  % Parameters for slice plot
  shading faceted;
  view([-42,22]);
  hc=colorbar;
  set(hc,'position',[0.932 0.3 0.02 0.6]);
  caxis([-10 10]);
  xlim([0 x_dim+2]);
  ylim([0 y_dim+2]);
  zlim([0 z_dim+2]);
  xlabel('x domain');
  label('y domain');
  zlabel('z domain');
  % Pause
  pause(step);
  % Get current frame from figure
  frame=getframe(figure(1));    
  if l==1
   [animated,cmap]=rgb2ind(frame.cdata,256,'nodither');
  else
   animated(:,:,1,l)=rgb2ind(frame.cdata,cmap,'nodither');
  end  
end
% Write final animated gif
imwrite(animated,cmap,'Heat_3D.gif','DelayTime',step,'LoopCount',inf); 
end
