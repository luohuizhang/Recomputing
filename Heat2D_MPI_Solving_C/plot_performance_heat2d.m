% Plot performances for MPI Heat2D %

% Loading input file
x=load('performances.txt');
% Get runtimes and number of processes
for i=1:3
runTime(1:9,i)=(x((i-1)*9+1:i*9,3));
end 
% Sequential time for each size
timeSeq(1:3) = runTime(1,1:3);
% Compute speedup values
for i=1:3
speedUp(1:9,i) = timeSeq(i)./runTime(1:9,i);
end
% Define histogram parameters
y= [1 2 4 8 16 32 64 128 256];
figure(1);
h=bar(log2(y),speedUp(1:9,1:3),'b'); % get histogram
set(gca, 'Xlim',[-1 9]);
set(gca, 'YLimMode', 'Auto');
set(gca,'xticklabel',{'1','2','4','8','16','32','64','128','256'});
set(h(1),'facecolor','b'); % use color name
set(h(2),'facecolor','r'); % or use RGB triple
set(h(3),'facecolor','g'); % or use a color defined in the help for plot
hPatch = findobj(h,'Type','patch');
set(hPatch,'facealpha',0.4); 
% Define legend and labels
h=legend('N=512^2','N=1024^2','N=2048^2');
rect = [0.66, 0.76, .15, .15];
set(h, 'Position', rect,'color','none');
hPatch = findobj(h,'Type','patch');
set(hPatch,'facealpha',0.4); 
xlabel('Number of Processes');
ylabel('SpeedUp');
title('Benchmark : SpeedUp vs Number of Processes');
