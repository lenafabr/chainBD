
%% read and visualize from file with many snapshots
data = dlmread('example1.dump.out','',0,1);

nbead = data(1,1);
nconfig = size(data,1)/(data(1,1)+1)

% look at last snapshot
beads = data(end-nbead+1:end,:);
plot3(beads(:,1),beads(:,2),beads(:,3),'.-')

%% show movie
b = nbead+1;
for cc = 1:nconfig
    beads = data(b*(cc-1)+2:b*cc,:);
    plot3(beads(:,1),beads(:,2),beads(:,3),'.-')
    %view([1,0,0])
    xlim([-10,10])
    ylim([-10,10])
    zlim([-10,10])
    drawnow
    pause(0.1)
end

%% look at MSD of end bead over time
data = dlmread('example1.out');

track = data(:,3:5);
tracklist = {track};

MSD = getMSD(tracklist);

%%
times = (1:length(MSD))*1e-3*100;
loglog(times,MSD,times,6*times,times,6*times.^0.5)
xlabel('time')
ylabel('MSD')