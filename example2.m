%% read and visualize from file with many snapshots
data = dlmread('example2.dump.out','',0,1);

nbead = data(1,1);
nconfig = size(data,1)/(data(1,1)+1)

%% show movie
b = nbead+1;
for cc = 1:nconfig
    beads = data(b*(cc-1)+2:b*cc,:);
    plot3(beads(:,1),beads(:,2),beads(:,3),'.-')
    view([1,0,0])
    axis equal
    drawnow
    pause(0.1)
end
