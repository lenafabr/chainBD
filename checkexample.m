% look at example snapshot
data = dlmread('example.dump.out','',0,1);
beads = data(2:end,:);
plot3(beads(:,1),beads(:,2),beads(:,3),'.-')

%% read and visualize from file with many snapshots
data = dlmread('example.dump.out','',0,1);

nbead = data(1,1);
nconfig = size(data,1)/(data(1,1)+1)

beads = data(end-nbead+1:end,:);
plot3(beads(:,1),beads(:,2),beads(:,3),'.-')