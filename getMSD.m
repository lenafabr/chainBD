function [MSD,ste,cnt,MSDindv,steindv,cntindv] = msdFromTracks(tracklist)
% calculate MSD for individual tracks and for a collection of tracks as a
% whole
% WARNING: does not properly account for nonconsecutive frames

allvals = [];


for tc = 1:length(tracklist)
    track = tracklist{tc}(:,1:3);
    
    lt = size(track,1);
    
    MSDindv{tc} = [];
    cntindv{tc} = [];
    steindv{tc} = [];
    
    for dc = 1:lt-1           
        if (dc>length(allvals)); allvals{dc}=[]; end
        diffs = track(dc+1:dc:end,:)-track(1:dc:end-dc,:);
        
        sqdiffs = diffs(:,1).^2 + diffs(:,2).^2+diffs(:,3).^2;

        % statistics for individual tracks
        MSDindv{tc}(dc) = mean(sqdiffs);
        cntindv{tc}(dc) = length(sqdiffs);
        steindv{tc}(dc) = std(sqdiffs)/sqrt(cntindv{tc}(dc));
        
        % put all msd for this timeinterval together
        allvals{dc} = [allvals{dc}; sqdiffs]; 
    end
end

length(allvals)
for dc = 1:length(allvals)
    MSD(dc) = mean(allvals{dc});
    stdval = std(allvals{dc}')';
    cnt(dc) = length(allvals{dc});
    ste(dc) = std(allvals{dc})/cnt(dc);
end