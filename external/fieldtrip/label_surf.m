function newsurf = label_surf(tets, tolabel, neighborlabel)


% LABEL_SURF label boundaries of given regions neighboring (an)other
% region(s)

% boundary labeling


facenb=faceneighbors(tets(:,1:4));
newsurf = [];

Lia = find(ismember(tets(:,5),tolabel));
for ii = 1:length(Lia)
    for jj = 1:4
        if facenb(Lia(ii),jj) == 0
            continue
        end
        if ismember(tets(facenb(Lia(ii),jj),5),neighborlabel)
            switch jj
                case 1
                    newsurf = [newsurf; tets(Lia(ii),[1,2,3,5])];
                case 2
                    newsurf = [newsurf; tets(Lia(ii),[1,2,4,5])];
                case 3
                    newsurf = [newsurf; tets(Lia(ii),[1,3,4,5])];
                case 4
                    newsurf = [newsurf; tets(Lia(ii),[2,3,4,5])];
            end
        end
    end
end