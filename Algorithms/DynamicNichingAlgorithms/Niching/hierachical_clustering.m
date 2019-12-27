function niching = hierachical_clustering(decs, K)
    matdis = pdist2(decs, decs);
    
    tree_cluster = linkage(matdis);
    
    idx = cluster(tree_cluster, 'maxclust', K);
    
    niching = cell(max(idx), 1);
    for i = 1:length(niching)
        niching{i} = find(idx == i);
    end
end