function niching = nearest_better_clustering(decs, fit, phi)
	[~, sort_idx] = sort(fit, 'descend');
    matdis = pdist2(decs(sort_idx, :), decs(sort_idx, :));
    
    n=length(matdis);
    nbc=zeros(n,3);
    nbc(1:n,1)=1:n;
    nbc(1,2)=-1;
    nbc(1,3)=0;
    for i=2:n
        [u,v]=min(matdis(i,1:i-1));
        nbc(i,2)=v;
        nbc(i,3)=u;
    end
    meandis=phi*mean(nbc(2:n,3));
    nbc(nbc(:,3)>meandis,2)=-1;
    nbc(nbc(:,3)>meandis,3)=0;
    seeds=nbc(nbc(:,2)==-1,1);
    m=zeros(n,2);
    m(1:n,1)=1:n;
    for i=1:n
        j=nbc(i,2);
        k=j;
        while j~=-1
            k=j;
            j=nbc(j,2);
        end
        if k==-1
            m(i,2)=i;
        else
            m(i,2)=k;
        end
    end

    % construct the result
    niching = cell(length(seeds), 1);
    for i=1:length(seeds)
        niching{i} = sort_idx(m(m(:, 2) == seeds(i), 1));
    end 
end

