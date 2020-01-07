function SI = search_information(A)

n = size(A,1);
mask_ut = triu(true(n,n),1);

A(A<eps) = eps;
A(mask_ut) = 0;

[~, SP, SPB] = get_shortest_path_lengths(1./A);
SI = get_information_shortest_paths_wei_und(A,SP,SPB,sum((A)),1);
