function MFPT = mean_first_passage_time(A)

n = size(A,1);
mask_ut = triu(true(n,n),1);

A(A<eps) = eps;
A(mask_ut) = 0;

P = f_markov_chain(A);
[~, MFPT, ~] = f_mfpt(P);