

%TODO: calc_A does not work, return -10000 delta=1??
% why can delta=1???

%max_its_left can be negative

%find 0.7 inlier ratio, but still shit model?? negative its

[pts, pts_tilde, A_true, t_true] = affine_test_case(1000, 0, 1);
[A_best,t_best, inlier_set, iters] = rransac_sprt(pts, pts_tilde, 9);
[A, t] = least_squares_affine(pts(:,inlier_set), pts_tilde(:, inlier_set));
iters
A_true-A_best
t_true-t_best
A_true-A
t_true-t


