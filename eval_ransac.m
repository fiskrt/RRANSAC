
iter_all = {}; times_all = {};
for d=0:3
    iterations = [];
    times = [];
    for outlier_ratio=0.9:0.1:0.9
        [pts, pts_tilde, A_true, t_true] = affine_test_case(2000, outlier_ratio, 1);
        tot_it = 0;
        tic;
        N = 10;
        for j=1:N
            [A_best,t_best, inlier_set, iters] = rransac_tdd(pts, pts_tilde, 9, d);
            tot_it = tot_it + iters;
        end
        fin = toc;
        iterations(end+1) = tot_it/N;
        times(end+1) = fin/N;
    end
    iter_all{end+1} = iterations
    times_all{end+1} = times
end
%% Plot values
figure;
for i=1:4
    semilogy(cell2mat(iter_all1(i)), 'DisplayName',"d="+ num2str(i-1))
    hold on;
end
xticklabels(linspace(0, 0.9, 10))
xlabel('Outlier ratio')
ylabel('Number of iterations')
legend('Location','NorthWest' )
title('RANSAC number of iterations')
figure;
for i=1:4
    semilogy(cell2mat(times_all1(i)), 'DisplayName',"d="+ num2str(i-1))
    hold on;
end
xticklabels(linspace(0, 0.9, 10))
xlabel('Outlier ratio')
ylabel('Running time [s]')
legend('Location','NorthWest' )
title('RANSAC running time')

%% From now on we use d=1.
% With and without local opt for inlier ratios 0:0.6
% Record: run-time, the number of iterations, the inlier ratio,
% number of inliers, and the average residual length among all true inliers
n_points = 1000;
%sd = 1;
outlier_ratio = 0.4;
sds = [0.5,1,2,4,8,16];
iterations = []; times = []; n_inliers = []; n_inlier_ratios = []; res_lengths = [];
%for outlier_ratio=0:0.1:0.6
for sd=sds
    [pts, pts_tilde, A_true, t_true] = affine_test_case(n_points, outlier_ratio, sd);
    n_true_inlier = n_points-ceil(n_points*outlier_ratio);
    tot_it = 0; tot_inl = 0; tot_res_len = 0; N = 10;
    tic;
    for j=1:N
        [A_best,t_best, inlier_set, iters] = rransac_tdd(pts, pts_tilde, (3*sd).^2, 1);
        %[A_best,t_best, inlier_set, iters] = rransac_sprt(pts, pts_tilde, (3*sd).^2);
        %[A_best, t_best] = least_squares_affine(pts(:,inlier_set), pts_tilde(:, inlier_set));
        tot_it = tot_it + iters;
        tot_inl = tot_inl + sum(inlier_set);
        tot_res_len = tot_res_len + mean(residual_lengths( ...
            A_best, t_best,pts(:, 1:n_true_inlier), pts_tilde(:, 1:n_true_inlier)));
    end
    fin = toc;
    iterations(end+1) = tot_it/N;
    times(end+1) = fin/N;
    n_inliers(end+1) = tot_inl/N;
    n_inlier_ratios(end+1) = tot_inl/(N*n_points);
    res_lengths(end+1) = tot_res_len/N;
end

%Tab = [transpose(0:0.1:0.6) iterations' times' n_inliers' n_inlier_ratios' res_lengths']
Tab = [sds' iterations' times' n_inliers' n_inlier_ratios' res_lengths']

%% For outlier ratio 0.4, we try to vary the noise, i.e sigma.
% threshold = (3*sigma)^2

[pts, pts_tilde, A_true, t_true] = affine_test_case(2000, 0.01, 1);

%tic;
summ = 0;
for i=1:1
    [A_best,t_best, inlier_set, iters] = rransac_tdd(pts, pts_tilde, 9, 1);
    summ = summ + iters;
end
%toc
summ

[pts, pts_tilde, A_true, t_true] = affine_test_case(2000, 0.95, 1);
%tic;
summ = 0;
for i=1:1
    [A_best,t_best, inlier_set, iters] = rransac_tdd(pts, pts_tilde, 9, 1);
    summ = summ + iters;
end
%toc
summ



