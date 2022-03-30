function [pts, pts_tilde, A_true, t_true] = affine_test_case(n_points, outlier_rate, noise_sd)
    n_outlier = ceil(n_points*outlier_rate);

    A_true = [2 5; 1 1];
    t_true = [10;-25];

    pts = rand(2,n_points-n_outlier);
    pts(1,:) = pts(1,:).*640;
    pts(2,:) = pts(2,:).*480;
    pts_tilde = A_true*pts + t_true;
    % Add Gaussian noise to points with mean 0 and std 'noise_sd'
    pts_tilde = pts_tilde + normrnd(0,noise_sd,2,n_points-n_outlier);

    % Set the magnitude of the outliers to be 10 times as big as
    % the inliers
    outlier_range = 100*max(pts(:));
    pts = [pts outlier_range.*rand(2,n_outlier)];
    pts_tilde = [pts_tilde outlier_range.*rand(2,n_outlier)];
end



