function [n_inlier_max, inliers_max] = local_opt(pts, pts_tilde, n_inliers, inliers, threshold)
    % Gets the maximum inlier count from a model refined using 
    % least squares on a samples of size min(n_inliers, 12).
    % The maximum is out of 10 tries.

    n_inlier_max = n_inliers;
    inliers_max = inliers;
    for k=1:10
        inlier_pts = pts(:,inliers);
        inlier_pts_tilde=pts_tilde(:,inliers);
        ind = randsample(size(inlier_pts,2), min(n_inliers, 12));
        [A_ls, t_ls] = least_squares_affine(inlier_pts(:,ind), inlier_pts_tilde(:,ind));
        inliers_lo = residual_lengths(A_ls, t_ls, pts, pts_tilde) <= threshold;
        n_inliers_lo = sum(inliers_lo);
        if n_inliers_lo > n_inlier_max
            n_inlier_max = n_inliers_lo;
            inliers_max = inliers_lo;
        end
    end
end