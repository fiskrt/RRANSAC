function [A_best,t_best, inlier_set, i] = rransac_tdd(pts, pts_tilde, threshold, d)
    % RRANSAC with the Tdd test using local refinement
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    max_inliers = 0; A_best = []; t_best = []; inlier_set = [];
    n=3;
    max_iters = 2500000;
    k=max_iters;
    for i = 1:max_iters
        ind = randsample(size(pts,2), n);
        [A, t] = estimate_affine(pts(:,ind), pts_tilde(:,ind));

        if ~Tdd_test(d, A, t, pts, pts_tilde, threshold)
            continue
        end

         % We passed the T_dd test! Count the inliers
        inliers = residual_lengths(A, t, pts, pts_tilde) <= threshold;
        n_inliers = sum(inliers);

        if (n_inliers > max_inliers)
            max_inliers = n_inliers;
            inlier_set = inliers;
            A_best = A;
            t_best = t;

            % Keep an estimate of the inlier ratio to determine the number
            % of iterations
            n_inliers_lo = n_inliers;
            %[n_inliers_lo, inlier_set] = local_opt(pts, pts_tilde, n_inliers, inliers, threshold);
            inlier_ratio = n_inliers_lo/(length(pts));
            k = log(0.05)/log(1-inlier_ratio.^(n+d));
        end
        % Run max k iterations
        if i>k
            break
        end
    end
end


% Performs the T(d,d) test as defined in [chum2002].
% I.e a subset of size d<<N (N total number of points) is
% tested for to see if the satisfy the threshold.
function pass = Tdd_test(d, A, t, pts, pts_tilde, threshold)
    % The T_{d,d} test. Evaluates the model (A,t) on d points.
    % If all d points are not consitent the test is not passed.

    if (d==0)
        pass = true;
        return
    end

    ind = randsample(size(pts,2), d);
    n_inliers = sum(residual_lengths(A, t, pts(:,ind), pts_tilde(:,ind)) <= threshold);
    if (n_inliers == d)
        pass = true;
    else
        pass = false;
    end
end