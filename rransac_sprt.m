function [A_best,t_best, inlier_set, i] = rransac_sprt(pts, pts_tilde, threshold)
    % RRANSAC with SPRT using local refinement
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % SPRT variables
    inlier_ratio = 0.2; delta_est = 0.01;
    epsilon_i = [inlier_ratio]; delta_i = [delta_est];
    A_i = [calculate_A(inlier_ratio, delta_est)];    
    k_i = [0]; % Number samples for test i
    
    % Other variables
    n=3; % minimal sample size
    test_num = 1;
    updated_test = false;
    n_rejected = 1;
    max_inliers = 0;
    max_iters = 100000;
    for i = 1:max_iters
        ind = randsample(size(pts,2), n);
        [A, t] = estimate_affine(pts(:,ind), pts_tilde(:,ind));
        k_i(test_num) = k_i(test_num) + 1;
       
        [pass, n_tested, n_consitent, inliers] = sprt_test(A, t, pts, pts_tilde, ...
            threshold, epsilon_i(test_num),delta_i(test_num), A_i(test_num));
        if ~pass % Model rejected
            % Re-estimate delta and keep epsilon the same
            % delta is the average # consistent in ALL rejected models 
            % incremental averaging results in
            n_rejected = n_rejected + 1;
            delta_est = delta_est*(n_rejected-1)/n_rejected + ...
                    n_consitent/(n_tested*n_rejected);
            if (delta_est > 1.05 * delta_i(test_num) || delta_est < 0.95 * delta_i(test_num))
            %if (abs(delta_est-delta_i(test_num))>0.05)
                updated_test = true;
                test_num = test_num + 1;
                delta_i(test_num) = delta_est;
                epsilon_i(test_num) = epsilon_i(test_num-1);
                A_i(test_num) = calculate_A(epsilon_i(test_num), delta_est);
                k_i(test_num) = 0;
            end
            continue;
        end
        % Model passed!

        if (n_consitent > max_inliers)
            max_inliers = n_consitent;
            inlier_set = inliers;
            A_best = A;
            t_best = t;

            % Estimate the inlier ratio using LO
            n_inliers_lo = local_opt(pts, pts_tilde, n_consitent, inliers, threshold);
            if n_inliers_lo/length(pts) > inlier_ratio
                % Found new inlier ratio with largest support so far
                % Design the next test
                updated_test = true;
                inlier_ratio = n_inliers_lo/length(pts);
                test_num = test_num + 1;
                delta_i(test_num) = delta_est;
                epsilon_i(test_num) = inlier_ratio;
                A_i(test_num) = calculate_A(inlier_ratio, delta_est);
                k_i(test_num) = 0;
            end
        end
        % Run max k iterations
        if updated_test
            updated_test = false;
            max_its_left = max_iters_sprt(inlier_ratio, epsilon_i, delta_i, k_i, A_i, test_num);
            if i> max_its_left
                break
            end
        end
    end
end

function [pass, n_tested, n_consitent, inliers] = sprt_test(A,t,pts, pts_tilde,threshold, epsilon, delta, A_decision)
    % Implementation of the SPRT test based on [chum, 2008].
    %
    % Returns: Model accepted or rejected, number of tested data points,
    % and the number of data points consistent with the model.
    % Note that log-likelihoods are used to prevent underflow.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    N = size(pts,2);
    ind = randsample(N, N);
    lhood_ratio = 0; n_tested = 0; n_consitent = 0; pass = true;
    inliers = false(1,N);
    
    for j=1:N
        n_tested = n_tested + 1;
        % Check if the jth point is an inlier
        % and compute likelihood ratio
        if (residual_lengths(A, t, pts(:,ind(j)), pts_tilde(:,ind(j))) <= threshold)
            % x_j = 1
            n_consitent = n_consitent + 1;
            inliers(ind(j)) = true;
            lhood_ratio = lhood_ratio + log(delta) - log(epsilon);
        else
            % x_j = 0
            lhood_ratio = lhood_ratio + log(1-delta) - log(1-epsilon);
            inliers(ind(j)) = false;
        end

        % Check with SPRT's decision threshold A
        if (lhood_ratio > log(A_decision))
            % Model rejected!
            pass = false;
            return;
        end
    end
    % Model passed since all points was tested without rejection.
end


function A_next = calculate_A(epsilon, delta)
% Calculate the SPRT decision threshold A
% Usually converges very fast
    t_m = 100; % (Time to estime affine trans.)/(time to evaluate 1 point)
    C = (1-delta)*log((1-delta)/(1-epsilon))+delta*log(delta/epsilon);
    A_0 = 1 + C*t_m;
    A_prev = A_0; A_next = -10000;
    while abs(A_next-A_prev) > 1e-3
        A_next = A_0 + log(A_prev);
        A_prev = A_next;
    end
end

function h_i = calc_h(eps, eps_i, del_i, l)
% Calculate h_i (eq. 14) using numerical method vpasolve
% TODO: add init table to increase robustness
    
    h_i = zeros(1, l);
    for i=1:l
        syms h
        eq = eps*(del_i/eps_i).^h + (1-eps)*((1-del_i)/(1-eps_i)).^h;
        % Don't include 0 to force non-zero solution.
        h_solution = vpasolve(eq==1, h, [0.001, 15]);
        if isempty(h_solution)
            h_solution = 0;
        end
        h_i(i) = h_solution;
    end
end

function k = max_iters_sprt(eps, eps_i, del_i, k_i, A_i, l)
% Get the number of iterations k_l (eq.16)
    eta_0 = 0.05;
    h_i = calc_h(eps, eps_i, del_i, l);
    P_g = eps.^3; % hardcoded 3...

    eta_l = 0;
    for i=1:l
        eta_l = eta_l + k_i(i)*log(1-P_g*(1-1/(A_i(i)^h_i(i))));
    end

    % The number of samples that needs to be drawn for l'th SPRT
    k = ceil((log(eta_0)-eta_l)/log(1-P_g/A_i(l)));
end

