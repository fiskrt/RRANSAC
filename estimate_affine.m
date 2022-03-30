function [A, t] = estimate_affine(pts, pts_tilde)
    % given three point corrospondences we wish to find the
    % affine transformation

    % TODO: make sure condition is not bad
    % normalize data, mean and std?
    
    M = [];
    b = [];
    for i=1:3
        M = [M;pts(:,i)' 0 0 1 0;
             0 0 pts(:,i)' 0 1];
        b = [b;pts_tilde(:, i)];
    end
    T = M\b;
    A = reshape(T(1:4), [2,2])';
    t = T(5:6);
end

