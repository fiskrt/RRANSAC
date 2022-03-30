function [A, t] = least_squares_affine(pts, pts_tilde)
    M = [];
    b = [];
    for i=1:length(pts)
        M = [M;pts(:,i)' 0 0 1 0;
             0 0 pts(:,i)' 0 1];
        b = [b;pts_tilde(:, i)];
    end
    T = M\b;
    A = reshape(T(1:4), [2,2])';
    t = T(5:6);
end

