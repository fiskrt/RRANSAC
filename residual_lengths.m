function res = residual_lengths(A, t, pts, pts_tilde)
    res = sum((A*pts+t-pts_tilde).^2, 1);
end

