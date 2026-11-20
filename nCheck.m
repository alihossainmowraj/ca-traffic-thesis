% Merged into initDist.m - no longer needs to be run.
global chk
n_chk = 1000;
chk = struct();
chk.vehicle_position = NaN(n_chk, 3);
chk.shft_state       = NaN(n_chk, 3);
for n = 1:size(P_var, 1)
    ID_r = P_var(n, 2);
    if ID_r <= n_chk && P_var(n, 9) > 0
        chk.vehicle_position(ID_r, :) = P_var(n, [9 10 7]);
        chk.shft_state(ID_r, :)       = P_var(n, [12 16 13]);
    end
end
fprintf('chk built: %d vehicles recorded at t = 0\n', sum(~isnan(chk.vehicle_position(:, 1))));