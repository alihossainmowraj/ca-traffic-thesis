%% =========================================================
%  init4.m  -  DYNAMIC INITIALIZATION OF THE CA TRAFFIC SIMULATION
%  Requires dataset3.m to have been run (variables in the base workspace).
%  Produces: CA, P_var, h_hist, P_entry, D_upk, D_mk, veh_exited, v_mean
%% =========================================================

%% ---------------- STAGE 0: initialization ----------------
if ~exist('rho', 'var'),    rho    = 0.05; end     % initial density (veh/cell)
if ~exist('Q_flow', 'var'), Q_flow = 0.5;  end     % initial flow

dt     = params.dt;
dx     = params.dx;
g_safe = params.g_safe;
p_min  = params.p_min;
p_max  = params.p_max;
l_down = params.l_down;

v_mean  = round(Q_flow / rho);
i_all   = 0;
ptr_veh = 1;
n_t     = size(P_veh, 1);
P_veh(:, 9) = 0;                                   % reset: all vehicles inactive

P_var  = zeros(0, 18);
h_hist = zeros(0, 2);

R = cell(M, max(K_count));
V = cell(M, max(K_count));
A = cell(M, max(K_count));
G = cell(M, max(K_count));
S = cell(M, max(K_count));
n_placed_lane = zeros(M, max(K_count));

%% ---------------- STAGE 1: vehicle placement ----------------
for m = M:-1:1
    L = L_seg(m);
    for k = 1:K_count(m)

        % ---- 1A: initialize the lane ----
        R{m, k} = zeros(2, L);
        V{m, k} = zeros(1, L);
        A{m, k} = zeros(1, L);
        S{m, k} = zeros(4, L);
        G{m, k} = double(F{m, k} == 0);            % 0 on forbidden cells

        l_f  = sum(F{m, k});
        n_mk = ceil(rho * (L - l_f));

        % continuation of this lane (row found by column 1, not by row index)
        row_c = find(C{m}(:, 1) == k, 1);
        if isempty(row_c)
            mu = -1;  k_down = 0;
        else
            mu = C{m}(row_c, 3);  k_down = C{m}(row_c, 2);
        end

        if mu == -1
            x = L - l_down - g_safe;               % terminating lane: reserve downstream end
        else
            x = L - g_safe;
        end

        % ---- 1B: placement loop ----
        v_placed = 0;
        for veh = 1:n_mk

            % (b) vehicle at the pointer
            ID_r  = ptr_veh;
            assert(P_veh(ID_r, 9) == 0, ...
                'init4: ptr_veh reached an active vehicle (row %d): not enough vehicles in P_veh.', ID_r);
            ID    = P_veh(ID_r, 2);
            l_ID  = P_veh(ID_r, 3);
            v_max = P_veh(ID_r, 4);
            d_max = P_veh(ID_r, 6);
            alpha_ID = P_veh(ID_r, 7);
            tau_ID   = P_veh(ID_r, 8);

            % (c) feasible front-end position (scan upstream)
            lane_full = false;
            while true
                if x - l_ID < 0               % rear would leave the segment
                    lane_full = true;
                    break;
                end
                y = max(v_mean * dt, x - (2 * l_ID + g_safe));
                y = max(1, min(y, x - l_ID + 1));  % window always covers the footprint
                if sum(G{m, k}(y:x)) == x - y + 1
                    break;
                end
                x = x - 1;
            end
            if lane_full
                break;                             % (d)/(f): nothing placed, no index consumed
            end

            % (e) front and rear ends
            i_new = i_all + 1;
            x_i   = x;
            r_i   = x_i - l_ID + 1;

            % (g) feasible velocity over the occupied cells
            v_sv = min(v_max, v_mean);
            for j = r_i:x_i
                while v_sv >= M_spd{m, k}(j)
                    v_sv = v_sv - 1;
                    if v_sv <= 0
                        v_sv = 0;
                        break;
                    end
                end
            end
            v_i = max(0, floor(v_sv));

            % (h) insert into the CA matrices
            if ~(isnumeric(i_new) && isscalar(i_new) && isfinite(i_new)) || ...
               ~(isnumeric(ID) && isscalar(ID) && isfinite(ID))
                warning('init4:BadIndex', 'Non-numeric R entry: i_all=%s, m=%d, k=%d', ...
                    num2str(i_new), m, k);
            end
            R{m, k}(1, r_i:x_i) = i_new;
            R{m, k}(2, r_i:x_i) = ID;
            G{m, k}(r_i:x_i)    = 0;
            V{m, k}(r_i:x_i)    = v_i;
            A{m, k}(r_i:x_i)    = 0;
            S{m, k}(1, r_i:x_i) = 1;
            S{m, k}(2, r_i:x_i) = double(v_i > 0);
            S{m, k}(3, r_i:x_i) = 0;
            S{m, k}(4, r_i:x_i) = 0;

            % (i) downstream-zone flag
            f_down_i = double(x_i >= L - l_down);

            % (j) perception range
            if f_down_i == 1 && mu == -1
                p_i = L - x_i;
            else
                p_star = ceil(v_mean * tau_ID + v_mean^2 / (2 * d_max * alpha_ID)) + g_safe;
                p_i    = min(p_max, max(p_min, p_star));
            end

            % (k) initial heading gap
            h_temp      = 0;
            reached_end = true;
            range_curr  = min(x_i + p_i, L);
            for j = x_i + 1 : range_curr
                if G{m, k}(j) == 0
                    reached_end = false;
                    break;
                end
                h_temp = h_temp + 1;
            end
            if mu ~= -1 && reached_end
                range_down = max(x_i + p_i - L, 0);
                if range_down > 0
                    if m < M && k_down >= 1
                        for j = 1 : min(range_down, size(G{m + 1, k_down}, 2))
                            if G{m + 1, k_down}(j) == 0
                                break;
                            end
                            h_temp = h_temp + 1;
                        end
                    elseif m == M
                        h_temp = h_temp + range_down;     % exit buffer is empty
                    end
                end
            end
            h_i = min(p_max, h_temp);

            % (l) P_var row and (m) heading-gap history
            new_row     = zeros(1, 18);
            new_row(1)  = i_new;
            new_row(2)  = ID_r;
            new_row(3)  = v_i / dt;
            new_row(4)  = v_i;
            new_row(5)  = p_i;
            new_row(6)  = h_i;
            new_row(7)  = x_i;
            new_row(8)  = r_i;
            new_row(9)  = m;
            new_row(10) = k;
            new_row(15) = f_down_i;
            P_var  = [P_var;  new_row];
            h_hist = [h_hist; h_i, h_i];

            % (n) mark active, commit the index
            P_veh(ID_r, 9) = 1;
            i_all = i_new;

            % (o) rear gap required before the next vehicle
            next_ID_r = mod(ptr_veh, n_t) + 1;
            kappa_i   = floor(P_veh(next_ID_r, 4) * P_veh(next_ID_r, 7) * dt / g_safe) * dx;
            b_i       = g_safe + v_mean * P_veh(next_ID_r, 8) + kappa_i;

            % (p) position pointer, (q) database pointer, (r) counter
            x       = r_i - ceil(b_i);
            ptr_veh = next_ID_r;
            v_placed = v_placed + 1;

            % (s) lane-full check
            if x <= 1
                break;
            end
        end
        n_placed_lane(m, k) = v_placed;
    end
end

%% ---------------- STAGE 2: CA struct ----------------
CA          = struct();
CA.R        = R;
CA.V        = V;
CA.A        = A;
CA.G        = G;
CA.S        = S;
CA.F        = F(1:M_phys, :);
CA.C        = C;
CA.M_spd    = M_spd;
CA.K_count  = K_count;
CA.L_seg    = L_seg;

%% ---------------- STAGE 3: entry standby queue ----------------
K1      = K_count(1);
P_entry = zeros(K1, 4);
for k = 1:K1
    row = find(P_veh(:, 9) == 0, 1);               % first inactive row (see note below)
    % circular search starting at ptr_veh, as in the algorithm:
    order = [ptr_veh:n_t, 1:ptr_veh-1];
    idx   = find(P_veh(order, 9) == 0, 1);
    if isempty(idx)
        warning('init4:NoVehicle', 'No inactive vehicle available for P_entry lane %d.', k);
        continue;
    end
    row = order(idx);
    P_entry(k, 1) = row;
    P_entry(k, 2) = P_veh(row, 3);                 % length
    P_entry(k, 3) = P_veh(row, 4);                 % v_max
    P_entry(k, 4) = P_veh(row, 7);                 % alpha
    P_veh(row, 9) = 1;                             % reserved as standby
    ptr_veh       = mod(row, n_t) + 1;
end

%% ---------------- STAGE 4: upstream snapshot (columns [i, v, h, x]) ----------------
D_upk = cell(1, K1);
for k = 1:K1
    rows = find(P_var(:, 9) == 1 & P_var(:, 10) == k);
    if isempty(rows)
        D_upk{k} = zeros(0, 4);
    else
        D_upk{k} = [P_var(rows, 1), P_var(rows, 4), P_var(rows, 6), P_var(rows, 7)];
    end
end

%% ---------------- STAGES 5-6: archive and exit counter ----------------
D_mk       = zeros(0, 9);
veh_exited = 0;

%% ---------------- INVARIANT CHECK ----------------
n_placed = size(P_var, 1);
assert(isequal(P_var(:, 1), (1:n_placed)'), 'init4: P_var(:,1) is not 1..n.');
assert(numel(unique(P_var(:, 2))) == n_placed, 'init4: duplicate vehicle rows in P_var.');
assert(sum(P_veh(:, 9) == 1) == n_placed + nnz(P_entry(:, 1)), ...
    'init4: P_veh active flags differ from placed + queued vehicles.');
for i = 1:n_placed
    m = P_var(i, 9);  k = P_var(i, 10);
    r = P_var(i, 8);  x = P_var(i, 7);
    assert(all(CA.R{m, k}(1, r:x) == i) && nnz(CA.R{m, k}(1, :) == i) == x - r + 1, ...
        'init4: footprint of vehicle %d does not match CA.R.', i);
end

%% ---------------- STAGE 7: summary ----------------
fprintf('\n--- init4 summary ---\n');
fprintf('Physical segments    : %d\n', M_phys);
fprintf('Lanes per segment    : %s\n', mat2str(K_count));
fprintf('Segment lengths      : %s\n', mat2str(L_seg));
fprintf('Vehicles in P_veh    : %d\n', n_t);
fprintf('Vehicles placed      : %d\n', n_placed);
fprintf('P_var size           : %d x %d\n', size(P_var, 1), size(P_var, 2));
fprintf('P_entry size         : %d x %d\n', size(P_entry, 1), size(P_entry, 2));
fprintf('v_mean               : %d\n', v_mean);
fprintf('CA fields            : %s\n', strjoin(fieldnames(CA)', ', '));
fprintf('Placed per lane (rows = segments):\n'); disp(n_placed_lane);