%% =========================================================
%  runSim2.m  -  Cellular Automata Traffic Simulation (rewrite)
%
%  Run order:  new scenario (once): dataset  ->  initDist   (writes the .mat file)
%              every simulation:    clear all  ->  runSim2  (reads the .mat file)
%  Optional:   set T_max in the base workspace before calling runSim2
%              (default 30 steps).
%
%  Variable names, CA fields and the P_var column layout are unchanged:
%    P_var columns:
%     1 i        2 ID_r     3 a        4 v        5 p        6 h
%     7 x        8 r        9 m       10 k       11 omega   12 delta
%    13 sigma   14 f_stbl  15 f_down  16 f_lshift 17 f_sshift 18 t_last
%
%  Frame convention: (m, k, r, x) are given in the frame of the segment
%  that holds the REAR of the vehicle. When x > L_seg(m) the vehicle is in
%  transfer (f_sshift = 1) and its front part occupies cells 1..x-L of the
%  continuation lane in segment m+1 (or the exit buffer when m = M).
%
%  "Extended lane": a lane (m,k) is read beyond its end L through its
%  continuation. Beyond the end of a dropping lane (mu = -1) there is a
%  wall (forbidden, not free). Beyond the last segment there is the empty
%  exit buffer. Upstream of cell 1 is treated as free.
%% =========================================================

function runSim2()

global chk

% ----- Initial state file (written by initDist, never modified here) -----
init_file = 'initDist.mat';   % <-- input file for this run

if ~exist(init_file, 'file')
    error('runSim2: initial-state file ''%s'' not found. Run dataset and initDist first.', init_file);
end
S = load(init_file);

req_vars = {'CA', 'P_var', 'P_veh', 'P_entry', 'D_upk', 'D_mk', 'h_hist', 'v_mean', 'params'};
for n = 1:numel(req_vars)
    if ~isfield(S, req_vars{n})
        error('runSim2: variable ''%s'' is missing from ''%s''.', req_vars{n}, init_file);
    end
end

CA       = S.CA;
P_var    = S.P_var;
P_veh    = S.P_veh;
P_entry  = S.P_entry;
D_upk    = S.D_upk;
D_mk     = S.D_mk;
h_hist   = S.h_hist;
v_mean   = S.v_mean;
params   = S.params;
if isfield(S, 'chk')
    chk = S.chk;                % log starts from the loaded state
end

if isfield(S, 'init_info') && isfield(S.init_info, 'created')
    fprintf('Initial state loaded from %s (created %s)\n', init_file, S.init_info.created);
else
    fprintf('Initial state loaded from %s\n', init_file);
end

M = numel(CA.L_seg);

% =========================================================
%  SIMULATION CONTROL
% =========================================================
T_max = 90;
if evalin('base', 'exist(''T_max'', ''var'')')
    T_max = evalin('base', 'T_max');
end
check_inv   = true;     % consistency check of P_var vs CA.R after each step
veh_exited  = 0;
veh_entered = 0;

% ---- Continuation lookup by column 1 of C (C itself is not changed) ----
%  CA.k_down{m}(k) : downstream lane of lane k (0 if the lane drops)
%  CA.mu{m}(k)     : continuation flag of lane k (0 straight, -1 drop)
CA.k_down = cell(1, M);
CA.mu     = cell(1, M);
for m = 1:M
    K = CA.K_count(m);
    CA.k_down{m} = zeros(1, K);
    CA.mu{m}     = -ones(1, K);
    for k = 1:K
        row = find(CA.C{m}(:, 1) == k, 1);
        if ~isempty(row)
            CA.k_down{m}(k) = CA.C{m}(row, 2);
            CA.mu{m}(k)     = CA.C{m}(row, 3);
        end
    end
end
CA.n_conflict = 0;      % number of times a footprint had to overwrite another vehicle

fprintf('Starting simulation ...\n');
fprintf('  T_max            : %d\n', T_max);
fprintf('  Initial vehicles : %d\n', sum(P_var(:, 9) > 0));
fprintf('  Active in P_veh  : %d\n', sum(P_veh(:, 9) == 1));

% =========================================================
%  MAIN LOOP
% =========================================================
for t_n = 1:T_max

    % ---- Steps A-C: process segments from downstream to upstream ----
    for m = M:-1:1
        [CA, P_var, P_veh, h_hist, veh_exited] = processSegment(t_n, m, ...
            CA, P_var, P_veh, h_hist, params, veh_exited);
    end

    % ---- Rebuild upstream snapshot of segment 1 ----
    D_upk = buildUpstreamSnapshot(CA, P_var);

    % ---- Step D: entry queue ----
    [CA, P_veh, P_var, P_entry, h_hist, veh_entered] = processSegmentEntry(t_n, ...
        D_upk, veh_entered, CA, P_var, P_veh, P_entry, h_hist, params, v_mean);

    % ---- Consistency check ----
    if check_inv
        n_bad = checkInvariants(CA, P_var);
        if n_bad > 0
            warning('runSim2:invariant', 't = %d : %d vehicles inconsistent with CA.R', t_n, n_bad);
        end
    end

    % ---- Progress display ----
    if mod(t_n, 10) == 0 || t_n == T_max
        fprintf('  t = %4d | on road = %5d | exited = %5d | entered = %5d | mean v = %5.2f\n', ...
            t_n, sum(P_var(:, 9) > 0), veh_exited, veh_entered, mean(P_var(P_var(:, 9) > 0, 4)));
    end
end

% =========================================================
%  POST-SIMULATION
% =========================================================
fprintf('\n--- Simulation complete ---\n');
fprintf('  Steps             : %d\n', T_max);
fprintf('  Vehicles exited   : %d\n', veh_exited);
fprintf('  Vehicles entered  : %d\n', veh_entered);
fprintf('  Vehicles on road  : %d\n', sum(P_var(:, 9) > 0));
fprintf('  Final P_var rows  : %d\n', size(P_var, 1));
fprintf('  Footprint conflicts: %d\n', CA.n_conflict);

% ----- Write back to base workspace -----
assignin('base', 'CA',          CA);
assignin('base', 'P_var',       P_var);
assignin('base', 'P_veh',       P_veh);
assignin('base', 'P_entry',     P_entry);
assignin('base', 'D_upk',       D_upk);
assignin('base', 'D_mk',        D_mk);
assignin('base', 'h_hist',      h_hist);
assignin('base', 'veh_exited',  veh_exited);
assignin('base', 'veh_entered', veh_entered);

fprintf('\nState written back to base workspace.\n');

end


%% =========================================================
%  processSegment  -  all vehicles whose rear is in segment m
%% =========================================================
function [CA, P_var, P_veh, h_hist, veh_exited] = processSegment(t_n, m, ...
    CA, P_var, P_veh, h_hist, params, veh_exited)

    % Unprocessed vehicles of this segment, most downstream first
    % (ties: faster vehicle first). Positions of unprocessed vehicles do
    % not change during the step, so the order is fixed once.
    rows = find(P_var(:, 9) == m & P_var(:, 18) ~= t_n);
    if isempty(rows)
        return;
    end
    [~, ord] = sortrows([-P_var(rows, 7), -P_var(rows, 4)]);
    rows = rows(ord);

    for idx = 1:numel(rows)
        i = rows(idx);
        [CA, P_var, P_veh, h_hist, veh_exited] = processVehicle(i, t_n, ...
            CA, P_var, P_veh, h_hist, params, veh_exited);
    end
end


%% =========================================================
%  processVehicle  -  one complete update of subject vehicle i
%% =========================================================
function [CA, P_var, P_veh, h_hist, veh_exited] = processVehicle(i, t_n, ...
    CA, P_var, P_veh, h_hist, params, veh_exited)

    svData = svDataBuild(i, P_var, P_veh);

    ID_r          = svData.ID_r;
    v_prev        = svData.v_prev;
    x_prev        = svData.x_prev;
    r_prev        = svData.r_prev;
    m_prev        = svData.m_prev;
    k_prev        = svData.k_prev;
    omega_prev    = svData.omega_prev;
    delta_prev    = svData.delta_prev;
    f_stbl_prev   = svData.f_stbl_prev;
    f_lshift_prev = svData.f_lshift_prev;

    M = numel(CA.L_seg);

    % =========================================================
    %  STEP 1 - Clear old occupancy (only cells owned by i)
    % =========================================================
    CA = clearOccupancy(CA, i, m_prev, k_prev, r_prev, x_prev, delta_prev, f_lshift_prev);

    % =========================================================
    %  STEP 2 - Target position (source frame)
    % =========================================================
    [tp, h_hist] = getTargetPosition(CA, P_var, params, h_hist, svData, i);

    a_new = tp.a_new;
    v_new = tp.v_new;
    p_new = tp.p_new;
    h_new = tp.h_new;
    x_new = tp.x_star;
    r_new = tp.r_star;

    % =========================================================
    %  STEP 3 - Lane-shift state machine
    % =========================================================
    k_new        = k_prev;
    delta_new    = 0;
    sigma_new    = 0;
    f_lshift_new = 0;
    f_stbl_new   = max(0, f_stbl_prev - 1);
    shift_started = false;

    if f_lshift_prev == 1
        % ---- Shift in progress ----
        sigma_new = getSigma(svData, tp, params);
        if sigma_new <= 0
            % shift completes: vehicle now belongs to the target lane
            k_try = k_prev + delta_prev;
            if k_try >= 1 && k_try <= CA.K_count(m_prev)
                k_new = k_try;
            end
            sigma_new  = 0;
            f_stbl_new = params.T_stbl;
        else
            % shift continues (dual-lane occupancy)
            delta_new    = delta_prev;
            f_lshift_new = 1;
            f_stbl_new   = 0;
        end
    else
        % ---- Possibly initiate a shift (uses omega from previous step) ----
        can_initiate = (omega_prev == 2) || (omega_prev == 1 && f_stbl_prev == 0);
        if can_initiate
            delta_star = getShiftDirection(CA, params, svData, tp);
            if delta_star ~= 0 && ...
               checkLaneShiftFeasibility(CA, params, svData, tp, delta_star)
                delta_new     = delta_star;
                f_lshift_new  = 1;
                sigma_new     = params.dy;       % remaining lateral distance
                f_stbl_new    = 0;
                shift_started = true;
            end
        end
    end

    % =========================================================
    %  STEP 4 - Safety net: the new footprint must be free
    % =========================================================
    free = isFootprintFree(CA, m_prev, k_new, r_new, x_new, delta_new, f_lshift_new);

    if ~free && shift_started
        % cancel the new shift and try straight
        delta_new = 0;  f_lshift_new = 0;  sigma_new = 0;
        f_stbl_new = max(0, f_stbl_prev - 1);
        free = isFootprintFree(CA, m_prev, k_new, r_new, x_new, delta_new, f_lshift_new);
    end

    while ~free && x_new > x_prev
        x_new = x_new - 1;
        r_new = r_new - 1;
        free  = isFootprintFree(CA, m_prev, k_new, r_new, x_new, delta_new, f_lshift_new);
    end

    if x_new < tp.x_star
        v_new = max(0, x_new - x_prev);
        a_new = round((v_new - v_prev) / params.dt);
    end

    if ~free
        CA.n_conflict = CA.n_conflict + 1;   % will overwrite: reported at the end
    end

    % =========================================================
    %  STEP 5 - Segment transfer / exit
    % =========================================================
    m_new = m_prev;
    L     = CA.L_seg(m_prev);

    if r_new > L
        if m_prev == M
            % ---- Rear has left the last segment: vehicle exits ----
            P_veh(ID_r, 9) = 0;
            P_var(i, 6:18) = 0;
            veh_exited     = veh_exited + 1;
            % recordChk(ID_r, t_n, 0, 0, 0, 0, 0, 0);
            return;
        end

        k_down = CA.k_down{m_prev}(k_new);
        if k_down == 0 && f_lshift_new == 1
            % rear leaves a dropping lane during a shift: the shift completes here
            k_try = k_new + delta_new;
            if k_try >= 1 && k_try <= CA.K_count(m_prev) && CA.k_down{m_prev}(k_try) > 0
                k_new = k_try;
                delta_new = 0;  f_lshift_new = 0;  sigma_new = 0;
                f_stbl_new = params.T_stbl;
                k_down = CA.k_down{m_prev}(k_new);
            end
        end
        if k_down == 0
            % should not happen (wall beyond a dropping lane): hold position
            warning('runSim2:transfer', 'Vehicle %d tried to leave dropping lane m=%d k=%d', i, m_prev, k_new);
            x_new = x_prev;  r_new = r_prev;  v_new = 0;
            a_new = round(-v_prev / params.dt);
        else
            m_new = m_prev + 1;
            k_new = k_down;
            x_new = x_new - L;
            r_new = r_new - L;
            if f_lshift_new == 1
                k_tgt = k_new + delta_new;
                if k_tgt < 1 || k_tgt > CA.K_count(m_new)
                    delta_new = 0;  f_lshift_new = 0;  sigma_new = 0;
                end
            end
        end
    end

    % =========================================================
    %  STEP 6 - Occupy new position
    % =========================================================
    CA = occupyPosition(CA, i, ID_r, m_new, k_new, r_new, x_new, ...
                        delta_new, f_lshift_new, v_new, a_new);

    % =========================================================
    %  STEP 7 - Flags
    % =========================================================
    L_new        = CA.L_seg(m_new);
    f_sshift_new = double(x_new > L_new);
    f_down_new   = double(f_sshift_new == 0 && x_new >= L_new - params.l_down);

    % =========================================================
    %  STEP 8 - Lane-shift trigger for the next step
    % =========================================================
    if f_lshift_new == 0
        omega_new = checkLaneShiftTrigger(CA, params, svData, tp, ...
                        m_new, k_new, x_new, v_new, f_down_new, f_sshift_new);
    else
        omega_new = 0;
    end

    % =========================================================
    %  STEP 9 - Commit to P_var
    % =========================================================
    P_var = updatePvar(i, t_n, a_new, v_new, p_new, h_new, ...
        x_new, r_new, m_new, k_new, ...
        omega_new, delta_new, sigma_new, ...
        f_stbl_new, f_down_new, f_lshift_new, f_sshift_new, P_var);
end


%% =========================================================
%  processSegmentEntry  -  insert queued vehicles at the upstream boundary
%% =========================================================
function [CA, P_veh, P_var, P_entry, h_hist, veh_entered] = processSegmentEntry(t_n, ...
    D_upk, veh_entered, CA, P_var, P_veh, P_entry, h_hist, params, v_mean)

    K1      = CA.K_count(1);
    L1      = CA.L_seg(1);
    dt      = params.dt;
    g_safe  = params.g_safe;
    l_entry = params.p_min;            % clearance required behind the first vehicle
    l_down  = params.l_down;
    p_min   = params.p_min;
    p_max   = params.p_max;
    n_t     = size(P_veh, 1);

    for k = 1:K1

        if P_entry(k, 1) == 0
            continue;
        end

        ID_r  = P_entry(k, 1);
        l_ID  = P_entry(k, 2);
        v_max = P_entry(k, 3);
        alpha = P_entry(k, 4);

        % ---- (a) gap feasibility: free cells from cell 1 ----
        j_block = find(CA.G{1, k} == 0, 1);
        if isempty(j_block)
            g_entry = L1;
        else
            g_entry = j_block - 1;
        end
        if g_entry < l_ID + l_entry + g_safe
            continue;
        end

        % ---- (b) entry speed ----
        D_lane = D_upk{k};
        if isempty(D_lane)
            v_avg = v_mean;
        else
            v_avg = mean(D_lane(:, 2));
        end
        g_front = g_entry - l_ID;                         % free cells ahead of the front
        v_entry = round(min(max(v_avg, alpha * v_max), v_max));
        v_entry = min(v_entry, floor(min(CA.M_spd{1, k}(1:l_ID))));
        v_entry = max(0, min(v_entry, g_front - g_safe));

        % ---- (c) position: rear at cell 1 ----
        r_new = 1;
        x_new = l_ID;
        i_new = size(P_var, 1) + 1;

        CA = occupyPosition(CA, i_new, ID_r, 1, k, r_new, x_new, 0, 0, v_entry, 0);

        tau_ID   = P_veh(ID_r, 8);
        d_max    = P_veh(ID_r, 6);
        p_star   = ceil(v_entry * tau_ID + v_entry^2 / (2 * alpha * d_max)) + g_safe;
        p_i      = min(p_max, max(p_min, p_star));
        h_i      = min(p_i, g_front);

        new_row     = zeros(1, 18);
        new_row(1)  = i_new;
        new_row(2)  = ID_r;
        new_row(3)  = 0;
        new_row(4)  = v_entry;
        new_row(5)  = p_i;
        new_row(6)  = h_i;
        new_row(7)  = x_new;
        new_row(8)  = r_new;
        new_row(9)  = 1;
        new_row(10) = k;
        new_row(15) = double(x_new >= L1 - l_down);
        new_row(18) = t_n;
        P_var(i_new, :)  = new_row;
        h_hist(i_new, :) = [h_i, h_i];
        % recordChk(ID_r, t_n, 1, k, x_new, 0, 0, 0);

        P_veh(ID_r, 9) = 1;
        veh_entered    = veh_entered + 1;

        % ---- (d) refill P_entry(k): next inactive row after ID_r ----
        P_entry(k, :) = 0;
        for offset = 1:n_t
            row = mod(ID_r - 1 + offset, n_t) + 1;
            if P_veh(row, 9) == 0
                P_entry(k, 1) = row;
                P_entry(k, 2) = P_veh(row, 3);     % length
                P_entry(k, 3) = P_veh(row, 4);     % v_max
                P_entry(k, 4) = P_veh(row, 7);     % alpha
                P_veh(row, 9) = 1;                 % reserved as standby
                break;
            end
        end
        if P_entry(k, 1) == 0
            fprintf('processSegmentEntry: no inactive vehicle for lane %d.\n', k);
        end
    end
end


%% =========================================================
%  buildUpstreamSnapshot  -  D_upk{k} = [i, v, h, x], upstream first
%% =========================================================
function D_upk = buildUpstreamSnapshot(CA, P_var)

    K1    = CA.K_count(1);
    D_upk = cell(1, K1);
    for k = 1:K1
        row = CA.R{1, k}(1, :);
        nz  = row(row > 0);
        if isempty(nz)
            D_upk{k} = zeros(0, 4);
        else
            ids = unique(nz, 'stable');
            ids = ids(:);
            D_upk{k} = [ids, P_var(ids, 4), P_var(ids, 6), P_var(ids, 7)];
        end
    end
end


%% =========================================================
%  svDataBuild
%% =========================================================
function svData = svDataBuild(i, P_var, P_veh)
    svData = struct();
    ID_r = P_var(i, 2);
    svData.ID_r          = ID_r;
    svData.l_ID          = P_veh(ID_r, 3);
    svData.v_max         = P_veh(ID_r, 4);
    svData.a_max         = P_veh(ID_r, 5);
    svData.d_max         = P_veh(ID_r, 6);
    svData.alpha         = P_veh(ID_r, 7);
    svData.tau           = P_veh(ID_r, 8);
    svData.phi           = P_veh(ID_r, 9);
    svData.a_prev        = P_var(i, 3);
    svData.v_prev        = P_var(i, 4);
    svData.p_prev        = P_var(i, 5);
    svData.h_prev        = P_var(i, 6);
    svData.x_prev        = P_var(i, 7);
    svData.r_prev        = P_var(i, 8);
    svData.m_prev        = P_var(i, 9);
    svData.k_prev        = P_var(i, 10);
    svData.omega_prev    = P_var(i, 11);
    svData.delta_prev    = P_var(i, 12);
    svData.sigma_prev    = P_var(i, 13);
    svData.f_stbl_prev   = P_var(i, 14);
    svData.f_down_prev   = P_var(i, 15);
    svData.f_lshift_prev = P_var(i, 16);
    svData.f_sshift_prev = P_var(i, 17);
end


%% =========================================================
%  updatePvar
%% =========================================================
function P_var = updatePvar(i, t_n, ...
    a_new, v_new, p_new, h_new, ...
    x_new, r_new, m_new, k_new, ...
    omega_new, delta_new, sigma_new, ...
    f_stbl, f_down_new, f_lshift_new, f_sshift_new, P_var)

    P_var(i, 3)  = a_new;
    P_var(i, 4)  = v_new;
    P_var(i, 5)  = p_new;
    P_var(i, 6)  = h_new;
    P_var(i, 7)  = x_new;
    P_var(i, 8)  = r_new;
    P_var(i, 9)  = m_new;
    P_var(i, 10) = k_new;
    P_var(i, 11) = omega_new;
    P_var(i, 12) = delta_new;
    P_var(i, 13) = sigma_new;
    P_var(i, 14) = f_stbl;          % already decremented by the caller
    P_var(i, 15) = f_down_new;
    P_var(i, 16) = f_lshift_new;
    P_var(i, 17) = f_sshift_new;
    P_var(i, 18) = t_n;
    % recordChk(P_var(i, 2), t_n, m_new, k_new, x_new, delta_new, f_lshift_new, sigma_new);
end


%% =========================================================
%  recordChk  -  log position / shift state into global chk

% function recordChk(ID_r, t_n, m_new, k_new, x_new, delta_new, f_lshift_new, sigma_new)
%     global chk
%     if isempty(chk) || ID_r > size(chk.vehicle_position, 1)
%         return;
%     end
%     c_pos = 3 * t_n + (1:3);
%     c_sft = 3 * t_n + (1:3);
%     if size(chk.vehicle_position, 2) < c_pos(end)
%         chk.vehicle_position(:, end+1:c_pos(end)) = NaN;
%     end
%     if size(chk.shft_state, 2) < c_sft(end)
%         chk.shft_state(:, end+1:c_sft(end)) = NaN;
%     end
%     chk.vehicle_position(ID_r, c_pos) = [m_new, k_new, x_new];
%     chk.shft_state(ID_r, c_sft)       = [delta_new, f_lshift_new, sigma_new];
% end

%% =========================================================


%% =========================================================
%  getTargetPosition  -  acceleration, velocity and target cells
%% =========================================================
function [tp, h_hist] = getTargetPosition(CA, P_var, params, h_hist, svData, i)

    % ---- Static properties ----
    l_ID  = svData.l_ID;
    v_max = svData.v_max;
    a_max = svData.a_max;
    d_max = svData.d_max;
    alpha = svData.alpha;

    % ---- Dynamic state ----
    a_prev        = svData.a_prev;
    v_prev        = svData.v_prev;
    x_prev        = svData.x_prev;
    r_prev        = svData.r_prev;
    m_prev        = svData.m_prev;
    k_prev        = svData.k_prev;
    delta_prev    = svData.delta_prev;
    f_lshift_prev = svData.f_lshift_prev;

    dt     = params.dt;
    g_safe = params.g_safe;

    % =========================================================
    %  STEP 1 - Heading gap and perception range
    % =========================================================
    [h_new, p_new] = getHeadingGap(CA, params, x_prev, m_prev, k_prev, svData);

    k_tgt = k_prev + delta_prev;
    shift_lane_ok = (f_lshift_prev == 1 && k_tgt >= 1 && k_tgt <= CA.K_count(m_prev));
    if shift_lane_ok
        % While shifting, the vehicle is leaving its own lane: a static
        % obstruction there (lane end, forbidden cells) no longer stops it,
        % but other vehicles in either lane still do.
        [h_tgt, p_tgt] = getHeadingGap(CA, params, x_prev, m_prev, k_tgt, svData);
        p_new = min(p_new, p_tgt);
        h_own = scanGapVehicles(CA, m_prev, k_prev, x_prev, p_new);
        h_new = min([h_own, h_tgt, p_new]);
    end

    % =========================================================
    %  STEP 2 - Gap history and gap rate
    % =========================================================
    h_prev2      = h_hist(i, 2);
    h_hist(i, 2) = h_hist(i, 1);
    h_hist(i, 1) = h_new;
    dh = (h_new - h_prev2) / (2 * dt);

    % =========================================================
    %  STEP 3 - Braking tolerance
    % =========================================================
    if dh <= 0
        brk_tol = -min(alpha * v_prev, 0.4 * v_prev);
    else
        brk_tol = 0;
    end

    % =========================================================
    %  STEP 4 - Acceleration decision
    % =========================================================
    if alpha > 0.5
        if dh >= 0
            a_dec = 1;
        elseif dh >= brk_tol
            if (h_new - g_safe) > (v_prev + a_prev * dt * dt)
                a_dec = 0;
            else
                a_dec = -1;
            end
        else
            a_dec = -1;
        end
    else
        if dh > 0
            a_dec = 1;
        elseif dh > brk_tol
            if (h_new - g_safe) > (v_prev + a_prev * dt * dt)
                a_dec = 0;
            else
                a_dec = -1;
            end
        else
            a_dec = -1;
        end
    end

    % ---- Free road: nothing seen within the perception range ----
    %  (h is capped by p, and p shrinks with speed, so without this rule a
    %   shrinking p looks like a shrinking gap and drivers brake to a stop)
    if h_new >= p_new
        a_dec = 1;
    end

    % ---- Start-up rule: a stopped vehicle with free space ahead moves off ----
    if v_prev == 0 && dh >= 0 && (h_new - g_safe) > 0
        a_dec = 1;
    end

    % =========================================================
    %  STEP 5 - Cooperative deceleration (low-sensitivity drivers)
    % =========================================================
    if alpha < 0.5 && a_dec >= 0 && v_prev > 0
        if needsCooperation(CA, P_var, params, svData, p_new, i)
            a_dec = -2;
        end
    end

    % =========================================================
    %  STEP 6 - Applied acceleration
    % =========================================================
    if a_dec == 1
        a_star = ceil(a_prev + alpha * (a_max - max(0, a_prev)));
    elseif a_dec == 0
        a_star = a_prev;
    elseif a_dec == -1
        if dh == 0
            dh_eff = -1;
        else
            dh_eff = dh;
        end
        u_urg  = abs(1 - brk_tol / dh_eff);
        a_star = -min(d_max, ceil(alpha * d_max * u_urg));
    else   % a_dec == -2
        a_star = a_prev - ceil((1 - alpha) * abs(a_prev));
    end
    a_star = max(-d_max, min(a_max, a_star));

    % =========================================================
    %  STEP 7 - Feasible velocity
    % =========================================================
    v_star  = max(0, min(v_prev + a_star * dt, v_max));
    v_star  = round(v_star);

    % ---- Speed limits over occupied cells and the reach ----
    x_reach = x_prev + ceil(v_star * dt);
    if shift_lane_ok
        M_spd_min = min(extLane(CA, 'M_spd', m_prev, k_tgt, r_prev, x_reach));
    else
        M_spd_min = min(extLane(CA, 'M_spd', m_prev, k_prev, r_prev, x_reach));
    end
    v_star = min(v_star, floor(M_spd_min));

    % ---- Safety: never advance more than the free gap minus the safety gap ----
    v_star = min(v_star, max(0, h_new - g_safe));

    v_new = max(0, v_star);
    a_new = round((v_new - v_prev) / dt);

    % =========================================================
    %  STEP 8 - Target positions (source frame)
    % =========================================================
    if f_lshift_prev == 1
        v_ref     = alpha * v_max / 2;
        theta_deg = 22.5 * (v_ref / (v_new + v_ref));
        theta_rad = theta_deg * pi / 180;
        d_long    = v_new * dt * cos(theta_rad);
        x_star = round(x_prev + d_long);
        r_star = x_star - round(l_ID * cos(theta_rad)) + 1;
    else
        x_star = x_prev + v_new * dt;
        r_star = x_star - l_ID + 1;
    end

    % =========================================================
    %  STEP 9 - Output
    % =========================================================
    tp = struct();
    tp.p_new         = round(p_new);
    tp.h_new         = round(h_new);
    tp.dh            = dh;
    tp.a_dec         = a_dec;
    tp.a_new         = a_new;
    tp.v_new         = v_new;
    tp.x_star        = x_star;
    tp.r_star        = r_star;
    tp.f_sshift_star = double(x_star > CA.L_seg(m_prev));
end


%% =========================================================
%  needsCooperation  -  is a neighbour about to merge into my lane?
%  True when an adjacent lane is obstructed ahead (forbidden cells or
%  the end of a dropping lane within p_new) AND a vehicle in that lane
%  is roughly alongside (front within g_safe behind my front or ahead,
%  speed >= mine - 1). Only a moving driver cooperates: a stopped one
%  holding still next to a stopped merger would block both forever.
%% =========================================================
function coop = needsCooperation(CA, P_var, params, svData, p_new, i)

    coop   = false;
    m      = svData.m_prev;
    k      = svData.k_prev;
    x      = svData.x_prev;
    v      = svData.v_prev;
    g_safe = params.g_safe;
    p_min  = params.p_min;
    L      = CA.L_seg(m);

    for side = [-1, 1]
        kn = k + side;
        if kn < 1 || kn > CA.K_count(m)
            continue;
        end

        F_ahead = extLane(CA, 'F', m, kn, x + 1, x + p_new);
        if ~any(F_ahead == 1)
            continue;
        end

        ids = extLane(CA, 'R', m, kn, max(1, x - p_min), x);
        ids = unique(ids(ids > 0 & ids ~= i));
        for n = 1:numel(ids)
            id = ids(n);
            if P_var(id, 9) == m
                x_n = P_var(id, 7);
            elseif P_var(id, 9) == m + 1
                x_n = P_var(id, 7) + L;
            else
                continue;
            end
            if x_n >= x - g_safe && P_var(id, 4) >= v - 1
                coop = true;
                return;
            end
        end
    end
end


%% =========================================================
%  getHeadingGap  -  free cells ahead of x_sv in extended lane (m_sv,k_sv)
%% =========================================================
function [h_new, p_new] = getHeadingGap(CA, params, x_sv, m_sv, k_sv, svData)

    p_max  = params.p_max;
    p_min  = params.p_min;
    g_safe = params.g_safe;

    v_prev = svData.v_prev;
    p_star = ceil(v_prev * svData.tau + (v_prev^2) / (2 * svData.alpha * svData.d_max)) + g_safe;
    p_new  = min(p_max, max(p_min, p_star));

    h_new = scanGap(CA, m_sv, k_sv, x_sv, p_new);
end


%% =========================================================
%  scanGap  -  number of consecutive free cells after x (up to p)
%% =========================================================
function h = scanGap(CA, m, k, x, p)
    g = extLane(CA, 'G', m, k, x + 1, x + p);
    j = find(g == 0, 1);
    if isempty(j)
        h = p;
    else
        h = j - 1;
    end
end


%% =========================================================
%  scanGapVehicles  -  like scanGap, but only other vehicles block
%% =========================================================
function h = scanGapVehicles(CA, m, k, x, p)
    ids = extLane(CA, 'R', m, k, x + 1, x + p);
    j   = find(ids > 0, 1);
    if isempty(j)
        h = p;
    else
        h = j - 1;
    end
end


%% =========================================================
%  getSigma  -  remaining lateral distance of a lane shift
%% =========================================================
function sigma_new = getSigma(svData, tp, params)

    v_new      = tp.v_new;
    sigma_prev = svData.sigma_prev;
    v_max      = svData.v_max;
    alpha      = svData.alpha;
    dt         = params.dt;

    if v_new > 0
        v_ref     = alpha * v_max / 2;
        theta_deg = 22.5 * (v_ref / (v_new + v_ref));   % angle decreases with speed
        theta_rad = theta_deg * pi / 180;
        d_lat     = v_new * dt * sin(theta_rad);
        sigma_new = max(0, sigma_prev - d_lat);
    else
        sigma_new = sigma_prev;
    end
end


%% =========================================================
%  getShiftDirection  -  delta = -1 (left), +1 (right), 0 (none)
%% =========================================================
function delta_new = getShiftDirection(CA, params, svData, tp)

    m_sv = svData.m_prev;
    k_sv = svData.k_prev;
    x_sv = tp.x_star;
    v_new = tp.v_new;

    p_star   = ceil(v_new * svData.tau + (v_new^2) / (2 * svData.alpha * svData.d_max)) + params.g_safe;
    p_latest = min(params.p_max, max(params.p_min, p_star));

    K_m = CA.K_count(m_sv);
    L   = CA.L_seg(m_sv);

    % ---- Boundary lanes: only one direction exists ----
    if k_sv == 1
        delta_new = 1;
        return;
    elseif k_sv == K_m
        delta_new = -1;
        return;
    end

    k_left  = k_sv - 1;
    k_right = k_sv + 1;

    % ---- Near the downstream edge: avoid dropping lanes ----
    if x_sv >= L - params.l_down && x_sv <= L
        mu_left  = CA.mu{m_sv}(k_left);
        mu_right = CA.mu{m_sv}(k_right);
        if mu_left == -1 && mu_right ~= -1
            delta_new = 1;   return;
        elseif mu_left ~= -1 && mu_right == -1
            delta_new = -1;  return;
        elseif mu_left == -1 && mu_right == -1
            delta_new = 0;   return;
        end
    end

    % ---- Obstructions (forbidden cells / walls) ahead in each lane ----
    %  Road geometry is known to the driver, so static obstructions are
    %  looked for over p_max (not only the speed-dependent range).
    F_left  = extLane(CA, 'F', m_sv, k_left,  x_sv + 1, x_sv + params.p_max);
    F_right = extLane(CA, 'F', m_sv, k_right, x_sv + 1, x_sv + params.p_max);
    j_left  = find(F_left  == 1, 1);
    j_right = find(F_right == 1, 1);
    left_blocked  = ~isempty(j_left);
    right_blocked = ~isempty(j_right);

    if left_blocked && ~right_blocked
        delta_new = 1;   return;
    elseif ~left_blocked && right_blocked
        delta_new = -1;  return;
    elseif left_blocked && right_blocked
        if j_right > j_left
            delta_new = 1;
        elseif j_left > j_right
            delta_new = -1;
        else
            delta_new = 0;
        end
        return;
    end

    % ---- Both free: choose the larger heading gap (right on tie) ----
    h_left  = scanGap(CA, m_sv, k_left,  x_sv, p_latest);
    h_right = scanGap(CA, m_sv, k_right, x_sv, p_latest);
    if h_right >= h_left
        delta_new = 1;
    else
        delta_new = -1;
    end
end


%% =========================================================
%  checkLaneShiftTrigger  -  omega for the next step
%  0 = no reason to shift, 1 = optional shift, 2 = mandatory shift
%% =========================================================
function omega_new = checkLaneShiftTrigger(CA, params, svData, tp, ...
    m_new, k_new, x_new, v_new, f_down_new, f_sshift_new)

    l_ID   = svData.l_ID;
    g_safe = params.g_safe;
    p_new  = tp.p_new;
    L      = CA.L_seg(m_new);
    K_m    = CA.K_count(m_new);

    omega_new = 0;

    % ---- CHECK 1: heading gap compressing and a neighbour lane offers more ----
    if tp.dh < 0
        h_own = scanGap(CA, m_new, k_new, x_new, p_new);
        h_adj = -Inf;
        if k_new - 1 >= 1
            h_adj = max(h_adj, scanGap(CA, m_new, k_new - 1, x_new, p_new));
        end
        if k_new + 1 <= K_m
            h_adj = max(h_adj, scanGap(CA, m_new, k_new + 1, x_new, p_new));
        end
        if h_adj > h_own
            omega_new = max(omega_new, 1);
        end
    end

    % ---- CHECK 2: own lane terminates at the downstream edge ----
    if f_down_new == 1 && f_sshift_new == 0 && CA.mu{m_new}(k_new) == -1
        if L - x_new - g_safe > 4 * v_new
            omega_new = max(omega_new, 1);
        else
            omega_new = 2;
        end
    end

    % ---- CHECK 3: static obstruction (forbidden cells) ahead ----
    p_star = ceil(v_new * svData.tau + (v_new^2) / (2 * svData.alpha * svData.d_max)) + g_safe;
    p_obs  = min(params.p_max, max(params.p_min, p_star));
    F_ahead = extLane(CA, 'F', m_new, k_new, x_new + 1, x_new + p_obs);
    j_obs   = find(F_ahead == 1, 1);
    if ~isempty(j_obs)
        d_obs = j_obs - 1;
        if d_obs > 4 * v_new - l_ID
            omega_new = max(omega_new, 1);
        else
            omega_new = 2;
        end
    end
end


%% =========================================================
%  checkLaneShiftFeasibility  -  is lane k_prev+delta_star free from
%  (r_prev - g_safe) to x_star (extended lane, so a transfer is covered)?
%% =========================================================
function feasible = checkLaneShiftFeasibility(CA, params, svData, tp, delta_star)

    feasible = false;
    m  = svData.m_prev;
    kt = svData.k_prev + delta_star;
    if kt < 1 || kt > CA.K_count(m)
        return;
    end

    % (1) target cells free, from g_safe behind my rear to my new front
    j1 = svData.r_prev - params.g_safe;
    j2 = tp.x_star;
    g  = extLane(CA, 'G', m, kt, j1, j2);
    if ~all(g == 1)
        return;
    end

    % (2) no static obstruction (forbidden cells / lane end) in the target
    %     lane within p_max ahead: never shift into a dead end
    F = extLane(CA, 'F', m, kt, j2 + 1, j2 + params.p_max);
    feasible = ~any(F == 1);
end


%% =========================================================
%  isFootprintFree  -  all cells of the planned footprint are free
%% =========================================================
function free = isFootprintFree(CA, m, k, r, x, delta, f_lshift)
    free  = true;
    lanes = footprintLanes(CA, m, k, delta, f_lshift);
    for n = 1:numel(lanes)
        x_end = footprintEnd(CA, m, lanes(n), r, x, f_lshift);
        g = extLane(CA, 'G', m, lanes(n), r, x_end);
        if ~all(g == 1)
            free = false;
            return;
        end
    end
end


%% =========================================================
%  footprintEnd  -  last cell of the footprint in one lane.
%  A shifting vehicle (dual occupancy) only occupies its lanes up to the
%  first static obstruction (forbidden cell / lane end) in each lane.
%% =========================================================
function x_end = footprintEnd(CA, m, k, r, x, f_lshift)
    x_end = x;
    if f_lshift == 1
        F = extLane(CA, 'F', m, k, r, x);
        j = find(F == 1, 1);
        if ~isempty(j)
            x_end = r + j - 2;
        end
    end
end


%% =========================================================
%  footprintLanes  -  lanes of segment m occupied by a vehicle
%% =========================================================
function lanes = footprintLanes(CA, m, k, delta, f_lshift)
    lanes = k;
    if f_lshift == 1
        kt = k + delta;
        if kt >= 1 && kt <= CA.K_count(m)
            lanes = [k, kt];
        end
    end
end


%% =========================================================
%  occupyPosition  -  write vehicle i into CA (source + downstream part)
%  S rows: 1 occupied, 2 moving, 3 lane-shift side (-1 source, +1 target),
%          4 transfer part (-1 source segment, +1 downstream segment)
%% =========================================================
function CA = occupyPosition(CA, i, ID, m, k, r, x, delta, f_lshift, v, a)

    M     = numel(CA.L_seg);
    L     = CA.L_seg(m);
    lanes = footprintLanes(CA, m, k, delta, f_lshift);
    in_tr = double(x > L);

    for n = 1:numel(lanes)
        kk = lanes(n);
        if numel(lanes) == 2
            s3 = 2 * (kk ~= k) - 1;          % -1 source lane, +1 target lane
        else
            s3 = 0;
        end
        x_end = footprintEnd(CA, m, kk, r, x, f_lshift);

        % ---- part inside segment m ----
        j1 = max(1, r);
        j2 = min(x_end, L);
        if j2 >= j1
            jj = j1:j2;
            CA.R{m, kk}(1, jj) = i;
            CA.R{m, kk}(2, jj) = ID;
            CA.V{m, kk}(jj)    = v;
            CA.A{m, kk}(jj)    = a;
            CA.G{m, kk}(jj)    = 0;
            CA.S{m, kk}(:, jj) = repmat([1; double(v > 0); s3; -in_tr], 1, numel(jj));
        end

        % ---- part in the downstream segment (exit buffer is not stored) ----
        if x_end > L && m < M
            kd = CA.k_down{m}(kk);
            if kd > 0 && CA.mu{m}(kk) ~= -1
                j2d = min(x_end - L, CA.L_seg(m + 1));
                j1d = max(1, r - L);
                if j2d >= j1d
                    jj = j1d:j2d;
                    CA.R{m + 1, kd}(1, jj) = i;
                    CA.R{m + 1, kd}(2, jj) = ID;
                    CA.V{m + 1, kd}(jj)    = v;
                    CA.A{m + 1, kd}(jj)    = a;
                    CA.G{m + 1, kd}(jj)    = 0;
                    CA.S{m + 1, kd}(:, jj) = repmat([1; double(v > 0); s3; 1], 1, numel(jj));
                end
            end
        end
    end
end


%% =========================================================
%  clearOccupancy  -  erase vehicle i (only cells it owns)
%% =========================================================
function CA = clearOccupancy(CA, i, m, k, r, x, delta, f_lshift)

    M     = numel(CA.L_seg);
    L     = CA.L_seg(m);
    lanes = footprintLanes(CA, m, k, delta, f_lshift);

    for n = 1:numel(lanes)
        kk = lanes(n);
        CA = clearCells(CA, i, m, kk, max(1, r), min(x, L));
        if x > L && m < M
            kd = CA.k_down{m}(kk);
            if kd > 0
                CA = clearCells(CA, i, m + 1, kd, max(1, r - L), min(x - L, CA.L_seg(m + 1)));
            end
        end
    end
end


function CA = clearCells(CA, i, m, k, j1, j2)
    if j2 < j1
        return;
    end
    jj  = j1:j2;
    own = jj(CA.R{m, k}(1, jj) == i);      % membership check
    if isempty(own)
        return;
    end
    CA.R{m, k}(:, own) = 0;
    CA.V{m, k}(own)    = 0;
    CA.A{m, k}(own)    = 0;
    CA.S{m, k}(:, own) = 0;
    CA.G{m, k}(own)    = double(CA.F{m, k}(own) == 0);
end


%% =========================================================
%  extLane  -  read field fld of the extended lane (m,k), cells j1..j2
%  fld: 'G' (1 free), 'F' (1 forbidden), 'R' (owner index), 'M_spd'
%% =========================================================
function out = extLane(CA, fld, m, k, j1, j2)

    if j2 < j1
        out = zeros(1, 0);
        return;
    end
    j   = j1:j2;
    n   = numel(j);
    M   = numel(CA.L_seg);
    L   = CA.L_seg(m);

    % values for: [free cell, wall]
    switch fld
        case 'G',     v_free = 1;   v_wall = 0;
        case 'F',     v_free = 0;   v_wall = 1;
        case 'R',     v_free = 0;   v_wall = 0;
        case 'M_spd', v_free = Inf; v_wall = Inf;
    end
    out = v_free * ones(1, n);

    % ---- inside segment m ----
    in_seg = (j >= 1) & (j <= L);
    if any(in_seg)
        out(in_seg) = readField(CA, fld, m, k, j(in_seg));
    end

    % ---- beyond the end of segment m ----
    beyond = (j > L);
    if ~any(beyond)
        return;
    end
    jb = j(beyond) - L;

    if m == M
        % exit buffer: empty, with its own speed limit
        if strcmp(fld, 'M_spd')
            kx  = CA.k_down{m}(k);
            if kx < 1, kx = k; end
            spd = CA.M_spd{M + 1, kx};
            out(beyond) = spd(min(jb, numel(spd)));
        end
    elseif CA.mu{m}(k) == -1 || CA.k_down{m}(k) == 0
        out(beyond) = v_wall;                 % lane drop: wall
    else
        kd   = CA.k_down{m}(k);
        L2   = CA.L_seg(m + 1);
        vals = v_free * ones(1, numel(jb));
        ok   = jb <= L2;
        if any(ok)
            vals(ok) = readField(CA, fld, m + 1, kd, jb(ok));
        end
        out(beyond) = vals;
    end
end


function vals = readField(CA, fld, m, k, jj)
    switch fld
        case 'G',     vals = CA.G{m, k}(jj);
        case 'F',     vals = CA.F{m, k}(jj);
        case 'R',     vals = CA.R{m, k}(1, jj);
        case 'M_spd', vals = CA.M_spd{m, k}(jj);
    end
end


%% =========================================================
%  checkInvariants  -  every vehicle on the road owns exactly its cells
%% =========================================================
function n_bad = checkInvariants(CA, P_var)
    n_bad = 0;
    on_road = find(P_var(:, 9) > 0);
    for n = 1:numel(on_road)
        i = on_road(n);
        m = P_var(i, 9);  k = P_var(i, 10);
        r = P_var(i, 8);  x = P_var(i, 7);
        x_end = footprintEnd(CA, m, k, r, x, P_var(i, 16));
        j1 = max(1, r);  j2 = min(x_end, CA.L_seg(m));
        if j2 >= j1 && ~all(CA.R{m, k}(1, j1:j2) == i)
            n_bad = n_bad + 1;
        end
    end
end