%% =========================================================
%  SIM8.M  —  Cellular Automata Traffic Simulation

function runSim()

% ----- Globals shared with local functions -----
global CA P_var P_veh P_entry D_upk D_mk h_hist params v_mean


% ----- Fetch from base workspace -----
CA       = evalin('base', 'CA');
P_var    = evalin('base', 'P_var');
P_veh    = evalin('base', 'P_veh');
P_entry  = evalin('base', 'P_entry');
D_upk    = evalin('base', 'D_upk');
D_mk     = evalin('base', 'D_mk');
h_hist   = evalin('base', 'h_hist');
v_mean   = evalin('base', 'v_mean');
params   = evalin('base', 'params');
% S_all    = evalin('base', 'S_all');
% K_count  = CA.K_count;
% L_seg    = CA.L_seg;
% M_spd    = CA.M_spd;
% F        = CA.F;
M        = size(CA.K_count, 2);

% =========================================================
%  SIMULATION CONTROL
% =========================================================
T_max      = 30;
veh_exited = 0;
veh_entered  = 0;

fprintf('Starting simulation ...\n');
fprintf('  T_max            : %d\n', T_max);
fprintf('  Initial vehicles : %d\n', size(P_var, 1));
fprintf('  Active in P_veh  : %d\n', sum(P_veh(:, 9) == 1));

% =========================================================
%  MAIN LOOP
% =========================================================
for t_n = 1:T_max

    % ---- Step A: Process most downstream segment (m = M) ----
    [CA, P_var, h_hist, veh_exited] = processSegmentLast(t_n, veh_exited, CA, P_var, P_veh, h_hist, params);

    % ---- Step B: Process interior segments (m = M-1 .. 2) ----
    for m = M-1 : -1 : 2
        [CA, P_var, h_hist] = processSegmentInterior(t_n, CA, P_var, P_veh, h_hist, params, m);
    end

    % ---- Step C: Process first segment (m = 1), rebuild D_upk ----
    [CA, P_var, h_hist, D_upk] = processSegmentFirst(t_n, CA, P_var, P_veh, h_hist, params);

    % ---- Step D: Process entry queue ----
    [CA, P_veh, P_var, P_entry, h_hist, veh_entered] = processSegmentEntry(t_n, ...
    D_upk, veh_entered, CA, P_var, P_veh, P_entry, h_hist, params, v_mean);

    % ---- Progress display ----
    if mod(t_n, 10) == 0
        fprintf('  t = %4d | active = %5d | exited = %5d\n  | entered = %5d', ...
            t_n, sum(P_veh(:, 9) == 1), veh_exited, veh_entered);
    end
end

% =========================================================
%  POST-SIMULATION
% =========================================================
fprintf('\n--- Simulation complete ---\n');
fprintf('  Steps            : %d\n', T_max);
fprintf('  Vehicles exited  : %d\n', veh_exited);
fprintf('  Vehicles entered  : %d\n', veh_entered);
fprintf('  Final P_var rows : %d\n', size(P_var, 1));
fprintf('  D_mk rows        : %d\n', size(D_mk, 1));


% ----- Write back to base workspace -----
assignin('base', 'CA',         CA);
assignin('base', 'P_var',      P_var);
assignin('base', 'P_veh',      P_veh);
assignin('base', 'P_entry',    P_entry);
assignin('base', 'D_upk',      D_upk);
assignin('base', 'D_mk',       D_mk);
assignin('base', 'h_hist',     h_hist);
assignin('base', 'veh_exited', veh_exited);

fprintf('\nState written back to base workspace.\n');

end


%% processSegmentLast
function [CA, P_var, h_hist, veh_exited] = processSegmentLast(t_n, veh_exited, CA, P_var, P_veh, h_hist, params)

    m = size(CA.L_seg, 2);
    L = CA.L_seg(m);

    while true

        % ---- Select next subject vehicle in segment M ----
        [i, svData] = getSubjectVehicle(m, CA, P_var, P_veh, t_n);
        if i == 0
            break;
        end

        % ---- Extract state ----
        ID_r          = svData.ID_r;
        x_prev        = svData.x_prev;
        r_prev        = svData.r_prev;
        k_prev        = svData.k_prev;
        delta_prev    = svData.delta_prev;
        f_stbl_prev   = svData.f_stbl_prev;
        f_down_prev   = svData.f_down_prev;
        f_lshift_prev = svData.f_lshift_prev;
        f_sshift_prev = svData.f_sshift_prev;

        % =====================================================
        %  CASE 1 — Interior vehicle: delegate to upstream handler
        % =====================================================
        if f_down_prev == 0 && f_sshift_prev == 0
            [CA, P_var, h_hist] = processVehicleUpstream(CA, P_var, h_hist, params, svData, i, t_n);
            continue;
        end

        % =====================================================
        %  CASES 2–5: Downstream-edge or exiting vehicle
        % =====================================================

        % ---- Compute target position ----
        [tp, h_hist] = getTargetPosition(CA, P_var, params, h_hist, svData, i);
        v_new  = tp.v_new;
        a_new  = tp.a_new;
        p_new  = tp.p_new;
        h_new  = tp.h_new;
        x_star = tp.x_star;
        r_star = tp.r_star;

        % ---- Clear old occupancy (cap at L) ----
        x_clear_end = min(x_prev, L);
        CA = clearOccupancy(CA, m, k_prev, r_prev, x_clear_end);

        if f_lshift_prev == 1
            k_adj_prev = k_prev + delta_prev;
            if k_adj_prev >= 1 && k_adj_prev <= CA.K_count(m)
                CA = clearOccupancy(CA, m, k_adj_prev, r_prev, x_clear_end);
            end
        end

        % =====================================================
        %  SUB-CASE 3 / 5-COMPLETE — Full exit: rear-end crosses L
        % =====================================================
        %  This covers:
        %    3. Instantaneous exit   (f_sshift_prev == 0, r_star > L)
        %    5. Exit completion      (f_sshift_prev == 1, r_star > L)
        if r_star > L
            P_veh(ID_r, 9) = 0;          % deactivate in database
            P_var(i, 6:18) = 0;          % clear dynamic fields, keep i
            veh_exited     = veh_exited + 1;
            continue;
        end

        % =====================================================
        %  SUB-CASE 2 / 4 / 5-CONTINUE — Partial or no crossing
        % =====================================================

        % ---- Determine exiting status for this step ----
        if f_sshift_prev == 1
            f_sshift_new = 1;            % already exiting: continue
        elseif x_star > L
            f_sshift_new = 1;            % front crossed: initiate exit
        else
            f_sshift_new = 0;            % stays fully inside M
        end

        % ---- Clamp positions for CA storage ----
        r_new = min(r_star, L);          % rear stays in M
        if r_new < 1
            r_new = 1;
        end
        x_new = min(x_star, L);          % front capped at L

        % ---- Downstream flag ----
        if f_sshift_new == 1
            f_down_new = 0;              % exiting, no longer edge-bound
        else
            f_down_new = double(x_new >= L - params.l_down);
        end

        % ---- Lane-shift handling ----
        if f_lshift_prev == 1

            sigma_new = getSigma(svData, tp, pa);

            if sigma_new <= 0
                % ---- Shift finalizes ----
                k_new = k_prev + delta_prev;
                if k_new < 1 || k_new > CA.K_count(m)
                    k_new = k_prev;
                end
                delta_new    = 0;
                f_lshift_new = 0;
                sigma_new    = 0;
                f_stbl_end   = params.T_stbl;
                omega_new    = 0;

                % occupyPosition(CA, P_var, tp, m_curr, k_curr, j_start, j_end, i)
                CA = occupyPosition(CA, P_var, tp, m, k_new, r_new, x_new, i);
                for j = r_new:x_new
                    CA.S{m, k_new}(:, j) = ...
                        [1; double(v_new>0); 0; -f_sshift_new];
                end
            else
                % ---- Shift continues (dual lane) ----
                k_new        = k_prev;
                delta_new    = delta_prev;
                f_lshift_new = 1;
                f_stbl_end   = 0;
                omega_new    = 0;
                k_tgt        = k_prev + delta_prev;

                % occupyPosition(CA, P_var, tp, m_curr, k_curr, j_start, j_end, i)
                CA = occupyPosition(CA, P_var, tp, m, k_prev, r_new, x_new, i);
                if k_tgt >= 1 && k_tgt <= CA.K_count(m)
                    CA = occupyPosition(CA, P_var, tp, m, k_tgt, r_new, x_new, i);
                    for j = r_new:x_new
                        CA.S{m, k_prev}(:, j) = [1; double(v_new>0); -1; -f_sshift_new];
                        CA.S{m, k_tgt }(:, j) = [1; double(v_new>0); +1; -f_sshift_new];
                    end
                else
                    for j = r_new:x_new
                        CA.S{m, k_prev}(:, j) = ...
                            [1; double(v_new>0); 0; -f_sshift_new];
                    end
                    f_lshift_new = 0;
                    delta_new    = 0;
                end
            end

        else
            % ---- No active shift → straight progression ----
            k_new        = k_prev;
            omega_new    = 0;
            delta_new    = 0;
            sigma_new    = 0;
            f_lshift_new = 0;
            f_stbl_end   = max(0, f_stbl_prev - 1);

            % occupyPosition(CA, P_var, tp, m_curr, k_curr, j_start, j_end, i)
            CA = occupyPosition(CA, P_var, tp, m, k_new, r_new, x_new, i);
            for j = r_new:x_new
                CA.S{m, k_new}(:, j) = [1; double(v_new>0); 0; -f_sshift_new];
            end
        end

        % =====================================================
        %  Commit to P_var
        % =====================================================
        P_var = updatePvar(i, t_n, a_new, v_new, p_new, h_new, ...
            x_new, r_new, m, k_new, ...
            omega_new, delta_new, sigma_new, ...
            f_stbl_end, f_down_new, f_lshift_new, f_sshift_new, P_var);
    end
end


%%  processSegmentInterior
function [CA, P_var, h_hist] = processSegmentInterior(t_n, CA, P_var, P_veh, h_hist, params, m)
    
    while true

        % ---- Select next subject vehicle in segment m ----
        [i, svData] = getSubjectVehicle(m, CA, P_var, P_veh, t_n);

        % ---- No more unprocessed vehicles in this segment ----
        if i == 0
            break;              % break: main loop of this function
        end

        % ---- Dispatch based on vehicle state ----
        if svData.f_down_prev == 1 || svData.f_sshift_prev == 1
            % Downstream-edge or extended-transfer vehicle
            [CA, P_var, h_hist] = processVehicleDownstream(CA, P_var, h_hist, params, svData, i, t_n);
        else
            % Interior vehicle (straight / lane-shift only)
            [CA, P_var, h_hist] = processVehicleUpstream(CA, P_var, h_hist, params, svData, i, t_n);
        end

    end  % end of main loop

end


%%  processSegmentFirst
function [CA, P_var, h_hist, D_upk] = processSegmentFirst(t_n, CA, P_var, P_veh, h_hist, params)

    m = 1;

    % =====================================================
    %  MAIN LOOP — process all vehicles in segment 1
    % =====================================================
    while true

        % ---- Select next subject vehicle in segment 1 ----
        [i, svData] = getSubjectVehicle(m, CA, P_var, P_veh, t_n);

        % ---- No more unprocessed vehicles in this segment ----
        if i == 0
            break;              % break: main loop of this function
        end

        % ---- Dispatch based on vehicle state ----
        if svData.f_down_prev == 1 || svData.f_sshift_prev == 1
            % Downstream-edge or extended-transfer vehicle
            [CA, P_var, h_hist] = processVehicleDownstream(CA, P_var, h_hist, params, svData, i, t_n);
        else
            % Interior vehicle (straight / lane-shift only)
            [CA, P_var, h_hist] = processVehicleUpstream(CA, P_var, h_hist, params, svData, i, t_n);
        end

        % Loop repeats; the handler has already set P_var(i, 18) = t_n,
        % so this vehicle will not be re-selected in the same time step.

    end  % end of main loop (segment 1 vehicle processing)

    % =====================================================
    %  REBUILD D_upk  —  per-lane snapshot of segment 1
    % =====================================================

    %  INITIATE D_upk  —  one field per lane of segment 1
    D_upk = struct();
    k1 = CA.K_count(1);
    for k = 1:k1
        D_upk.(sprintf('k_%d', k)) = zeros(0, 4);
    end
    
    
    %  POPULATE D_upk  —  scan each lane from upstream to downstream
    for k = 1:k1
    
        L = size(CA.R{1, k}, 2);
        fld = sprintf('k_%d', k);
    
        j = 1;
        while j <= L
    
            i_here = CA.R{1, k}(1, j);
    
            if i_here ~= 0
                v_i = P_var(i_here, 4);
                h_i = P_var(i_here, 6);
                x_i = P_var(i_here, 7);
    
                D_upk.(fld) = [D_upk.(fld); i_here, v_i, h_i, x_i];
    
                j = x_i + h_i;      % skip vehicle + known empty run
            else
                j = j + 1;
            end
        end

     end  % end of lane loop

end


%%  processSegmentEntry
function [CA, P_veh, P_var, P_entry, h_hist, veh_entered] = processSegmentEntry(t_n, ...
    D_upk, veh_entered, CA, P_var, P_veh, P_entry, h_hist, params, v_mean)


    %  CONSTANTS
    K_count = CA.K_count;
    K1      = K_count(1);
    L1      = CA.L_seg(1);

    dt      = params.dt;
    g_safe  = params.g_safe;
    l_entry = params.p_min;
    l_down  = params.l_down;

    n_t = size(P_veh, 1);

    % =====================================================
    %  LANE LOOP — for each lane of segment 1
    % =====================================================
    for k = 1:K1

        % ---- Skip lane if no vehicle is queued ----
        if P_entry(k, 1) == 0
            continue;               % continue: lane loop
        end

        ID_r  = P_entry(k, 1);
        l_ID  = P_entry(k, 2);
        v_max = P_entry(k, 3);
        alpha = P_entry(k, 4);

        % =================================================
        %  (a) ENTRY FEASIBILITY
        % =================================================

        % ---- Extract the lane's snapshot ----
        fld    = sprintf('k_%d', k);
        D_lane = D_upk.(fld);

        % ---- Gap feasibility ----
        %  If the lane has any vehicle, the most upstream one (first row)
        %  has its rear-end at r_first = x_first - l_first + 1. The gap
        %  from cell 1 to the rear-end is therefore r_first - 1.
        %  If the lane is empty, the entire segment length is free.
        if isempty(D_lane)
            g_entry = L1;
        else
            i_first = D_lane(1, 1);
            x_first = D_lane(1, 4);
            IDr_fst = P_var(i_first, 2);
            l_first = P_veh(IDr_fst, 3);
            r_first = x_first - l_first + 1;
            g_entry = r_first - 1;
        end

        if g_entry < l_ID + l_entry + g_safe
            continue;               % continue: lane loop — insufficient gap
        end

        % ---- Speed feasibility ----
        if isempty(D_lane)
            v_avg = v_mean;
        else
            v_avg = mean(D_lane(:, 2));
        end

        v_entry = max(v_avg, alpha * v_max);
        v_entry = round( min(v_entry, v_max) );

        if v_entry > g_entry - g_safe
            continue;               % continue: lane loop — would overshoot gap
        end

        % =================================================
        %  (b) ENTRY POSITION
        % =================================================
        x_new = round(v_entry * dt) - l_entry;
        r_new = x_new - l_ID + 1;

        if r_new < 1
            continue;               % continue: lane loop — vehicle doesn't fit
        end

        % ---- Verify target cells are free ----
        cells_free = true;
        for j = r_new:x_new
            if CA.R{1, k}(1, j) ~= 0
                cells_free = false;
                break;              % break: free-cell check loop
            end
        end

        if ~cells_free
            continue;               % continue: lane loop — occupancy conflict
        end

        % =================================================
        %  (c) COMMIT ENTRY
        % =================================================

        % ---- Assign a new dynamic index ----
        i_new = size(P_var, 1) + 1;

        % ---- Build target-position struct for occupyPosition ----
        tp = struct();
        tp.v_new         = v_entry;
        tp.a_new         = round(v_entry / dt);
        tp.p_new         = g_entry + g_safe;
        tp.h_new         = g_entry;
        tp.x_star        = x_new;
        tp.r_star        = r_new;
        tp.m_star        = 1;
        tp.f_sshift_star = 0;

        % ---- Occupy CA cells ----
        % occupyPosition(CA, P_var, tp, m_curr, k_curr, j_start, j_end, i)
        CA = occupyPosition(CA, P_var, tp, 1, k, r_new, x_new, i_new);
        for j = r_new:x_new
            CA.S{1, k}(:, j) = [1; (v_entry > 0); 0; 0];
        end

        % ---- Append a new row to P_var ----
        P_var(i_new, :)  = zeros(1, 18);
        P_var(i_new,  1) = i_new;
        P_var(i_new,  2) = ID_r;
        P_var(i_new,  3) = round(v_entry / dt);
        P_var(i_new,  4) = v_entry;
        P_var(i_new,  5) = g_entry + g_safe;
        P_var(i_new,  6) = g_entry;
        P_var(i_new,  7) = x_new;
        P_var(i_new,  8) = r_new;
        P_var(i_new,  9) = 1;
        P_var(i_new, 10) = k;
        P_var(i_new, 11) = 0;
        P_var(i_new, 12) = 0;
        P_var(i_new, 13) = 0;
        P_var(i_new, 14) = 0;
        P_var(i_new, 15) = double(x_new >= L1 - l_down);
        P_var(i_new, 16) = 0;
        P_var(i_new, 17) = 0;
        P_var(i_new, 18) = t_n;

        % ---- Extend heading-gap history ----
        h_hist(i_new, :) = [g_entry, g_entry];

        % ---- Mark vehicle as active ----
        P_veh(ID_r, 9) = 1;
        veh_entered = veh_entered + 1;

        % =================================================
        %  (d) REFILL P_entry(k)
        % =================================================
        P_entry(k, :) = 0;

        i_ptr = mod(i_new, n_t) + 1;
        found = false;

        for offset = 0:n_t-1
            r_star = mod(i_ptr - 1 + offset, n_t) + 1;
            if P_veh(r_star, 9) == 0
                P_entry(k, 1) = r_star;
                P_entry(k, 2) = P_veh(r_star, 3);      % l_ID
                P_entry(k, 3) = P_veh(r_star, 4);      % v_max
                P_entry(k, 4) = P_veh(r_star, 7);      % alpha
                P_veh(r_star, 9) = 1;
                found = true;
                break;              % break: circular scan loop
            end
        end

        if ~found
            fprintf('processSegmentEntry: no inactive vehicle for lane %d.\n', k);
        end

    end  % end of lane loop

end

%%  processVehicleUpstream
function [CA, P_var, h_hist] = processVehicleUpstream(CA, P_var, h_hist, params, svData, i, t_n)

    % =========================================================
    %  STEP 1 — Extract state
    % =========================================================
    ID_r          = svData.ID_r;
    alpha         = svData.alpha;
    v_max         = svData.v_max;
    l_ID          = svData.l_ID;
    a_prev        = svData.a_prev;
    v_prev        = svData.v_prev;
    x_prev        = svData.x_prev;
    r_prev        = svData.r_prev;
    k_prev        = svData.k_prev;
    m_prev        = svData.m_prev;
    omega_prev    = svData.omega_prev;
    delta_prev    = svData.delta_prev;
    sigma_prev    = svData.sigma_prev;
    f_stbl_prev   = svData.f_stbl_prev;
    f_lshift_prev = svData.f_lshift_prev;

    L      = CA.L_seg(m_prev);
    K_m    = CA.K_count(m_prev);
    l_down = params.l_down;
    M_spd  = CA.M_spd;

    % =========================================================
    %  STEP 2 — Defaults for _new
    % =========================================================
    m_new        = m_prev;
    k_new        = k_prev;
    omega_new    = omega_prev;
    delta_new    = delta_prev;
    sigma_new    = sigma_prev;
    f_stbl_new   = max(0, f_stbl_prev - 1);
    f_down_new   = 0;
    f_lshift_new = f_lshift_prev;
    f_sshift_new = 0;

    % =========================================================
    %  STEP 3 — Clear old occupancy
    % =========================================================
    CA = clearOccupancy(CA, m_prev, k_prev, r_prev, x_prev);

    if f_lshift_prev == 1
        k_adj_prev = k_prev + delta_prev;
        if k_adj_prev >= 1 && k_adj_prev <= K_m
            CA = clearOccupancy(CA, m_prev, k_adj_prev, r_prev, x_prev);
        end
    end

    % =========================================================
    %  STEP 4 — Compute target position
    % =========================================================
    [tp, h_hist] = getTargetPosition(CA, P_var, params, h_hist, svData, i);

    a_new = tp.a_new;
    v_new = tp.v_new;
    p_new = tp.p_new;
    h_new = tp.h_new;
    x_new = round(tp.x_star);
    r_new = round(tp.r_star);

    if x_new >= L - l_down
        f_down_new = 1;
    end

    % =========================================================
    %  STEP 5 — Determine shift status
    %  Case A: an active shift is in progress
    %  Case B: no active shift; possibly initiate one
    % =========================================================

    if f_lshift_prev == 1

        % -----------------------------------------------------
        %  CASE A — active shift
        % -----------------------------------------------------
        sigma_new = getSigma(svData, tp, params);

        if sigma_new <= 0

            % ---- Shift completes ----
            k_new = k_prev + delta_prev;
            if k_new < 1 || k_new > K_m
                k_new = k_prev;
            end
            delta_new    = 0;
            sigma_new    = 0;
            f_lshift_new = 0;
            f_stbl_new   = params.T_stbl;
            omega_new    = omega_prev;

            CA = occupyPosition(CA, P_var, tp, m_prev, k_new, r_new, x_new, i);
            for j = r_new : x_new                     % LOOP: state write (finalized)
                CA.S{m_prev, k_new}(:, j) = [1; double(v_new>0); 0; 0];
            end
        else

            % ---- Shift continues (dual lane) ----
            k_new        = k_prev;
            delta_new    = delta_prev;
            f_lshift_new = 1;
            f_stbl_new   = 0;
            omega_new    = 0;
            k_target        = k_prev + delta_prev;

            % occupyPosition(CA, P_var, tp, m_curr, k_curr, j_start, j_end, i)
            CA = occupyPosition(CA, P_var, tp, m_prev, k_prev, r_new, x_new, i);
            if k_target >= 1 && k_target <= K_m
                CA = occupyPosition(CA, P_var, tp, m_prev, k_target, r_new, x_new, i);
                for j = r_new : x_new                 % LOOP: state write (dual lane)
                    CA.S{m_prev, k_prev}(:, j) = [1; double(v_new>0); -1; 0];
                    CA.S{m_prev, k_target }(:, j) = [1; double(v_new>0);  1; 0];
                end
            else
                % Fallback: target invalid → release source
                for j = r_new : x_new                 % LOOP: state write (fallback)
                    CA.S{m_prev, k_prev}(:, j) = [1; double(v_new>0); 0; 0];
                end
                f_lshift_new = 0;
                delta_new    = 0;
            end
        end

    else

        % -----------------------------------------------------
        %  CASE B — no active shift; consider initiating one
        % -----------------------------------------------------
        k_new           = k_prev;
        omega_star      = checkLaneShiftTrigger(svData, tp, CA, params, i, ...
                                                h_hist, f_down_new);
        omega_new       = max(omega_star, omega_prev);
        can_initiate    = (omega_new == 2) || ...
                          (omega_new == 1 && f_stbl_prev == 0);

        delta_star   = 0;
        f_lshift_new = 0;
        k_target        = k_new;
        sigma_new    = 0;

        if omega_new ~= 0 && can_initiate

            delta_star = getShiftDirection(CA, params, svData, tp);

            if delta_star ~= 0
                feasible = checkLaneShiftFeasibility(CA, params, svData, tp);

                if feasible
                    k_target = k_new + delta_star;
                    if k_target >= 1 && k_target <= K_m
                        sigma_new    = params.dy;         % lateral units
                        f_lshift_new = 1;
                        f_stbl_new   = 0;
                        delta_new    = delta_star;
                    else
                        delta_star    = 0;
                        k_target        = k_new;
                        delta_new    = 0;
                    end
                end
            end
        end

        % ---- Occupy source lane (always) ----
        CA = occupyPosition(CA, P_var, tp, m_prev, k_prev, r_new, x_new, i);

        if f_lshift_new == 1
            % ---- Occupancy: dual lane ----
            % occupyPosition(CA, P_var, tp, m_curr, k_curr, j_start, j_end, i)
            CA = occupyPosition(CA, P_var, tp, m_prev, k_target, r_new, x_new, i);
            for j = r_new : x_new                     % LOOP: state write (dual lane initiation)
                CA.S{m_prev, k_prev}(:, j) = [1; double(v_new>0); -1; 0];
                CA.S{m_prev, k_target }(:, j) = [1; double(v_new>0);  1; 0];
            end
            k_new = k_prev;
        else
            % ---- Occupancy: straight ----
            for j = r_new : x_new                     % LOOP: state write (straight)
                CA.S{m_prev, k_prev}(:, j) = [1; double(v_new>0); 0; 0];
            end
            delta_new = 0;
            sigma_new = 0;
            % Re-evaluate omega on the straight-progressed state
            omega_new = checkLaneShiftTrigger(svData, tp, CA, params, i, ...
                                              h_hist, f_down_new);
        end
    end

    % =========================================================
    %  STEP 8 — Commit to P_var
    % =========================================================
    P_var = updatePvar(i, t_n, a_new, v_new, p_new, h_new, ...
                       x_new, r_new, m_new, k_new, ...
                       omega_new, delta_new, sigma_new, ...
                       f_stbl_new, f_down_new, f_lshift_new, f_sshift_new, P_var);

end

%% processVehicleDownstream
function [CA, P_var, h_hist] = processVehicleDownstream(CA, P_var, h_hist, params, svData, i, t_n)

    % =========================================================
    %  EXTRACT ALL VEHICLE DATA
    % =========================================================
    ID_r         = svData.ID_r;
    alpha        = svData.alpha;
    v_max        = svData.v_max;

    a_prev       = svData.a_prev;
    v_prev       = svData.v_prev;
    p_prev       = svData.p_prev;
    h_prev       = svData.h_prev;
    x_prev       = svData.x_prev;
    r_prev       = svData.r_prev;
    k_prev       = svData.k_prev;
    m_prev       = svData.m_prev;
    omega_prev   = svData.omega_prev;
    delta_prev   = svData.delta_prev;
    sigma_prev   = svData.sigma_prev;
    f_stbl_prev  = svData.f_stbl_prev;
    f_down_prev  = svData.f_down_prev;
    f_lshift_prev = svData.f_lshift_prev;
    f_sshift_prev = svData.f_sshift_prev;

    % =========================================================
    %  TARGET VALUES FROM tp
    % =========================================================
    [tp, h_hist] = getTargetPosition(CA, P_var, params, h_hist, svData, i);
    x_star = tp.x_star;
    r_star = tp.r_star;
    v_new  = tp.v_new;
    a_new  = tp.a_new;
    p_new  = tp.p_new;
    h_new  = tp.h_new;

    L = size(CA.R{m_prev, k_prev}, 2);

    % =========================================================
    %  DEFAULTS FOR _new VARIABLES
    % =========================================================
    % m_new        = m_prev;
    % k_new        = k_prev;
    % omega_new    = omega_prev;
    delta_new    = delta_prev;
    sigma_new    = sigma_prev;
    f_stbl       = max(0, f_stbl_prev - 1);
    f_lshift_new = f_lshift_prev;
    % f_sshift_new = f_sshift_prev;
    % f_down_new   = 1;

    % =========================================================
    %  CASE 1: VEHICLE ALREADY IN EXTENDED TRANSFER (f_sshift_prev == 1)
    % =========================================================
    if f_sshift_prev == 1

        % ---- Determine downstream lane ----
        k_down = 0;
        if k_prev <= size(CA.C{m_prev}, 1)
            k_down = CA.C{m_prev}(k_prev, 2);
        end
        if k_down == 0
            warning('processDownstream: invalid k_down for m=%d k=%d', m_prev, k_prev);
            return;
        end

        % ---- Clear source-segment occupancy (rear part) ----
        
        CA = clearOccupancy(CA, m_prev, k_prev, r_prev, L);
        if f_lshift_prev == 1
            k_adj_prev = k_prev + delta_prev;
            if k_adj_prev >= 1 && k_adj_prev <= CA.K_count(m_prev)
                CA = clearOccupancy(CA, m_prev, k_adj_prev, r_prev, L);
            end
        end

        % ---- Clear downstream-segment occupancy (front part) ----
        CA = clearOccupancy(CA, m_prev+1, k_down, 1, x_prev);
        if f_lshift_prev == 1
            k_down_adj = k_down + delta_prev;
            if k_down_adj >= 1 && k_down_adj <= CA.K_count(m_prev+1)
                CA = clearOccupancy(CA, m_prev+1, k_down_adj, 1, x_prev);
            end
        end

        if r_star > L
            % =====================================================
            %  SUBCASE 1a: TRANSFER COMPLETES
            % =====================================================
            m_new        = m_prev + 1;
            k_new        = k_down;
            x_new        = x_star - L;
            r_new        = r_star - L;
            f_down_new   = double(x_new >= size(CA.R{m_new, k_new}, 2) - params.l_down);
            f_sshift_new = 0;
            omega_new    = omega_prev;

            if f_lshift_prev == 1
                sigma_new = getSigma(svData, tp, params);

                if sigma_new <= 0
                    % ---- shift finalizes ----
                    k_new = k_down + delta_prev;
                    if k_new < 1 || k_new > CA.K_count(m_new)
                        k_new = k_down;
                    end
                    delta_new    = 0;
                    f_lshift_new = 0;
                    sigma_new    = 0;
                    f_stbl       = params.T_stbl;

                    CA = occupyPosition(CA, P_var, tp, m_new, k_new, r_new, x_new, i);
                    for j = r_new:x_new
                        CA.S{m_new, k_new}(:, j) = [1; double(v_new>0); 0; 0];
                    end
                else
                    % ---- shift continues in downstream segment ----
                    k_new        = k_down;
                    delta_new    = delta_prev;
                    f_lshift_new = 1;
                    f_stbl       = 0;

                    k_target = k_new + delta_new;
                    CA = occupyPosition(CA, P_var, tp, m_new, k_new, r_new, x_new, i);
                    if k_target >= 1 && k_target <= CA.K_count(m_new)
                        CA = occupyPosition(CA, P_var, tp, m_new, k_target, r_new, x_new, i);
                        for j = r_new:x_new
                            CA.S{m_new, k_new}(:, j)    = [1; double(v_new>0); -1; 0];
                            CA.S{m_new, k_target}(:, j) = [1; double(v_new>0);  1; 0];
                        end
                    else
                        for j = r_new:x_new
                            CA.S{m_new, k_new}(:, j) = [1; double(v_new>0); 0; 0];
                        end
                        f_lshift_new = 0; delta_new = 0;
                    end
                end

            else
                % ---- no lane shift active → straight in downstream ----
                k_new = k_down;
                CA = occupyPosition(CA, P_var, tp, m_new, k_new, r_new, x_new, i);
                for j = r_new:x_new
                    CA.S{m_new, k_new}(:, j) = [1; double(v_new>0); 0; 0];
                end
            end

        else
            % =====================================================
            %  SUBCASE 1b: TRANSFER CONTINUES
            % =====================================================
            m_new        = m_prev;
            k_new        = k_prev;
            x_new        = x_star;
            r_new        = r_star;
            f_down_new   = 0;
            f_sshift_new = 1;
            omega_new    = omega_prev;

            if f_lshift_prev == 1
                sigma_new = getSigma(svData, tp, params);

                if sigma_new <= 0
                    k_new = k_prev + delta_prev;
                    if k_new < 1 || k_new > CA.K_count(m_prev)
                        k_new = k_prev;
                    end
                    delta_new    = 0;
                    f_lshift_new = 0;
                    sigma_new    = 0;
                    f_stbl       = params.T_stbl;

                    CA = occupyPosition(CA, P_var, tp, m_prev, k_new, r_new, L, i);
                    for j = r_new:L
                        CA.S{m_prev, k_new}(:, j) = [1; double(v_new>0); 0; -1];
                    end
                    if x_new > L
                        CA = occupyPosition(CA, P_var, tp, m_prev+1, k_down, 1, x_new-L, i);
                        for j = 1:(x_new-L)
                            CA.S{m_prev+1, k_down}(:, j) = [1; double(v_new>0); 0; 1];
                        end
                    end
                else
                    % keep both lanes occupied
                    k_new        = k_prev;
                    delta_new    = delta_prev;
                    f_lshift_new = 1;
                    f_stbl       = 0;
                    k_target     = k_prev + delta_prev;

                    CA = occupyPosition(CA, P_var, tp, m_prev, k_prev, r_new, L, i);
                    if k_target >= 1 && k_target <= CA.K_count(m_prev)
                        CA = occupyPosition(CA, P_var, tp, m_prev, k_target, r_new, L, i);
                        for j = r_new:L
                            CA.S{m_prev, k_prev}(:, j)   = [1; double(v_new>0); -1; -1];
                            CA.S{m_prev, k_target}(:, j) = [1; double(v_new>0);  1; -1];
                        end
                    end
                    if x_new > L
                        CA = occupyPosition(CA, P_var, tp, m_prev+1, k_down, 1, x_new-L, i);
                        for j = 1:(x_new-L)
                            CA.S{m_prev+1, k_down}(:, j) = [1; double(v_new>0); -1; 1];
                        end
                        if f_lshift_new == 1 && ...
                           (k_down + delta_new) >= 1 && (k_down + delta_new) <= CA.K_count(m_prev+1)
                            CA = occupyPosition(CA, P_var, tp, m_prev+1, k_down+delta_new, 1, x_new-L, i);
                            for j = 1:(x_new-L)
                                CA.S{m_prev+1, k_down+delta_new}(:, j) = [1; double(v_new>0); 1; 1];
                            end
                        end
                    end
                end

            else % case 1 & f_lshift = 0
                k_new = k_prev;
                CA = occupyPosition(CA, P_var, tp, m_prev, k_new, r_new, L, i);
                for j = r_new:L
                    CA.S{m_prev, k_new}(:, j) = [1; double(v_new>0); 0; -1];
                end
                if x_new > L
                    CA = occupyPosition(CA, P_var, tp, m_prev+1, k_down, 1, x_new-L, i);
                    for j = 1:(x_new-L)
                        CA.S{m_prev+1, k_down}(:, j) = [1; double(v_new>0); 0; 1];
                    end
                end

                % ---- Lane shift initiation using omega_prev ----
                omega_new    = checkLaneShiftTrigger(svData, tp, CA, params, i, ...
                               h_hist, f_down_new);
                can_initiate = (omega_new == 2) || (omega_new == 1 && f_stbl_prev == 0);

                if omega_new ~= 0 
                    if can_initiate
                        delta_new = getShiftDirection(CA, params, svData, tp);
    
                        if delta_new ~= 0
                            feasible = checkLaneShiftFeasibility(svData, tp, CA, params);
    
                            if feasible
                                k_target  = k_new + delta_new;
                                v_ref     = alpha * v_max / 2;
                                theta_deg = 22.5 * (v_ref / (v_new + v_ref));
                                theta_rad = theta_deg * pi / 180;
                                sigma_new = params.dy / tan(theta_rad);
    
                                f_lshift_new = 1;
                                f_stbl       = 0;
    
                                if k_target >= 1 && k_target <= CA.K_count(m_new)
                                    CA = occupyPosition(CA, P_var, tp, m_new, k_target, r_new, x_new, i);
                                    for j = r_new:x_new
                                        CA.S{m_new, k_new}(:, j)    = [1; double(v_new>0); -1; 0];
                                        CA.S{m_new, k_target}(:, j) = [1; double(v_new>0);  1; 0];
                                    end
                                else
                                    f_lshift_new = 0; 
                                    sigma_new = 0;
                                end
                            end
                        end
                    end
                else % if omega_new = 0
                    omega_new = checkLaneShiftTrigger(svData, tp, CA, params, i, ...
                                              h_hist, f_down_new);
                end
            end
        end

    % =========================================================
    %  CASE 2: VEHICLE NOT YET IN TRANSFER (f_sshift_prev == 0)
    % =========================================================
    else  % f_sshift_prev == 0

        % ---- Clear source-segment occupancy ----
        
        CA = clearOccupancy(CA, m_prev, k_prev, r_prev, L);
        if f_lshift_prev == 1
            k_adj_prev = k_prev + delta_prev;
            if k_adj_prev >= 1 && k_adj_prev <= CA.K_count(m_prev)
                CA = clearOccupancy(CA, m_prev, k_adj_prev, r_prev, L);
            end
        end

        if r_star > L
            % =====================================================
            %  SUBCASE 2a: INSTANTANEOUS TRANSFER / EXIT
            % =====================================================
            m_new        = m_prev + 1;
            k_new        = CA.C{m_prev}(k_prev, 2);
            x_new        = x_star - L;
            r_new        = r_star - L;
            f_down_new   = double(x_new >= size(CA.R{m_new, k_new}, 2) - params.l_down);
            f_sshift_new = 0;
            omega_new    = omega_prev;

            if f_lshift_prev == 1
                sigma_new = getSigma(svData, tp, params);
                if sigma_new <= 0
                    k_new = k_new + delta_prev;
                    if k_new < 1 || k_new > CA.K_count(m_new)
                        k_new = CA.C{m_prev}(k_prev, 2);
                    end
                    delta_new    = 0;
                    f_lshift_new = 0;
                    sigma_new    = 0;
                    f_stbl       = params.T_stbl;

                    CA = occupyPosition(CA, P_var, tp, m_new, k_new, r_new, x_new, i);
                    for j = r_new:x_new
                        CA.S{m_new, k_new}(:, j) = [1; double(v_new>0); 0; 0];
                    end
                else
                    k_new        = CA.C{m_prev}(k_prev, 2);
                    delta_new    = delta_prev;
                    f_lshift_new = 1;
                    f_stbl       = 0;

                    CA = occupyPosition(CA, P_var, tp, m_new, k_new, r_new, x_new, i);
                    CA = occupyPosition(CA, P_var, tp, m_new, k_new+delta_new, r_new, x_new, i);
                    for j = r_new:x_new
                        CA.S{m_new, k_new}(:, j) = [1; double(v_new>0); -1; 0];
                    end
                    for j = r_new:x_new
                        CA.S{m_new, k_new+delta_new}(:, j) = [1; double(v_new>0); 1; 0];
                    end
                end
            else  % case 2 & f_lshift = 0

                % ---- Straight in downstream ----
                k_new = CA.C{m_prev}(k_prev, 2);
                CA = occupyPosition(CA, P_var, tp, m_new, k_new, r_new, x_new, i);
                for j = r_new:x_new
                    CA.S{m_new, k_new}(:, j) = [1; double(v_new>0); 0; 0];
                end

                % ---- Lane shift initiation using omega_prev ----
                omega_new    = checkLaneShiftTrigger(svData, tp, CA, params, i, ...
                                              h_hist, f_down_new);
                can_initiate = (omega_new == 2) || (omega_new == 1 && f_stbl_prev == 0);

                if omega_new ~= 0 
                    if can_initiate
                        delta_new = getShiftDirection(CA, params, svData, tp);
    
                        if delta_new ~= 0
                            feasible = checkLaneShiftFeasibility(svData, tp, CA, params);
    
                            if feasible
                                k_target  = k_new + delta_new;
                                v_ref     = alpha * v_max / 2;
                                theta_deg = 22.5 * (v_ref / (v_new + v_ref));
                                theta_rad = theta_deg * pi / 180;
                                sigma_new = params.dy / tan(theta_rad);
    
                                f_lshift_new = 1;
                                f_stbl       = 0;
    
                                if k_target >= 1 && k_target <= CA.K_count(m_new)
                                    CA = occupyPosition(CA, P_var, tp, m_new, k_target, r_new, x_new, i);
                                    for j = r_new:x_new
                                        CA.S{m_new, k_new}(:, j)    = [1; double(v_new>0); -1; 0];
                                        CA.S{m_new, k_target}(:, j) = [1; double(v_new>0);  1; 0];
                                    end
                                else
                                    f_lshift_new = 0; 
                                    sigma_new = 0;
                                end
                            end
                        end
                    end
                else
                    omega_new = checkLaneShiftTrigger(svData, tp, CA, params, i, ...
                                              h_hist, f_down_new);
                end
            end

        elseif x_star > L && r_star <= L
            % =====================================================
            %  SUBCASE 2b: INITIATE EXTENDED TRANSFER
            % =====================================================
            m_new        = m_prev;
            k_new        = k_prev;
            x_new        = x_star;
            r_new        = r_star;
            f_down_new   = 0;
            f_sshift_new = 1;
            omega_new    = omega_prev;

            k_down = CA.C{m_prev}(k_prev, 2);

            if f_lshift_prev == 1
                sigma_new = getSigma(svData, tp, params);

                if sigma_new <= 0
                    k_new = k_prev + delta_prev;
                    if k_new < 1 || k_new > CA.K_count(m_prev)
                        k_new = k_prev;
                    end
                    delta_new    = 0;
                    f_lshift_new = 0;
                    sigma_new    = 0;
                    f_stbl       = params.T_stbl;

                    CA = occupyPosition(CA, P_var, tp, m_prev, k_new, r_new, L, i);
                    for j = r_new:L
                        CA.S{m_prev, k_new}(:, j) = [1; double(v_new>0); 0; -1];
                    end
                    CA = occupyPosition(CA, P_var, tp, m_prev+1, k_down, 1, x_new-L, i);
                    for j = 1:(x_new-L)
                        CA.S{m_prev+1, k_down}(:, j) = [1; double(v_new>0); 0; 1];
                    end
                else
                    k_new        = k_prev;
                    delta_new    = delta_prev;
                    f_lshift_new = 1;
                    f_stbl       = 0;

                    CA = occupyPosition(CA, P_var, tp, m_prev, k_prev, r_new, L, i);
                    CA = occupyPosition(CA, P_var, tp, m_prev, k_prev+delta_new, r_new, L, i);
                    for j = r_new:L
                        CA.S{m_prev, k_prev}(:, j)           = [1; double(v_new>0); -1; -1];
                        CA.S{m_prev, k_prev+delta_new}(:, j) = [1; double(v_new>0);  1; -1];
                    end
                    CA = occupyPosition(CA, P_var, tp, m_prev+1, k_down, 1, x_new-L, i);
                    for j = 1:(x_new-L)
                        CA.S{m_prev+1, k_down}(:, j) = [1; double(v_new>0); -1; 1];
                    end
                    if (k_down + delta_new) >= 1 && (k_down + delta_new) <= CA.K_count(m_prev+1)
                        CA = occupyPosition(CA, P_var, tp, m_prev+1, k_down+delta_new, 1, x_new-L, i);
                        for j = 1:(x_new-L)
                            CA.S{m_prev+1, k_down+delta_new}(:, j) = [1; double(v_new>0); 1; 1];
                        end
                    end
                end
            else
                k_new = k_prev;

                CA = occupyPosition(CA, P_var, tp, m_prev, k_new, r_new, L, i);
                for j = r_new:L
                    CA.S{m_prev, k_new}(:, j) = [1; double(v_new>0); 0; -1];
                end
                CA = occupyPosition(CA, P_var, tp, m_prev+1, k_down, 1, x_new-L, i);
                for j = 1:(x_new-L)
                    CA.S{m_prev+1, k_down}(:, j) = [1; double(v_new>0); 0; 1];
                end

                % ---- Lane shift initiation using omega_prev (in source segment) ----
                omega_new    = checkLaneShiftTrigger(svData, tp, CA, params, i, ...
                                              h_hist, f_down_new);
                can_initiate = (omega_new == 2) || (omega_new == 1 && f_stbl_prev == 0);

                if omega_new ~= 0 && can_initiate
                    delta_new = getShiftDirection(CA, params, svData, tp);

                    if delta_new ~= 0
                        feasible = checkLaneShiftFeasibility(svData, tp, CA, params);

                        if feasible
                            k_target  = k_new + delta_new;
                            v_ref     = alpha * v_max / 2;
                            theta_deg = 22.5 * (v_ref / (v_new + v_ref));
                            theta_rad = theta_deg * pi / 180;
                            sigma_new = params.dy / tan(theta_rad);

                            f_lshift_new = 1;
                            f_stbl       = 0;

                            if k_target >= 1 && k_target <= CA.K_count(m_prev)
                                CA = occupyPosition(CA, P_var, tp, m_prev, k_target, r_new, L, i);
                                for j = r_new:L
                                    CA.S{m_prev, k_new}(:, j)    = [1; double(v_new>0); -1; -1];
                                    CA.S{m_prev, k_target}(:, j) = [1; double(v_new>0);  1; -1];
                                end
                            else
                                f_lshift_new = 0; delta_new = 0; sigma_new = 0;
                            end
                        else
                            omega_new = 0; 
                            delta_new = 0; 
                            sigma_new = 0;
                        end
                    else
                        omega_new = 0; 
                        delta_new = 0; 
                        sigma_new = 0;
                    end
                else
                    delta_new = 0;
                    sigma_new = 0;
                    omega_new = checkLaneShiftTrigger(svData, tp, CA, params, i, ...
                                              h_hist, f_down_new);
                end
            end

        else
            % =====================================================
            %  SUBCASE 2c: STAYS WITHIN SEGMENT
            % =====================================================
            m_new        = m_prev;
            k_new        = k_prev;
            x_new        = x_star;
            r_new        = r_star;
            f_down_new   = 1;
            f_sshift_new = 0;
            omega_new    = omega_prev;

            if f_lshift_prev == 1
                sigma_new = getSigma(svData, tp, params);

                if sigma_new <= 0
                    % ---- shift completes ----
                    k_new = k_prev + delta_prev;
                    if k_new < 1 || k_new > CA.K_count(m_prev)
                        k_new = k_prev;
                    end
                    delta_new    = 0;
                    f_lshift_new = 0;
                    sigma_new    = 0;
                    f_stbl       = params.T_stbl;

                    CA = occupyPosition(CA, P_var, tp, m_prev, k_new, r_new, x_new, i);
                    for j = r_new:x_new
                        CA.S{m_prev, k_new}(:, j) = [1; double(v_new>0); 0; 0];
                    end
                else
                    % ---- shift continues ----
                    k_new        = k_prev;
                    delta_new    = delta_prev;
                    f_lshift_new = 1;
                    f_stbl       = 0;
                    k_target     = k_prev + delta_prev;

                    CA = occupyPosition(CA, P_var, tp, m_prev, k_prev, r_new, x_new, i);
                    if k_target >= 1 && k_target <= CA.K_count(m_prev)
                        CA = occupyPosition(CA, P_var, tp, m_prev, k_target, r_new, x_new, i);
                        for j = r_new:x_new
                            CA.S{m_prev, k_prev}(:, j)   = [1; double(v_new>0); -1; 0];
                            CA.S{m_prev, k_target}(:, j) = [1; double(v_new>0);  1; 0];
                        end
                    else
                        for j = r_new:x_new
                            CA.S{m_prev, k_prev}(:, j) = [1; double(v_new>0); 0; 0];
                        end
                        f_lshift_new = 0; delta_new = 0;
                    end
                end

            else
                % ---- Not currently shifting ----
                k_new = k_prev;

                % ---- Straight progression first ----
                CA = occupyPosition(CA, P_var, tp, m_prev, k_prev, r_new, x_new, i);
                for j = r_new:x_new
                    CA.S{m_prev, k_prev}(:, j) = [1; double(v_new>0); 0; 0];
                end

                % ---- Lane shift initiation using omega_prev ----
                omega_new    = omega_prev;
                can_initiate = (omega_new == 2) || (omega_new == 1 && f_stbl_prev == 0);

                if omega_new ~= 0 && can_initiate
                    delta_new = getShiftDirection(CA, params, svData, tp);

                    if delta_new ~= 0
                        feasible = checkLaneShiftFeasibility(svData, tp, CA, params);

                        if feasible
                            k_target  = k_new + delta_new;
                            v_ref     = alpha * v_max / 2;
                            theta_deg = 22.5 * (v_ref / (v_new + v_ref));
                            theta_rad = theta_deg * pi / 180;
                            sigma_new = params.dy / tan(theta_rad);

                            f_lshift_new = 1;
                            f_stbl       = 0;

                            if k_target >= 1 && k_target <= CA.K_count(m_prev)
                                CA = occupyPosition(CA, P_var, tp, m_prev, k_target, r_new, x_new, i);
                                for j = r_new:x_new
                                    CA.S{m_prev, k_prev}(:, j)   = [1; double(v_new>0); -1; 0];
                                    CA.S{m_prev, k_target}(:, j) = [1; double(v_new>0);  1; 0];
                                end
                            else
                                f_lshift_new = 0; delta_new = 0; sigma_new = 0;
                            end
                        else
                            omega_new = 0; 
                            delta_new = 0; 
                            sigma_new = 0;
                        end
                    else
                        omega_new = 0; 
                        delta_new = 0;
                        sigma_new = 0;
                    end
                else
                    delta_new = 0;
                    sigma_new = 0;
                    omega_new = checkLaneShiftTrigger(svData, tp, CA, params, i, ...
                                              h_hist, f_down_new);
                end
            end
        end
    end

    % =========================================================
    %  UPDATE P_var
    % =========================================================
    P_var = updatePvar(i, t_n, a_new, v_new, p_new, h_new, ...
        x_new, r_new, m_new, k_new, omega_new, delta_new, sigma_new, ...
        f_stbl, f_down_new, f_lshift_new, f_sshift_new, P_var);

end

%% getSubjectVehicle 
function [i, svData] = getSubjectVehicle(m_curr, CA, P_var, P_veh, t_n)
% getSubjectVehicle Subject Vehicle Selection
    R   = CA.R;
    L   = size(R{m_curr,1}, 2);
    K_m = CA.K_count(m_curr);

    i      = 0;
    svData = struct();          % always defined, even when no SV is found

    %% NEED TO OPTIMIZE THIS L %%
    for j = L:-1:1
        CV = [];
        n  = 0;
        for k = K_m:-1:1
            if R{m_curr,k}(1, j) == 0
                continue;
            else
                i_temp = R{m_curr,k}(1, j);
            end

            if P_var(i_temp, 9) ~= m_curr
                continue;
            end

            if P_var(i_temp, 18) == t_n
                continue;
            end
            if P_var(i_temp, 10) ~= k
                continue;
            end
            n = n + 1;
            CV(n) = i_temp;
        end

        if n == 0
            continue;
        elseif n == 1
            i      = CV(1);
            svData = svDataBuild(i, P_var, P_veh);
            return;
        else
            i_best = CV(1);
            v_best = P_var(i_best, 4);
            for idx = 2:length(CV)
                i_curr = CV(idx);
                v_curr = P_var(i_curr, 4);
                if v_curr > v_best
                    i_best = i_curr;
                    v_best = v_curr;
                end
            end
            i      = i_best;
            svData = svDataBuild(i, P_var, P_veh);
            return;
        end
    end
end


%% Build svData struct
function svData = svDataBuild(i, P_var, P_veh)
    % if i == 0
    %     fprintf ("SV is not found in step %d", t_n)
    %     return;
    % end
    svData = struct();
    ID_r = P_var(i, 2);
    svData.ID_r = ID_r;
    svData.l_ID = P_veh(ID_r, 3);
    svData.v_max = P_veh(ID_r,4);
    svData.a_max = P_veh(ID_r,5);
    svData.d_max = P_veh(ID_r,6);
    svData.alpha = P_veh(ID_r, 7);
    svData.tau = P_veh(ID_r, 8);
    svData.phi = P_veh(ID_r, 9);
    svData.a_prev = P_var(i, 3);
    svData.v_prev = P_var(i, 4);
    svData.p_prev = P_var(i, 5);
    svData.h_prev = P_var(i, 6);
    svData.x_prev = P_var(i, 7);
    svData.r_prev = P_var(i, 8);
    svData.m_prev = P_var(i, 9);
    svData.k_prev = P_var(i, 10);
    svData.omega_prev = P_var(i, 11);
    svData.delta_prev = P_var(i, 12);
    svData.sigma_prev = P_var(i, 13);
    svData.f_stbl_prev = P_var(i, 14);
    svData.f_down_prev = P_var(i, 15);
    svData.f_lshift_prev = P_var(i, 16);
    svData.f_sshift_prev = P_var(i, 17);
end


%% update P_var
function P_var = updatePvar(i, t_n, ...
    a_new, v_new, p_new, h_new, ...
    x_new, r_new, m_new, k_new, ...
    omega_new, delta_new, sigma_new, ...
    f_stbl, f_down_new, f_lshift_new, f_sshift_new, P_var)

    P_var(i,3)  = a_new;
    P_var(i,4)  = v_new;
    P_var(i,5)  = p_new;
    P_var(i,6)  = h_new;
    P_var(i,7)  = x_new;
    P_var(i,8)  = r_new;
    P_var(i,9)  = m_new;
    P_var(i,10) = k_new;
    P_var(i,11) = omega_new;
    P_var(i,12) = delta_new;
    P_var(i,13) = sigma_new;
    P_var(i,14) = f_stbl-1;
    P_var(i,15) = f_down_new;
    P_var(i,16) = f_lshift_new;
    P_var(i,17) = f_sshift_new;
    P_var(i,18) = t_n;
end


%%  getTargetPosition
function [tp, h_hist] = getTargetPosition(CA, P_var, params, h_hist, svData, i)

    % =========================================================
    %  STEP 1 — Extract state
    % =========================================================

    % ---- Static properties ----
    l_ID  = svData.l_ID;
    v_max = svData.v_max;
    a_max = svData.a_max;
    d_max = svData.d_max;
    alpha = svData.alpha;
    tau   = svData.tau;

    % ---- Dynamic state (source frame) ----
    a_prev        = svData.a_prev;
    v_prev        = svData.v_prev;
    x_prev        = svData.x_prev;
    r_prev        = svData.r_prev;
    m_prev        = svData.m_prev;
    k_prev        = svData.k_prev;
    delta_prev    = svData.delta_prev;
    f_down_prev   = svData.f_down_prev;
    f_lshift_prev = svData.f_lshift_prev;
    f_sshift_prev = svData.f_sshift_prev;

    % ---- Parameters ----
    dt     = params.dt;
    g_safe = params.g_safe;
    p_min  = params.p_min;

    M     = size(CA.R, 1);
    K_cnt = CA.K_count;

    % =========================================================
    %  STEP 2 — Effective frame
    % =========================================================
    eff   = getEffectivePosition(CA, svData);
    x_eff = eff.x_eff;
    m_eff = eff.m_eff;
    k_eff = eff.k_eff;

    % =========================================================
    %  STEP 3 — Heading gap and perception range
    % =========================================================
    [h_new, p_new] = getHeadingGap(CA, params, x_eff, m_eff, k_eff, svData);

    if f_lshift_prev == 1
        k_tgt = k_eff + delta_prev;
        if k_tgt >= 1 && k_tgt <= K_cnt(m_eff)
            [h_tgt, p_tgt] = getHeadingGap(CA, params, x_eff, m_eff, k_tgt, svData);
            h_new = min(h_new, h_tgt);
            p_new = min(p_new, p_tgt);
        end
    end

    % =========================================================
    %  STEP 4 — Update gap history; compute Δh
    % =========================================================
    h_prev2      = h_hist(i, 2);
    h_hist(i, 2) = h_hist(i, 1);
    h_hist(i, 1) = h_new;

    dh = (h_new - h_prev2) / (2 * dt);

    % =========================================================
    %  STEP 5 — Braking tolerance
    % =========================================================
    if dh <= 0
        brk_tol = -min(alpha * v_prev, 0.4 * v_prev);
    else
        brk_tol = 0;
    end

    % =========================================================
    %  STEP 6 — Acceleration decision (basic)
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

    % =========================================================
    %  STEP 7 — Cooperative deceleration override
    % =========================================================
    if alpha < 0.5

        L_prev = size(CA.R{m_prev, k_prev}, 2);

        % =====================================================
        %  LEFT side (k_prev − 1)
        % =====================================================
        k_left = k_prev - 1;
        if k_left >= 1

            obs_left          = false;
            obs_left_distance = inf;

            % ---- 7L.1: forbidden scan in current segment ----
            for j = x_prev+1 : min(L_prev, x_prev + p_new)              % LOOP: forbidden scan (left, current)
                if CA.F{m_prev, k_left}(j) == 1
                    obs_left          = true;
                    obs_left_distance = j - x_prev;
                    break;                                                   % break: forbidden scan (left, current)
                end
            end

            % ---- 7L.2: lane drop at downstream edge ----
            if f_down_prev == 1 && (x_prev + 4*v_prev*dt) > L_prev && ...
                    k_left <= size(CA.C{m_prev}, 1)

                mu_left = CA.C{m_prev}(k_left, 3);

                if mu_left == -1
                    obs_left          = true;
                    obs_left_distance = min(obs_left_distance, L_prev - x_prev);

                elseif (m_prev + 1) <= M
                    k_down_left = CA.C{m_prev}(k_left, 2);
                    if k_down_left > 0 && k_down_left <= K_cnt(m_prev + 1)
                        L_down_left = size(CA.R{m_prev + 1, k_down_left}, 2);
                        for j = 1 : min(x_prev + p_new - L_prev, L_down_left)     % LOOP: forbidden scan (left, downstream)
                            if CA.F{m_prev + 1, k_down_left}(j) == 1
                                obs_left          = true;
                                obs_left_distance = min(obs_left_distance, L_prev + j - x_prev);
                                break;                                                % break: forbidden scan (left, downstream)
                            end
                        end
                    end
                end
            end

            % ---- 7L.3: search for merging vehicle ----
            if obs_left

                veh_found = false;

                if f_sshift_prev == 0

                    % Current segment of k_left
                    j0 = max(1, x_prev - p_min);
                    for j = j0 : x_prev                                       % LOOP: merging-vehicle scan (left, current)
                        if CA.R{m_prev, k_left}(1, j) ~= 0
                            i_left = CA.R{m_prev, k_left}(1, j);
                            v_left = P_var(i_left, 4);
                            x_left = P_var(i_left, 7);
                            if x_left >= x_prev - g_safe && v_left >= v_prev - 1
                                veh_found = true;
                                break;                                            % break: merging-vehicle scan (left, current)
                            end
                        end
                    end

                    % Downstream segment of k_left
                    if ~veh_found && (x_prev + p_min) > L_prev && ...
                            k_left <= size(CA.C{m_prev}, 1) && ...
                            CA.C{m_prev}(k_left, 3) ~= -1

                        k_down_left = CA.C{m_prev}(k_left, 2);
                        if k_down_left > 0 && (m_prev + 1) <= M
                            L_down_left = size(CA.R{m_prev + 1, k_down_left}, 2);
                            for j = 1 : min(x_prev + p_min - L_prev, L_down_left)   % LOOP: merging-vehicle scan (left, downstream)
                                if CA.R{m_prev + 1, k_down_left}(1, j) ~= 0
                                    i_left = CA.R{m_prev + 1, k_down_left}(1, j);
                                    v_left = P_var(i_left, 4);
                                    if v_left >= v_prev - 1
                                        veh_found = true;
                                        break;                                            % break: merging-vehicle scan (left, downstream)
                                    end
                                end
                            end
                        end
                    end

                else
                    % Extended transfer: only downstream segment matters
                    if k_left <= size(CA.C{m_prev}, 1)
                        mu_left = CA.C{m_prev}(k_left, 3);
                        if mu_left ~= -1 && (m_prev + 1) <= M
                            k_down_left = CA.C{m_prev}(k_left, 2);
                            if k_down_left > 0 && k_down_left <= K_cnt(m_prev + 1)
                                L_down_left = size(CA.R{m_prev + 1, k_down_left}, 2);
                                for j = 1 : min(x_eff + p_min, L_down_left)        % LOOP: merging-vehicle scan (left, downstream only)
                                    if CA.R{m_prev + 1, k_down_left}(1, j) ~= 0
                                        i_left = CA.R{m_prev + 1, k_down_left}(1, j);
                                        v_left = P_var(i_left, 4);
                                        if v_left >= v_prev - 1
                                            veh_found = true;
                                            break;                                    % break: merging-vehicle scan (left, downstream only)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                if veh_found && obs_left_distance > (x_prev + 4*v_prev*dt)
                    a_dec = -2;
                end
            end
        end

        % =====================================================
        %  RIGHT side (k_prev + 1) — symmetric to LEFT
        % =====================================================
        k_right = k_prev + 1;
        if k_right <= K_cnt(m_prev)

            obs_right          = false;
            obs_right_distance = inf;

            % ---- 7R.1: forbidden scan in current segment ----
            for j = x_prev+1 : min(L_prev, x_prev + p_new)              % LOOP: forbidden scan (right, current)
                if CA.F{m_prev, k_right}(j) == 1
                    obs_right          = true;
                    obs_right_distance = j - x_prev;
                    break;                                                   % break: forbidden scan (right, current)
                end
            end

            % ---- 7R.2: lane drop at downstream edge ----
            if f_down_prev == 1 && (x_prev + 4*v_prev*dt) > L_prev && ...
                    k_right <= size(CA.C{m_prev}, 1)

                mu_right = CA.C{m_prev}(k_right, 3);

                if mu_right == -1
                    obs_right          = true;
                    obs_right_distance = min(obs_right_distance, L_prev - x_prev);

                elseif (m_prev + 1) <= M
                    k_down_right = CA.C{m_prev}(k_right, 2);
                    if k_down_right > 0 && k_down_right <= K_cnt(m_prev + 1)
                        L_down_right = size(CA.R{m_prev + 1, k_down_right}, 2);
                        for j = 1 : min(x_prev + p_new - L_prev, L_down_right)    % LOOP: forbidden scan (right, downstream)
                            if CA.F{m_prev + 1, k_down_right}(j) == 1
                                obs_right          = true;
                                obs_right_distance = min(obs_right_distance, L_prev + j - x_prev);
                                break;                                                % break: forbidden scan (right, downstream)
                            end
                        end
                    end
                end
            end

            % ---- 7R.3: search for merging vehicle ----
            if obs_right

                veh_found = false;

                if f_sshift_prev == 0

                    j0 = max(1, x_prev - p_min);
                    for j = j0 : x_prev                                       % LOOP: merging-vehicle scan (right, current)
                        if CA.R{m_prev, k_right}(1, j) ~= 0
                            i_right = CA.R{m_prev, k_right}(1, j);
                            v_right = P_var(i_right, 4);
                            x_right = P_var(i_right, 7);
                            if x_right >= x_prev - g_safe && v_right >= v_prev - 1
                                veh_found = true;
                                break;                                            % break: merging-vehicle scan (right, current)
                            end
                        end
                    end

                    if ~veh_found && (x_prev + p_min) > L_prev && ...
                            k_right <= size(CA.C{m_prev}, 1) && ...
                            CA.C{m_prev}(k_right, 3) ~= -1

                        k_down_right = CA.C{m_prev}(k_right, 2);
                        if k_down_right > 0 && (m_prev + 1) <= M
                            L_down_right = size(CA.R{m_prev + 1, k_down_right}, 2);
                            for j = 1 : min(x_prev + p_min - L_prev, L_down_right)   % LOOP: merging-vehicle scan (right, downstream)
                                if CA.R{m_prev + 1, k_down_right}(1, j) ~= 0
                                    i_right = CA.R{m_prev + 1, k_down_right}(1, j);
                                    v_right = P_var(i_right, 4);
                                    if v_right >= v_prev - 1
                                        veh_found = true;
                                        break;                                            % break: merging-vehicle scan (right, downstream)
                                    end
                                end
                            end
                        end
                    end

                else
                    if k_right <= size(CA.C{m_prev}, 1)
                        mu_right = CA.C{m_prev}(k_right, 3);
                        if mu_right ~= -1 && (m_prev + 1) <= M
                            k_down_right = CA.C{m_prev}(k_right, 2);
                            if k_down_right > 0 && k_down_right <= K_cnt(m_prev + 1)
                                L_down_right = size(CA.R{m_prev + 1, k_down_right}, 2);
                                for j = 1 : min(x_eff + p_min, L_down_right)        % LOOP: merging-vehicle scan (right, downstream only)
                                    if CA.R{m_prev + 1, k_down_right}(1, j) ~= 0
                                        i_right = CA.R{m_prev + 1, k_down_right}(1, j);
                                        v_right = P_var(i_right, 4);
                                        if v_right >= v_prev - 1
                                            veh_found = true;
                                            break;                                    % break: merging-vehicle scan (right, downstream only)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                if veh_found && obs_right_distance > (x_prev + 4*v_prev*dt)
                    a_dec = -2;
                end
            end
        end
    end

    % =========================================================
    %  STEP 8 — Applied acceleration
    % =========================================================
    if a_dec == 1
        a_star = a_prev + alpha * (a_max - max(0, a_prev));
    elseif a_dec == 0
        a_star = a_prev;
    elseif a_dec == -1
        u_urg  = 1 - brk_tol / max(dh, -1);
        a_star = -min(d_max, ceil(alpha * d_max * u_urg));
    elseif a_dec == -2
        a_star = a_prev - ceil((1 - alpha) * abs(a_prev));
    end

    a_star = max(-d_max, min(a_max, a_star));

    % =========================================================
    %  STEP 9 — Feasible velocity
    % =========================================================
    v_star = max(0, min(v_prev + a_star * dt, v_max));

    L_src = CA.L_seg(m_prev);

    % ---- Source segment scan ----
    if f_sshift_prev == 0
        x_reach = round(x_prev + v_star * dt);
        j_src   = r_prev : min(x_reach, L_src);
    else
        j_src = r_prev : L_src;
    end

    M_spd_min = min(CA.M_spd{m_prev, k_prev}(j_src));

    if delta_prev ~= 0
        k_tgt_src = k_prev + delta_prev;
        if k_tgt_src >= 1 && k_tgt_src <= K_cnt(m_prev)
            if ~any(CA.F{m_prev, k_tgt_src}(j_src) == 1)
                cand      = min(CA.M_spd{m_prev, k_tgt_src}(j_src));
                M_spd_min = min(M_spd_min, cand);
            end
        end
    end

    % ---- Downstream segment scan (recompute forward reach after source scan) ----
    if f_sshift_prev == 0
        x_reach = round(x_prev + v_star * dt);
    else
        x_reach = round(x_eff + v_star * dt);     % in downstream frame
    end

    if x_reach > L_src && m_prev < M

        mu = CA.C{m_prev}(k_prev, 3);

        if mu ~= -1

            m_down = m_prev + 1;
            k_down = CA.C{m_prev}(k_prev, 2);

            if k_down > 0 && k_down <= K_cnt(m_down)

                L_down = CA.L_seg(m_down);

                if f_sshift_prev == 1
                    j_down = 1 : min(x_reach, L_down);
                else
                    j_down = 1 : min(x_reach - L_src, L_down);
                end

                if ~isempty(j_down)

                    M_spd_min = min(M_spd_min, min(CA.M_spd{m_down, k_down}(j_down)));

                    if delta_prev ~= 0
                        k_tgt_down = k_down + delta_prev;
                        if k_tgt_down >= 1 && k_tgt_down <= K_cnt(m_down)
                            if ~any(CA.F{m_down, k_tgt_down}(j_down) == 1)
                                cand      = min(CA.M_spd{m_down, k_tgt_down}(j_down));
                                M_spd_min = min(M_spd_min, cand);
                            end
                        end
                    end
                end
            end
        end
    end

    if v_star > M_spd_min
        v_star = M_spd_min;
    end

    v_new = round(v_star);
    a_new = round((v_new - v_prev) / dt);

    % =========================================================
    %  STEP 10 — Target positions (source frame)
    % =========================================================
    if f_lshift_prev == 1
        v_ref     = alpha * v_max / 2;
        theta_deg = 22.5 * (v_ref / (v_new + v_ref));
        theta_rad = theta_deg * pi / 180;
        d_long    = v_new * dt * cos(theta_rad);

        x_star = round(x_prev + d_long);
        r_star = x_star - round(l_ID * cos(theta_rad)) + 1;
    else
        x_star = round(x_prev + v_new * dt);
        r_star = x_star - l_ID + 1;
    end

    % =========================================================
    %  STEP 11 — Segment-shift flag
    % =========================================================
    f_sshift_star = 0;

    if k_prev <= size(CA.C{m_prev}, 1)
        mu = CA.C{m_prev}(k_prev, 3);
    else
        mu = -1;
    end

    if x_star > L_src && r_star <= L_src && mu ~= -1
        f_sshift_star = 1;
    end

    % =========================================================
    %  STEP 12 — Assemble output
    % =========================================================
    tp = struct();
    tp.p_new         = round(p_new);
    tp.h_new         = round(h_new);
    tp.a_dec         = a_dec;
    tp.a_new         = a_new;
    tp.v_new         = v_new;
    tp.x_star        = x_star;
    tp.r_star        = r_star;
    tp.f_sshift_star = f_sshift_star;

end


%%  getEffectivePosition
function eff_pos = getEffectivePosition(CA, svData)

    % ---- Extract source-frame state ----
    m_prev        = svData.m_prev;
    k_prev        = svData.k_prev;
    x_prev        = svData.x_prev;
    f_sshift_prev = svData.f_sshift_prev;

    % ---- Default: source frame is the effective frame ----
    m_eff = m_prev;
    k_eff = k_prev;
    x_eff = x_prev;

    % ---- Override when the front-end has crossed the segment boundary ----
    if f_sshift_prev == 1

        L_prev = size(CA.R{m_prev, k_prev}, 2);

        m_eff = m_prev + 1;
        k_eff = CA.C{m_prev}(k_prev, 2);
        x_eff = x_prev - L_prev;

        % r_eff stays in the source frame — do not transform it.
    end

    eff_pos = struct();
    eff_pos.x_eff = x_eff;
    eff_pos.r_eff = svData.r_prev;
    eff_pos.m_eff = m_eff;
    eff_pos.k_eff = k_eff;
end


%% occPos
function CA = occupyPosition(CA, P_var, tp, m_curr, k_curr, j_start, j_end, i)
    
    for j = j_start:j_end
        CA.R{m_curr, k_curr}(1, j)  = i;
        CA.R{m_curr, k_curr}(2, j)  = P_var(i, 2);
        CA.V{m_curr, k_curr}(j)     = tp.v_new;
        CA.A{m_curr, k_curr}(j)     = tp.a_new;
        CA.G{m_curr, k_curr}(j)     = 0;
    end
    
end


%% clearOcc
function CA = clearOccupancy(CA, m_curr, k_curr, j_start, j_end)

    L  = CA.L_seg(m_curr);
    j_start = max(1, j_start);
    j_end = min(L, j_end);

    for j = j_start:j_end
        CA.R{m_curr, k_curr}(:, j)  = 0;
        CA.V{m_curr, k_curr}(j)     = 0;
        CA.A{m_curr, k_curr}(j)     = 0;
        CA.G{m_curr, k_curr}(j)     = 1;
        CA.S{m_curr, k_curr}(:, j)  = 0;
    end
end


%%  getShiftDirection
function delta_new = getShiftDirection(CA, params, svData, tp)

    % =========================================================
    %  STEP 1 — Extract state
    % =========================================================

    % ---- From svData ----
    m_prev   = svData.m_prev;
    k_prev   = svData.k_prev;
    tau      = svData.tau;
    alpha    = svData.alpha;
    d_max    = svData.d_max;

    % ---- From tp (post-step state) ----
    v_new         = tp.v_new;
    x_star        = tp.x_star;
    f_sshift_star = tp.f_sshift_star;

    % ---- From params ----
    p_max  = params.p_max;
    p_min  = params.p_min;
    g_safe = params.g_safe;

    M     = size(CA.R, 1);
    K_cnt = CA.K_count;

    delta_new = svData.delta_prev;

    % =========================================================
    %  STEP 2 — Determine effective (segment, lane, x)
    % =========================================================

    L_prev = CA.L_seg(m_prev);

    if f_sshift_star == 1
        m_sv = m_prev + 1;
        k_sv = CA.C{m_prev}(k_prev, 2);
        x_sv = x_star - L_prev;      % downstream-frame position
    else
        m_sv = m_prev;
        k_sv = k_prev;
        x_sv = x_star;
    end

    if k_sv == 0
        delta_new = 0;
        return;
    end

    % =========================================================
    %  STEP 3 — Boundary lane check
    % =========================================================

    k_lanes = K_cnt(m_sv);

    if k_sv == 1
        delta_new = 1;
        return;
    elseif k_sv == k_lanes
        delta_new = -1;
        return;
    end

    k_left  = k_sv - 1;
    k_right = k_sv + 1;

    % =========================================================
    %  STEP 4 — Continuation check (only when f_down_new == 1)
    % =========================================================

    f_down_new = double(x_sv >= CA.L_seg(m_sv) - params.l_down);

    if f_down_new == 1

        if k_left <= size(CA.C{m_sv}, 1)
            mu_left = CA.C{m_sv}(k_left, 3);
        else
            mu_left = 0;
        end

        if k_right <= size(CA.C{m_sv}, 1)
            mu_right = CA.C{m_sv}(k_right, 3);
        else
            mu_right = 0;
        end

        if mu_left == -1 && mu_right ~= -1
            delta_new = 1;
            return;
        elseif mu_left ~= -1 && mu_right == -1
            delta_new = -1;
            return;
        elseif mu_left == -1 && mu_right == -1
            delta_new = 0;
            return;
        end
    end

    % =========================================================
    %  STEP 5 — Perception range (at the new position)
    % =========================================================
    p_star   = ceil(v_new * tau + (v_new^2) / (2 * alpha * d_max)) + g_safe;
    p_latest = min(p_max, max(p_min, p_star));

    % =========================================================
    %  STEP 6 — Obstruction scan (left and right)
    % =========================================================

    L_eff = CA.L_seg(m_sv);
    j_end = min(x_sv + p_latest, L_eff);

    % -----------------------------------------------------
    %  Left lane: current segment portion
    % -----------------------------------------------------
    left_blocked = false;
    d_obs_left   = p_latest;

    for j = x_sv+1 : j_end                        % LOOP: obstruction scan (left, current)
        if CA.F{m_sv, k_left}(j) == 1
            left_blocked = true;
            d_obs_left   = j - x_sv;
            break;                                 % break: obstruction scan (left, current)
        end
    end

    % -----------------------------------------------------
    %  Left lane: downstream portion
    % -----------------------------------------------------
    if ~left_blocked && (x_sv + p_latest) > L_eff && m_sv < M

        if k_left <= size(CA.C{m_sv}, 1)
            mu_left = CA.C{m_sv}(k_left, 3);
        else
            mu_left = -1;
        end

        if mu_left ~= -1
            k_down_left = CA.C{m_sv}(k_left, 2);
            if k_down_left > 0 && k_down_left <= K_cnt(m_sv + 1)
                L_next = CA.L_seg(m_sv + 1);
                for j = 1 : min(x_sv + p_latest - L_eff, L_next)      % LOOP: obstruction scan (left, downstream)
                    if CA.F{m_sv+1, k_down_left}(j) == 1
                        left_blocked = true;
                        d_obs_left   = (L_eff - x_sv) + j;
                        break;                                    % break: obstruction scan (left, downstream)
                    end
                end
            end
        end
    end

    % -----------------------------------------------------
    %  Right lane: current segment portion
    % -----------------------------------------------------
    right_blocked = false;
    d_obs_right   = p_latest;

    for j = x_sv+1 : j_end                        % LOOP: obstruction scan (right, current)
        if CA.F{m_sv, k_right}(j) == 1
            right_blocked = true;
            d_obs_right   = j - x_sv;
            break;                                 % break: obstruction scan (right, current)
        end
    end

    % -----------------------------------------------------
    %  Right lane: downstream portion
    % -----------------------------------------------------
    if ~right_blocked && (x_sv + p_latest) > L_eff && m_sv < M

        if k_right <= size(CA.C{m_sv}, 1)
            mu_right = CA.C{m_sv}(k_right, 3);
        else
            mu_right = -1;
        end

        if mu_right ~= -1
            k_down_right = CA.C{m_sv}(k_right, 2);
            if k_down_right > 0 && k_down_right <= K_cnt(m_sv + 1)
                L_next = CA.L_seg(m_sv + 1);
                for j = 1 : min(x_sv + p_latest - L_eff, L_next)      % LOOP: obstruction scan (right, downstream)
                    if CA.F{m_sv+1, k_down_right}(j) == 1
                        right_blocked = true;
                        d_obs_right   = (L_eff - x_sv) + j;
                        break;                                    % break: obstruction scan (right, downstream)
                    end
                end
            end
        end
    end

    % =========================================================
    %  STEP 7 — Decide based on obstruction status
    % =========================================================

    if left_blocked && ~right_blocked
        delta_new = 1;
        return;
    elseif ~left_blocked && right_blocked
        delta_new = -1;
        return;
    elseif left_blocked && right_blocked
        if d_obs_right > d_obs_left
            delta_new = 1;
        elseif d_obs_left > d_obs_right
            delta_new = -1;
        else
            delta_new = 0;
        end
        return;
    end

    % =========================================================
    %  STEP 8 — Both lanes free → compare heading gaps
    % =========================================================
    [h_left,  ~] = getHeadingGap(CA, params, x_sv, m_sv, k_left,  svData);
    [h_right, ~] = getHeadingGap(CA, params, x_sv, m_sv, k_right, svData);

    if h_right >= h_left
        delta_new = 1;
    else
        delta_new = -1;
    end

    % existency check
    k_check = k_sv + delta_new;
    if k_check < 1 || k_check > CA.K_count(m_prev)
        delta_new = 0;
    end

end


%%  checkLaneShiftTrigger
function omega_new = checkLaneShiftTrigger(svData, tp, CA, params, i, h_hist, f_down_new)

    % -------- Extract needed values --------
    v_new         = tp.v_new;
    h_new         = tp.h_new;
    f_sshift_star = tp.f_sshift_star;

    l_ID          = svData.l_ID;
    alpha         = svData.alpha;
    tau           = svData.tau;
    d_max         = svData.d_max;
    m_prev        = svData.m_prev;
    k_prev        = svData.k_prev;
    omega_prev    = svData.omega_prev;
    omega_new     = omega_prev;

    p_max  = params.p_max;
    p_min  = params.p_min;
    g_safe = params.g_safe;
    dt     = params.dt;

    L       = CA.L_seg(m_prev);
    mu      = CA.C{m_prev}(k_prev, 3);

    % -------- Determine effective (segment, lane, x) --------
    %  The effective frame identifies where the front-end currently is
    %  (or will be) at the next step, matching getEffectivePosition's
    %  convention.
    if f_sshift_star ~= 1
        m_star  = m_prev;
        k_star  = k_prev;
        x_new   = tp.x_star;
        k_left  = k_prev - 1;
        k_right = k_prev + 1;
    else
        m_star  = m_prev + 1;
        k_star  = CA.C{m_prev}(k_prev, 2);
        x_new   = tp.x_star - L;               % downstream-frame position
        k_left  = k_star - 1;
        k_right = k_star + 1;
    end

    k_lanes = CA.K_count(m_star);


    % =========================================================
    %  CHECK 1 — HEADING GAP COMPRESSION
    % =========================================================
    k_left_exist  = (k_left  >= 1);
    k_right_exist = (k_right <= k_lanes);
    dh = (h_new - h_hist(i, 2)) / (2 * dt);

    if dh < 0

        h_left  = -Inf;
        h_right = -Inf;

        if k_right_exist
            [h_right, ~] = getHeadingGap(CA, params, x_new, m_star, k_right, svData);
        end

        if k_left_exist
            [h_left, ~] = getHeadingGap(CA, params, x_new, m_star, k_left, svData);
        end

        if h_left >= h_new || h_right >= h_new
            omega_new = max(omega_prev, 1);
        end
    end

    % =========================================================
    %  CHECK 2 — LANE TERMINATION AT DOWNSTREAM EDGE
    % =========================================================
    if f_down_new == 1 && mu == -1 && f_sshift_star == 0
        if L - x_new - g_safe > 4 * v_new
            omega_new = 1;
        else
            omega_new = 2;
        end
    end

    % =========================================================
    %  CHECK 3 — STATIC OBSTRUCTION AHEAD
    % =========================================================

    % -------- Effective segment and lane --------
    if f_sshift_star == 1
        m_eff = m_star;                        % = m_prev + 1
        k_eff = k_star;                        % = C{m_prev}(k_prev, 2)
    else
        m_eff = m_prev;
        k_eff = k_prev;
    end

    % -------- Effective segment length --------
    L_eff = CA.L_seg(m_eff);

    % -------- Perception range for obstruction scan --------
    if f_down_new == 1 && mu == -1 && f_sshift_star ~= 1
        p_obs = L - x_new;
    else
        p_star = ceil(v_new * tau + (v_new^2) / (2 * alpha * d_max)) + g_safe;
        p_obs  = min(p_max, max(p_min, p_star));
    end


    % Set scan range --------
    j_end = x_new + p_obs;
    d_obs = 0;

    if j_end <= L_eff

        for j = x_new+1 : j_end                    % LOOP: obstruction scan (current segment)
            if CA.F{m_eff, k_eff}(j) == 1
                break;                             % break: obstruction scan (current segment)
            end
            d_obs = d_obs + 1;
        end

    else

        % -------- current segment portion --------
        d_curr = 0;
        for j = x_new+1 : L_eff                    % LOOP: obstruction scan (current portion)
            if CA.F{m_eff, k_eff}(j) == 1
                break;                             % break: obstruction scan (current portion)
            end
            d_curr = d_curr + 1;
        end

        % -------- downstream segment portion --------
        d_down = 0;

        if m_eff < size(CA.R, 1)                   % not the last segment
            k_down = CA.C{m_eff}(k_eff, 2);
            if k_down > 0 && (m_eff + 1) <= size(CA.F, 1) && ...
                    k_down <= size(CA.F{m_eff+1}, 2)
                L_next = size(CA.F{m_eff+1, k_down}, 2);
                for j = 1 : min(j_end - L_eff, L_next)    % LOOP: obstruction scan (downstream portion)
                    if CA.F{m_eff+1, k_down}(j) == 1
                        break;                             % break: obstruction scan (downstream portion)
                    end
                    d_down = d_down + 1;
                end
            end
        else
            % Last physical segment — beyond is exit buffer (empty)
            d_down = min(j_end - L_eff, p_max);
        end

        d_obs = d_curr + d_down;
    end

    % =========================================================
    %  SETTING OMEGA
    % =========================================================
    if d_obs < (j_end - x_new)
        if d_obs > 4 * v_new - l_ID
            omega_new = max(omega_new, 1);
        else
            omega_new = 2;
        end
    end

end


%%  checkLaneShiftFeasibility
function feasible = checkLaneShiftFeasibility(CA, params, svData, tp)

    feasible = false;

    % ---- Extract inputs ----
    m_prev        = svData.m_prev;
    k_prev        = svData.k_prev;
    r_prev        = svData.r_prev;
    delta_prev    = svData.delta_prev;
    f_sshift_prev = svData.f_sshift_prev;
    f_sshift_star = tp.f_sshift_star;
    x_star        = ceil(tp.x_star);
    g_safe        = params.g_safe;

    L_src = CA.L_seg(m_prev);
    M     = size(CA.R, 1);

    % ---- Target lane in source frame ----
    k_target_src = k_prev + delta_prev;
    if k_target_src < 1 || k_target_src > CA.K_count(m_prev)
        return;
    end

    % =========================================================
    %  CASE A — Vehicle is or will be in extended transfer
    % =========================================================
    if (f_sshift_prev == 1 || f_sshift_star == 1) && m_prev < M

        m_down = m_prev + 1;

        % ---- Continuation check for current lane ----
        mu     = CA.C{m_prev}(k_prev, 3);
        if mu == -1
            return;         % current lane drops; no shift target exists
        end
        k_down = CA.C{m_prev}(k_prev, 2);
        if k_down == 0
            return;         % defensive: mu inconsistent with mapping
        end

        k_target_down = k_down + delta_prev;
        if k_target_down < 1 || k_target_down > CA.K_count(m_down)
            return;
        end

        % ---- A1: source-segment part (rear half) ----
        j_start = max(1, r_prev - g_safe);
        for j = j_start : L_src             % LOOP: source occupancy
            if CA.G{m_prev, k_target_src}(j) == 0
                return;
            end
        end

        % ---- A2: downstream-segment part (front half) ----
        x_eff = x_star - L_src;
        if x_eff < 1
            % Invariant violation: extended transfer implies front past L_src.
            return;
        end

        L_down     = size(CA.R{m_down, k_target_down}, 2);
        j_end_down = min(x_eff, L_down);

        for j = 1 : j_end_down               % LOOP: downstream occupancy
            if CA.G{m_down, k_target_down}(j) == 0
                return;
            end
        end

        feasible = true;
        return;
    end

    % =========================================================
    %  CASE B — Vehicle stays within the current segment
    % =========================================================
    if x_star <= L_src

        j_start = max(1, r_prev - g_safe);
        for j = j_start : x_star             % LOOP: single-segment occupancy
            if CA.G{m_prev, k_target_src}(j) == 0
                return;
            end
        end

        feasible = true;
        return;
    end

    % =========================================================
    %  CASE C — Vehicle crosses the physical boundary at m_prev == M
    % =========================================================
    if m_prev == M && x_star > L_src
        % Last physical segment: boundary is the exit; nothing downstream.
        j_start = max(1, r_prev - g_safe);
        j_end   = min(x_star, L_src);
        for j = j_start : j_end             % LOOP: single-segment occupancy
            if CA.G{m_prev, k_target_src}(j) == 0
                return;
            end
        end

        feasible = true;
        return;
    end

    % =========================================================
    %  Defensive fallback (unreachable if invariants hold)
    % =========================================================
    warning('checkLaneShiftFeasibility:unhandledState', ...
        'm_prev=%d f_sshift_prev=%d f_sshift_star=%d x_star=%d L_src=%d', ...
        m_prev, f_sshift_prev, f_sshift_star, x_star, L_src);

end


%%  getHeadingGap
function [h_new, p_new] = getHeadingGap(CA, params, x_sv, m_sv, k_sv, svData)

    % =========================================================
    %  COMPUTE PERCEPTION RANGE
    % =========================================================
    p_max  = params.p_max;
    p_min  = params.p_min;
    g_safe = params.g_safe;

    v_prev = svData.v_prev;
    tau    = svData.tau;
    alpha  = svData.alpha;
    d_max  = svData.d_max;

    p_star = ceil( v_prev * tau + (v_prev^2) / (2 * alpha * d_max) ) + g_safe;
    p_new  = min( p_max , max( p_min , p_star ) );

    L_sv = CA.L_seg(m_sv); % size(CA.R{m_sv, k_sv}, 2);

    % ---- Continuation flag for this lane ----
    mu_sv = CA.C{m_sv}(k_sv, 3);

    % ---- Terminating-lane cap (only applies to the subject vehicle) ----
    if svData.f_down_prev == 1 && svData.f_sshift_prev == 0 && mu_sv == -1
        p_new = min( p_new , max(0, L_sv - x_sv) );
    end

    % =========================================================
    %  CURRENT-SEGMENT SCAN
    % =========================================================
    h_curr   = 0;
    scan_end = min( x_sv + p_new , L_sv );

    for j = x_sv + 1 : scan_end              % LOOP: forward scan (current segment)
        if CA.G{m_sv, k_sv}(j) == 0
            break;                           % break: forward scan — occupied/forbidden
        end
        h_curr = h_curr + 1;
    end

    % =========================================================
    %  DOWNSTREAM SCAN (only if perception overflows)
    % =========================================================
    h_down   = 0;
    overflow = (x_sv + p_new) - L_sv;

    if overflow > 0 && mu_sv ~= -1

        if m_sv == size(CA.R, 1)
            % Last physical segment: exit buffer is empty
            h_down = min(overflow, p_max);

        else
            k_down = CA.C{m_sv}(k_sv, 2);

            if k_down > 0 && ...
               (m_sv + 1) <= size(CA.R, 1) && ...
               k_down <= size(CA.G{m_sv + 1}, 2)

                L_down        = size(CA.G{m_sv + 1, k_down}, 2);
                scan_end_down = min(overflow, L_down);

                for j = 1 : scan_end_down    % LOOP: forward scan (downstream segment)
                    if CA.G{m_sv + 1, k_down}(j) == 0
                        break;               % break: forward scan — occupied/forbidden
                    end
                    h_down = h_down + 1;
                end
            end
        end
    end

    % =========================================================
    %  CAP AND RETURN
    % =========================================================
    h_new = min(p_max, h_curr + h_down);

end


%% getSigma
function sigma_new = getSigma(svData, tp, params)

    v_new      = tp.v_new;
    sigma_prev = svData.sigma_prev;
    v_max      = svData.v_max;
    alpha      = svData.alpha;
    dt         = params.dt;

    if v_new > 0
        % Reference speed (scales with driver aggressiveness and capability)
        v_ref = alpha * v_max / 2;

        % Shifting angle (max 22.5°, decreases with speed)
        theta_deg = 22.5 * (v_ref / (v_new + v_ref));
        theta_rad = theta_deg * pi / 180;

        % Lateral progress this step
        d_lat = v_new * dt * sin(theta_rad);

        % Decrement remaining distance (floored at 0)
        sigma_new = max(0, sigma_prev - d_lat);

    else
        % Stationary vehicle makes no lateral progress
        sigma_new = sigma_prev;
    end
end