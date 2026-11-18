% =========================================================
%  STAGE 1: STATIC SETUP
% =========================================================

%% -------------------------------------------------------
%  1.1 GLOBAL CONTROL PARAMETERS
%% -------------------------------------------------------

dt      = 1;       % time step duration (seconds)
dx      = 1;       % cell length (feet)
dy      = 10;      % lane width (feet)
g_safe  = 6*dx;    % safety gap (cells)
p_min   = 20;      % minimum perception range (cells)
p_max   = 150;     % maximum perception range (cells)
l_down  = p_max;   % downstream edge zone length (cells)
T_stbl  = 5;       % stability countdown initial value (time steps)
u_brk   = 9;       % max speed at speed breaker (cells/step)
u_curv  = 27;      % max speed at road curvature (cells/step)

params = struct(...
    'dt', dt, ...
    'dx', dx, ...
    'dy', dy, ...
    'g_safe', g_safe, ...
    'p_min', p_min, ...
    'p_max', p_max, ...
    'l_down', l_down, ...
    'T_stbl', T_stbl, ...
    'u_brk', u_brk, ...
    'u_curv', u_curv);

%% -------------------------------------------------------
%  1.2 NETWORK TOPOLOGY
%% -------------------------------------------------------

% Physical segments
S_all  = [1, 2, 3, 4];
M_phys = length(S_all);          % number of physical segments
M      = M_phys;                 % M is used for loops over physical segments
K_count      = [4, 5, 3, 4];     % K(m) = lane count in segment m
L_seg  = [1262, 1735, 3200, 2092]; % segment lengths in cells

% Length = p_max cells (sufficient for perception range)
l_exit = p_max;                  % length of exit buffer (cells)
K_exit = K_count(end);           % same number of lanes as last physical segment

% Total segments including M_exit (for matrix allocation)
M_total = M_phys + 1;            % 4 physical + 1 exit buffer
K_all   = [K_count, K_exit];     % lane counts including exit buffer
l_all   = [L_seg, l_exit];       % lengths including exit buffer

fprintf('Physical segments: %d\n', M_phys);
fprintf('Exit buffer length: %d cells\n', l_exit);

%% -------------------------------------------------------
%  1.3 CONTINUATION MATRICES
%% -------------------------------------------------------

C = cell(1, M_phys);

% Segment 1 -> 2 : 4 lanes -> 5 lanes, expansion on right
C{1} = [1 1  0;
        2 2  0;
        3 3  0;
        4 4  0;
        0 5  1];

% Segment 2 -> 3 : 5 lanes -> 3 lanes, drops on left and right
C{2} = [1 0 -1;
        2 1  0;
        3 2  0;
        4 3  0;
        5 0 -1];

% Segment 3 -> 4 : 3 lanes -> 4 lanes, expansion on left
C{3} = [0 1  1;
        1 2  0;
        2 3  0;
        3 4  0];

% Segment 4 -> M_exit : 4 lanes, straight to exit buffer
C{4} = [1 1  0;
        2 2  0;
        3 3  0;
        4 4  0];


%% -------------------------------------------------------
%  1.4 SPEED-LIMIT AND FORBIDDEN-ZONE MATRICES
%% -------------------------------------------------------

% Allocate for physical segments + M_exit
max_lanes = max(K_all);
M_spd = cell(M_total, max_lanes);
F     = cell(M_total, max_lanes);

% --- Speed breaker cell locations (physical segments only) ---
sb = cell(M_total, max_lanes);
for k = 1:K_count(1)
    sb{1,k} = [1233, 1234, 1235, 1236];
end
for k = 1:K_count(2)
    sb{2,k} = [25, 26, 27, 28, 1372, 1373, 1374, 1375];
end
for k = 1:K_count(3)
    sb{3,k} = [];
end
for k = 1:K_count(4)
    sb{4,k} = [110, 111, 112, 278, 279, 280, 2090, 2091, 2092];
end
% M_exit has no speed breakers
for k = 1:K_exit
    sb{M_total, k} = [];
end

% --- Road curvature cell locations (physical segments only) ---
rc = cell(M_total, max_lanes);
for m = 1:M_phys
    for k = 1:K_count(m)
        rc{m,k} = [];
    end
end
for k = 1:K_count(2)
    rc{2,k} = 1477:1654;   % segment 2, all lanes
end
% M_exit has no curvatures
for k = 1:K_exit
    rc{M_total, k} = [];
end

% --- Inaccessible domains (physical segments only) ---
Inacc = cell(M_total, max_lanes);
for m = 1:M_phys
    for k = 1:K_count(m)
        Inacc{m,k} = [];
    end
end
Inacc{2,5} = 50:1411;      % segment 2, lane 5 only
% M_exit has no inaccessible domains
for k = 1:K_exit
    Inacc{M_total, k} = [];
end

% --- Build M_spd and F ---
for m = 1:M_total
    L = l_all(m);
    for k = 1:K_all(m)
        
        % Base speed for this lane
        if m <= M_phys
            % Physical segment: compute base speed
            v_base_kmph = max(60, 30 + 10*k);
            v_base = ceil(v_base_kmph * 1000 / 1097.28);
        else
            % M_exit: fixed speed limit 55 cells/step
            v_base = 55;
        end
        
        % Initialize matrices
        M_spd{m,k} = ones(1, L) * v_base;
        F{m,k}     = zeros(1, L);
        
        % Apply speed breaker constraints (physical segments only)
        if m <= M_phys
            for sb_pos = sb{m,k}
                if sb_pos < 1 || sb_pos > L, continue; end
                M_spd{m,k}(sb_pos) = u_brk;
                for nn = 1:20
                    idx = sb_pos - nn;
                    if idx < 1, break; end
                    M_spd{m,k}(idx) = min(M_spd{m,k}(idx), ...
                        u_brk + nn * floor((v_base - u_brk) / 20));
                end
                for nn = 1:15
                    idx = sb_pos + nn;
                    if idx > L, break; end
                    M_spd{m,k}(idx) = min(M_spd{m,k}(idx), ...
                        u_brk + nn * floor((v_base - u_brk) / 15));
                end
            end
        end
        
        % Apply road curvature constraints (physical segments only)
        if m <= M_phys
            for rc_pos = rc{m,k}
                if rc_pos < 1 || rc_pos > L, continue; end
                M_spd{m,k}(rc_pos) = min(M_spd{m,k}(rc_pos), u_curv);
                for nn = 1:20
                    idx = rc_pos - nn;
                    if idx < 1, break; end
                    M_spd{m,k}(idx) = min(M_spd{m,k}(idx), ...
                        u_curv + nn * floor((v_base - u_curv) / 20));
                end
                for nn = 1:15
                    idx = rc_pos + nn;
                    if idx > L, break; end
                    M_spd{m,k}(idx) = min(M_spd{m,k}(idx), ...
                        u_curv + nn * floor((v_base - u_curv) / 15));
                end
            end
        end
        
        % Apply inaccessible domain (physical segments only)
        if m <= M_phys
            for id_pos = Inacc{m,k}
                if id_pos < 1 || id_pos > L, continue; end
                F{m,k}(id_pos)     = 1;
                M_spd{m,k}(id_pos) = 0;
                for nn = 1:20
                    idx = id_pos - nn;
                    if idx < 1, break; end
                    ref = max(1, idx - 1);
                    M_spd{m,k}(idx) = min(M_spd{m,k}(idx), ...
                        M_spd{m,k}(ref) + nn * floor((v_base - M_spd{m,k}(ref)) / 20));
                end
                for nn = 1:15
                    idx = id_pos + nn;
                    if idx > L, break; end
                    ref = max(1, id_pos - 1);
                    M_spd{m,k}(idx) = min(M_spd{m,k}(idx), ...
                        M_spd{m,k}(ref) + nn * floor((v_base - M_spd{m,k}(ref)) / 15));
                end
            end
        end
        
    end % lane loop
end % segment loop

%% =========================================================
%  STAGE 2: DYNAMIC SETUP
%% =========================================================

%  ALLOCATE CORE MATRIX CELL ARRAYS

R = cell(M, max(K_count));     % Road matrix      : 2 x l_seg(m)
V = cell(M, max(K_count));     % Speed matrix     : 1 x l_seg(m)
A = cell(M, max(K_count));     % Acceleration     : 1 x l_seg(m)
G = cell(M, max(K_count));     % Gap matrix       : 1 x l_seg(m)
S = cell(M, max(K_count));     % State matrix     : 4 x l_seg(m)
M_ext = cell(M_total, max(K_all));   % NOTE: will be used later

%% =========================================================
%  STAGE 3: GENERATE VEHICLE DATABASE P_veh
%% =========================================================

% --- Class definitions ---
vehType  = ["CAR","BUS","CNG","MICROBUS","SUV","TRUCK","BIKE"];
vehCode  = [10, 20, 30, 40, 50, 60, 70];
perCent  = [25, 25, 15, 15,  5,  5, 10];
n_class = size(perCent, 2);
n_t     = 3000;

% Parameter ranges [min, max] — one row per class
l_range     = [15 18; 35 50; 10 12; 20 28; 16 20; 40 60;  6  8];
vmax_range  = [46 55; 23 28; 18 23; 28 37; 46 55; 18 23; 55 64];
amax_range  = [ 3  5;  1  2;  2  3;  2  3;  3  5;  1  2;  4  6];
dmax_range  = [ 4  6;  2  4;  3  5;  3  5;  4  6;  2  3;  5  8];
alpha_range = [0.30 0.55; 0.60 0.85; 0.20 0.40; 0.40 0.65;
               0.35 0.60; 0.20 0.45; 0.60 0.80];
tau_range   = [1.0 2.0; 1.5 2.5; 2.0 3.0; 1.5 2.5;
               1.0 2.0; 2.0 3.5; 0.8 1.5];

% ============================================================
% STAGE 4: Allocate P_veh
% ============================================================

% col 1: serial number
P_veh = zeros(n_t, 9);
P_veh(:,1) = (1:n_t)';     

% --- Assign class codes proportionally then shuffle ---
counts      = round((perCent / 100) * n_t);
counts(end) = n_t - sum(counts(1:end-1));   % ensure exact total
ordered     = repelem(vehCode, counts);
shuf_idx    = randperm(n_t);
P_veh(:,2)  = ordered(shuf_idx);           % temporary class codes

% --- Populate parameters row by row ---
for nv = 1:n_t

    class_code = P_veh(nv, 2);
    cls        = find(vehCode == class_code, 1);

    % Col 3: length
    P_veh(nv,3) = randi([l_range(cls,1),    l_range(cls,2)]);

    % Col 4: v_max
    P_veh(nv,4) = randi([vmax_range(cls,1),  vmax_range(cls,2)]);

    % Col 5: a_max
    P_veh(nv,5) = randi([amax_range(cls,1),  amax_range(cls,2)]);

    % Col 6: d_max — enforce d_max >= a_max
    d_lo = max(dmax_range(cls,1), P_veh(nv,5));
    d_hi = dmax_range(cls,2);
    if d_lo > d_hi
        d_lo = dmax_range(cls,1);
    end
    P_veh(nv,6) = randi([d_lo, d_hi]);

    % Col 7: alpha (2 decimal places)
    P_veh(nv,7) = round(alpha_range(cls,1) + ...
                  rand * (alpha_range(cls,2) - alpha_range(cls,1)), 2);

    % Col 8: tau (1 decimal place)
    P_veh(nv,8) = round(tau_range(cls,1) + ...
                  rand * (tau_range(cls,2)  - tau_range(cls,1)), 1);

    % Col 9: phi_ID = 0 (inactive at start)
    P_veh(nv,9) = 0;

    % Col 2: overwrite with unique 6-digit ID
    P_veh(nv,2) = (vehCode(cls) * 10000) + nv;

end % nv loop

fprintf('P_veh generated: %d vehicles.\n', n_t);

% --- Summary ---
fprintf('Class distribution:\n');
for c = 1:n_class
    id_min = vehCode(c) * 10000 + 1;
    id_max = vehCode(c) * 10000 + n_t;
    actual = sum(P_veh(:,2) >= id_min & P_veh(:,2) <= id_max);
    fprintf('  %s (code %d): %d vehicles\n', vehType(c), vehCode(c), actual);
end

fprintf('\nStage 1 complete.\n');

% clear up unnecessary variables from workspace
% clear c cls actual class_code d_hi d_lo  id_max id_min id_pos k M n_class 
% clear idx max_lanes nn nv ordered perCent rc rc_pos ref sb sb_pos shuf_idx
% clear u_curv u_brk v_base v_base_kmph vehCode vehType 
% clear amax_range dmax_range vmax_range l_range tau_range alpha_range


