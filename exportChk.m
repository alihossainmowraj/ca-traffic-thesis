%% exportChk.m  -  write chk to an Excel file
%  Run after runSim2.   Needs MATLAB R2019a or newer (writecell).

global chk

file_name       = 'chk_log.xlsx';
keep_empty_rows = false;     % true: keep all 1000 rows, even vehicles never on the road

if isempty(chk)
    error('chk is empty. Run initChk before runSim2.');
end

if exist(file_name, 'file')
    delete(file_name);       % avoid stale columns from an earlier, longer run
end

% ---- rows to export ----
n_chk = size(chk.vehicle_position, 1);
ID_r  = (1:n_chk)';
if keep_empty_rows
    rows = ID_r;
else
    rows = ID_r(any(~isnan(chk.vehicle_position), 2));
end

% =========================================================
%  Sheet 1: vehicle_position   [m k x] per time step
% =========================================================
n_t   = size(chk.vehicle_position, 2) / 3;
names = {'m', 'k', 'x'};
hdr   = cell(1, 3 * n_t);
for t = 0:n_t-1
    for c = 1:3
        hdr{3*t + c} = sprintf('t%d_%s', t, names{c});
    end
end
data = num2cell(chk.vehicle_position(rows, :));
data(cellfun(@(v) isnan(v), data)) = {[]};          % NaN -> empty cell
out  = [{'ID_r'}, hdr; num2cell(rows), data];
writecell(out, file_name, 'Sheet', 'vehicle_position');

% =========================================================
%  Sheet 2: shft_state   [f_lshift sigma] per time step
% =========================================================
n_t   = size(chk.shft_state, 2) / 2;
names = {'flshift', 'sigma'};
hdr   = cell(1, 2 * n_t);
for t = 0:n_t-1
    for c = 1:2
        hdr{2*t + c} = sprintf('t%d_%s', t, names{c});
    end
end
data = num2cell(round(chk.shft_state(rows, :) * 100) / 100);   % sigma to 2 decimals
data(cellfun(@(v) isnan(v), data)) = {[]};
out  = [{'ID_r'}, hdr; num2cell(rows), data];
writecell(out, file_name, 'Sheet', 'shft_state');

fprintf('chk exported to %s  (%d vehicles, %d time steps incl. t0)\n', ...
    file_name, numel(rows), n_t);