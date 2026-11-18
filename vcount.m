v_count = size (P_veh, 1);
for nn = 1 : v_count
    id_max = 0;
    if P_veh(nn, 9) == 1
        id_max = nn;
    end
end
