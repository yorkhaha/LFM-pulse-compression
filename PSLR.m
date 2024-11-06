function pslr = PSLR(s_dB)
    pks      = findpeaks(s_dB);
    pks_sort = sort(pks,'descend');
    % pks_max  = find(pks == max(pks));
    % pslr     = pks(pks_max-1);
    pslr     = pks_sort(2);
end