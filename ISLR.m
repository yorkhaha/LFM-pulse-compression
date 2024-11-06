function islr = ISLR(s,t)
    s_norm = abs(s)/max(abs(s));
    pos_max = find(s_norm == max(s_norm));
    pos_left = islocalmin(s_norm(1:pos_max));
    t_left = t(pos_left);
    pos_left = find(t == t_left(end));
    pos_right = islocalmin(s_norm(pos_max:end));
    t_right = t(pos_max:end);
    t_right = t_right(pos_right);
    pos_right = find(t == t_right(1));
    s_norm_main = s_norm(pos_left:pos_right);
    P_main = sum(s_norm_main.^2);
    % P_main = sum(s_norm(pos_left:pos_right).^2);
    P_total = sum(s_norm.^2);
    islr = 10*log10((P_total-P_main)/P_main);
end