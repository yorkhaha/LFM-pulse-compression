function irw = IRW(s_dB,trc)
    pos_max = find(s_dB == max(s_dB));
    pos_left = find(abs(s_dB(1:pos_max)+3) == min(abs(s_dB(1:pos_max)+3)));
    pos_right = find(abs(s_dB(pos_max:end)+3) == min(abs(s_dB(pos_max:end)+3)));
    pos_right = pos_max+pos_right-1;
    irw = trc(pos_right)-trc(pos_left);
end