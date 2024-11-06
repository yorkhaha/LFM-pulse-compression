%% LFM信号脉压
% 2024/10/27 York Yuan
clear;clc;close all;

%% ---------------------------- 参数设置 ----------------------------------
RCS       = 1;                          % 反射系数
R_target  = [15,16]*1e3;                % 目标距离
alpha_os  = [0.8,1.0,1.2,1.4];          % 过采样率
TBP       = [25,50,100,200,400];        % 时间带宽积
fc        = 5e5;                        % 载频
Tr        = 30e-6;                      % 脉宽
Br        = 30e6;                       % 带宽
tau_start = 0.1e-3;                     % 接收开始时间
tau_win   = 0.04e-3;                    % 接收窗长
% c       = 299792458;                    % 光速
c         = 3e8;                        % 光速

%% --------------------------- LFM信号分析 --------------------------------
%% 基本参数
B  = TBP(3)/Tr;
K  = B/Tr;                              % 调斜率
fs = alpha_os(3)*B;                     % 采样率
dt = 1/fs;                              % 采样间隔
t  = -Tr/2:dt:Tr/2-dt;                  % 时间
N  = length(t);                         % 采样点数
f  = (-N/2:N/2-1)/N*fs;                 % 频率

%% LFM仿真
s = exp(1j*pi*K*t.^2);                  % LFM信号

figure;
plot(t,real(s));
xlabel('\it时间/秒');ylabel('\it幅度');title('LFM实部');
figure;
plot(t,imag(s));
xlabel('\it时间/秒');ylabel('\it幅度');title('LFM虚部');

%% LFM频谱
S = fftshift(fft(ifftshift(s)));        % LFM频谱

figure;
plot(f,abs(S));
xlabel('\it频率/赫兹');ylabel('\it幅度');title('LFM幅度谱');

%% 不同过采样率下的DFT结果
for i = 1:length(alpha_os)
    fsi = B*alpha_os(i);
    dti = 1/fsi;
    ti  = -Tr/2:dti:Tr/2-dti;
    Ni  = length(ti);
    fi  = (-Ni/2:Ni/2-1)/Ni*fsi;
    si  = exp(1j*pi*K*ti.^2);
    Si  = fftshift(fft(ifftshift(si)));
    figure;
    plot(fi,abs(Si));
    xlabel('\it频率/赫兹');ylabel('\it幅度');
    title(['过采样率为',num2str(alpha_os(i)),'的LFM幅度谱']);
end

%% 不同TBP的LFM信号的频谱
for i = 1:length(TBP)
    Bi  = TBP(i)/Tr;
    Ki  = Bi/Tr;
    fsi = alpha_os(3)*Bi;
    dti = 1/fsi;
    ti  = -Tr/2:dti:Tr/2-dti;
    Ni  = length(ti);
    fi  = (-Ni/2:Ni/2-1)/Ni*fsi;
    si  = exp(1j*pi*Ki*ti.^2);
    Si  = fftshift(fft(ifftshift(si)));
    figure;
    plot(fi,abs(Si));
    xlabel('\it频率/赫兹');ylabel('\it幅度');
    title(['TBP为',num2str(TBP(i)),'的LFM幅度谱']);
end

%% --------------------------- 脉冲压缩仿真 ------------------------------
%% 实现无误差的脉冲压缩
h = conj(s);                                                % 匹配滤波器
H = fftshift(fft(ifftshift(h)));                            % 滤波器频谱
Sc = S.*H;                                                  % 时域卷积频域相乘
Sc = ifftshift(Sc);
Sc = [Sc(1:end/2),zeros(1,29*length(Sc)),Sc(end/2+1:end)];  % 30倍上采样
frc = (-30*N/2:30*N/2-1)/(30*N)*30*fs;                      % 频率
sc = fftshift(ifft(Sc));                                    % 脉压信号
dtrc = 1/(30*fs);                                           % 时间间隔
trc = linspace(-Tr/2,Tr/2,30*N);                            % 时间

figure;
plot(frc,abs(fftshift(Sc)));
xlabel('\it频率/赫兹');ylabel('\it幅度');
title('脉冲压缩后的幅度谱');

figure;
plot(trc,sc);
xlabel('\it时间/秒');ylabel('\it幅度');axis tight;
title('脉冲压缩信号');

sc_dB = 20*log10(abs(sc)/max(abs(sc)));                     % 对数化
disp("无误差脉冲压缩：");
irw   = IRW(sc_dB,trc)                                      % 主瓣宽度
pslr  = PSLR(sc_dB)                                         % 峰值旁瓣比
islr  = ISLR(sc,trc)                                        % 积分旁瓣比

%% 观察频域加窗的影响
% 矩形窗
Sc_uniform_weighting    = Sc;
sc_uniform_weighting    = fftshift(ifft(Sc_uniform_weighting));
sc_uniform_weighting_dB = 20*log10(abs(sc_uniform_weighting)/max(abs(sc_uniform_weighting)));

figure;
plot(trc,sc_uniform_weighting);
xlabel('\it时间/秒');ylabel('\it幅度');axis tight;
title('矩形脉压');

disp("矩形窗脉冲压缩：");
irw  = IRW(sc_uniform_weighting_dB,trc)
pslr = PSLR(sc_uniform_weighting_dB)
islr = ISLR(sc_uniform_weighting,trc)

% 汉宁窗
B_left                  = length(frc)/2+1-round(B/2/(fs/N));
B_right                 = length(frc)/2+1+round(B/2/(fs/N));
H_hanning               = [zeros(1,B_left-1),(hann(B_right-B_left+1)).',zeros(1,length(frc)-B_right)];
Sc_hanning_weighting    = fftshift(Sc).*H_hanning;
sc_hanning_weighting    = fftshift(ifft(ifftshift(Sc_hanning_weighting)));
sc_hanning_weighting_dB = 20*log10(abs(sc_hanning_weighting)/max(abs(sc_hanning_weighting)));

figure;
plot(trc,sc_hanning_weighting);
xlabel('\it时间/秒');ylabel('\it幅度');axis tight;
title('汉宁窗脉压');

disp("汉宁窗脉冲压缩：");
irw  = IRW(sc_hanning_weighting_dB,trc)
pslr = PSLR(sc_hanning_weighting_dB)
islr = ISLR(sc_hanning_weighting,trc)

% 汉明窗
H_hamming               = [zeros(1,B_left-1),(hamming(B_right-B_left+1)).',zeros(1,length(frc)-B_right)];
Sc_hamming_weighting    = fftshift(Sc).*H_hamming;
sc_hamming_weighting    = fftshift(ifft(ifftshift(Sc_hamming_weighting)));
sc_hamming_weighting_dB = 20*log10(abs(sc_hamming_weighting)/max(abs(sc_hamming_weighting)));

figure;
plot(trc,sc_hamming_weighting);
xlabel('\it时间/秒');ylabel('\it幅度');axis tight;
title('汉明窗脉压');

disp("汉明窗脉冲压缩：");
irw  = IRW(sc_hamming_weighting_dB,trc)
pslr = PSLR(sc_hamming_weighting_dB)
islr = ISLR(sc_hamming_weighting,trc)

% 布莱克曼窗
H_blackman               = [zeros(1,B_left-1),(blackman(B_right-B_left+1)).',zeros(1,length(frc)-B_right)];
Sc_blackman_weighting    = fftshift(Sc).*H_blackman;
sc_blackman_weighting    = fftshift(ifft(ifftshift(Sc_blackman_weighting)));
sc_blackman_weighting_dB = 20*log10(abs(sc_blackman_weighting)/max(abs(sc_blackman_weighting)));

figure;
plot(trc,sc_blackman_weighting);
xlabel('\it时间/秒');ylabel('\it幅度');axis tight;
title('布莱克曼窗脉压');

disp("布莱克曼窗脉冲压缩：");
irw  = IRW(sc_blackman_weighting_dB,trc)
pslr = PSLR(sc_blackman_weighting_dB)
islr = ISLR(sc_blackman_weighting,trc)

% 凯瑟窗
H_kaiser               = [zeros(1,B_left-1),(kaiser(B_right-B_left+1,2.5)).',zeros(1,length(frc)-B_right)];
Sc_kaiser_weighting    = fftshift(Sc).*H_kaiser;
sc_kaiser_weighting    = fftshift(ifft(ifftshift(Sc_kaiser_weighting)));
sc_kaiser_weighting_dB = 20*log10(abs(sc_kaiser_weighting)/max(abs(sc_kaiser_weighting)));

figure;
plot(trc,sc_kaiser_weighting);
xlabel('\it时间/秒');ylabel('\it幅度');axis tight;
title('凯瑟窗脉压');

disp("凯瑟窗脉冲压缩：");
irw  = IRW(sc_kaiser_weighting_dB,trc)
pslr = PSLR(sc_kaiser_weighting_dB)
islr = ISLR(sc_kaiser_weighting,trc)

%% ------------------------- LFM回波仿真 -----------------------------------
%% 回波生成
fs = Br*alpha_os(3);
Nr = ceil((tau_win+Tr/2)*fs)+1;                         % 时间窗点数
% Nr = ceil(tau_win*fs)+1;
tr = linspace(tau_start-Tr/2,tau_start+tau_win,Nr);     % 时间窗
% tr = linspace(tau_start,tau_start+tau_win,Nr);
sr = zeros(1,Nr);                                       % 接收信号
K  = Br/Tr;
for i = 1:length(R_target)
    taui = 2*R_target(i)/c;                             % 时延
    ti = tr-taui;
    phasei = pi*K*ti.^2+2*pi*fc*ti;                     % 相位
    sr = sr+RCS*exp(1j*phasei).*(abs(ti)<=Tr/2);
end
figure;
plot(tr,abs(sr));
xlabel('\it时间/秒');ylabel('\it幅度');axis tight;
title('原始信号');

%% 完成脉冲压缩
sr  = sr.*exp(-1j*2*pi*fc*tr);                                  % 解调
thr = (-(Nr-1)/2:(Nr-1)/2)/Nr*(tau_win+Tr/2);                   % 滤波器时间
% thr = (-(Nr-1)/2:(Nr-1)/2)/Nr*tau_win;
hrc = ifftshift(conj(exp(1j*pi*K*thr.^2)).*(abs(thr)<=Tr/2));   % 匹配滤波器
src = ifft(ifftshift(fftshift(fft(sr)).*fftshift(fft(hrc))));   % 脉压信号

figure;
plot(c*tr/2,abs(src));
xlabel('\it距离/米');ylabel('\it幅度');axis tight;
title('脉冲压缩');