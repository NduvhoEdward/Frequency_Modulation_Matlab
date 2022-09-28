%{
    1490804
    Nduvho E. Ramashia
    Practical Lab 2.
    Frequency Modulation
    22-Sep-2022
%}

clear; close all;

%% Initializaions
    fs = 1*10^6;                  % Sampling frequency
    T = 1/fs;                     % Sampling period
    L = 100000;                   % Length of signal
    t = 0:T:0.1; %(0:L-1)*T;                % Time vector

    fc = 1000; 
    Ac = 1;
    Bf = 5;
    
    m_t1 = 2*sinc(100*t) + 10.*t;
    m_t2 = 2*sinc(100*t) + (1 - 10.*t);

    m_t = m_t1.*(heaviside(t)-heaviside(t-0.05)) + ...
          m_t2.*(heaviside(t-0.05)-heaviside(t-0.1)); 
    
%% Message time and frequency representations
    plot(t, m_t); 
    grid;
    xlim([0 0.1])
%%
    M_f = fft(m_t);
    fshift = (-L/2: L/2 )*(fs/L);
    ushift = fftshift(M_f);
    plot(fshift,abs(ushift));
    title('Double-Sided Amplitude Spectrum of m(t)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([-2000 2000])

%%

%% Modulated Signal
    W = obw(m_t,fs); %1500;
    m_max = max(m_t);
    kf = Bf*W/m_max;
    intg_m = cumtrapz(t, m_t);
    u_t = Ac*cos(2*pi*fc*t + 2*pi*kf*intg_m);

%%
    plot(t, u_t);
    grid;
    xlim([0 0.1])
%%
    U_f = fft(u_t);
    fshift = (-L/2: L/2)*(fs/L);
    ushift = fftshift(U_f);
    plot(fshift,abs(ushift));
    title('Double-Sided Amplitude Spectrum of u(t)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([-1500 1500])


%%
    fi_t = fc + kf*m_t;
    plot(t,fi_t);
    grid on;

%%
    bwf = obw(fi_t,fs);
%%  

%%

%%
    %2*pi*fc*t + 2*pi*kf*intg_m
    cos_angle = acos(u_t);
    theta = cos_angle - 2*pi*fc*t;
%%

%%
    m_out = diff(theta);
    m_out = m_out/max(m_out);
    plot(t(:,1:length(m_out)),m_out);
    grid on; 
    %xlim([0 0.03])
%%
    Q = unwrap(m_out);
    plot(t(:,1:length(Q)),Q)

%%

%%

