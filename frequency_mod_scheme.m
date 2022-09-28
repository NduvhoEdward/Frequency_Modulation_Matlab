%{
    1490804
    Nduvho E. Ramashia
    Practical Lab 2.
    Frequency Modulation
    22-Sep-2022
%}

%% sub.Q1 -- Message Signal
    clear; close all;
    fs = 1*10^6;                  % Sampling frequency
    T = 1/fs;                     % Sampling period
    L = 100000;                   % Length of signal
    t = (0:L)*T;                  % Time vector

    fc = 1000; 
    Ac = 1;
    Bf = 5;
    
    m_t1 = 2*sinc(100*t) + 10.*t;
    m_t2 = 2*sinc(100*t) + (1 - 10.*t);

    m_t = m_t1.*(heaviside(t)-heaviside(t-0.05)) + ...
          m_t2.*(heaviside(t-0.05)-heaviside(t-0.1)); 
    
    figure(1);
    plot(t, m_t); 
    title('Massage Signal');
    xlabel('Time');
    ylabel('Amplitude');
    xlim([0 0.1]);
    grid;

    M_f = fft(m_t);
    fshift = (-L/2: L/2 )*(fs/L);
    ushift = fftshift(M_f);
    figure(2);
    plot(fshift,abs(ushift));
    title('Double-Sided Amplitude Spectrum of m(t)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([-200 200]);
    grid on;

%% Sub.Q2 -- Modulated Signal. 
    W = obw(m_t,fs); %1500;
    m_max = max(m_t);
    kf = Bf*W/m_max;
    intg_m = cumtrapz(t, m_t);
    u_t = Ac*cos(2*pi*fc*t + 2*pi*kf*intg_m);
    
    figure(3);
    plot(t, u_t);
    grid;
    title('Frequency Modulated Signal');
    xlabel('Time');
    ylabel('Amplitude');
    xlim([0 0.1])

    U_f = fft(u_t);
    fshift = (-L/2: L/2)*(fs/L);
    ushift = fftshift(U_f);
    figure(4);
    plot(fshift,abs(ushift));
    title('Amplitude Spectrum of u(t)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([-1500 1500])
    grid on;

%% Sub.Q3 -- Instanteneous Frequency
    fi_t = fc + kf*m_t;
    bwfi = obw(fi_t,fs);
    figure(5);
    plot(t,fi_t);
    grid on;
    title('Instanteneous Frequency');
    xlabel('Time');
    ylabel('Amplitude');
    grid on;

%% Sub.Q4 -- Demodulated Signal
    cos_angle = acos(u_t);
    theta = cos_angle - 2*pi*fc*t;
    m_out = diff(theta)/T;
    m_out(end+1)=m_out(end);
    figure(6);
    plot(t,m_out, 'DisplayName','Freq. and Amp. Modulated sig.');
    hold on;

    m_out = envelope(m_out,10,'rms');
    plot(t,m_out, 'DisplayName','AM Envelope')

    title('Demodulated Signal');
    xlabel('Time');
    ylabel('Amplitude');  

    legend;
    hold off;
    grid on;    
%% Sub.Q5 -- Demodulated Noisy Signal
    n_t = sqrt(0.05)*randn(size(t));
    y_t = u_t + n_t;
    figure(7);
    plot(t, y_t);
    grid on;
    title('Noisy Received Signal');
    xlabel('Time');
    ylabel('Amplitude');  

    cos_angle2 = acos(y_t);
    theta2 = real(cos_angle2) - 2*pi*fc*t;
    m_out2 = diff(theta2);
    m_out2 = m_out2/T;
    m_out2(end+1)=m_out2(end);
    figure(8);
    plot(t,m_out2,'DisplayName','Noisy Freq. & Amp. Modulated sig.'); 
    hold on;

    m_out2 = envelope(m_out2,100,'rms');
    plot(t,m_out2,'DisplayName','Noisy signal Envelope')
    grid on;
    hold off;
    xlim([0 0.01])

    title('Demodulated Signal');
    xlabel('Time');
    ylabel('Amplitude');  
%%