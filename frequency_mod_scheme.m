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
    t = (0:L-1)*T;                % Time vector

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
    fshift = (-L/2: L/2 -1)*(fs/L);
    ushift = fftshift(M_f);
    plot(fshift,abs(ushift));
    title('Double-Sided Amplitude Spectrum of m(t)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([-200 200])

%%

%%
    
%%    
   
