%{
    1490804
    Nduvho E. Ramashia
    Practical Lab 2.
    Frequency Modulation
    22-Sep-2022
%}

clear; close all;

%% Initializaions
    Fs = 1*10^6;                  % Sampling frequency
    T = 1/Fs;                     % Sampling period
    L = 100000;                   % Length of signal
    t = (0:L-1)*T;                % Time vector

    fc = 1000; 
    Ac = 1;
    Bf = 5;

    m_t1 = 2*sin(pi*100*t)./(pi*100*t) + 10.*t;
    m_t2 = 2*sin(pi*100*t)./(pi*100*t) + (1 - 10.*t);

    m_t = m_t1.*(heaviside(t)-heaviside(t-0.05)) + ...
          m_t2.*(heaviside(t-0.05)-heaviside(t-0.1)); 
    
    c_t = Ac*cos(2*pi*fc*t);

    plot(t, m_t); 
    grid

%% 
    




