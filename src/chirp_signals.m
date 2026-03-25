%%
% Nume si prenume: Gorgan Anda Raluca 
%

clearvars
clc

%%IDNTIFICAREA SISTEMELOR DE ORDIN II folosind SEMNALE CHIRP
%% Magic numbers (replace with received numbers)
m = 5; 
n = 5; 

%% Process data (fixed, do not modify)
a1 = 2*(0.15+(m+n/20)/30)*(1000+n*300);
a2 = (1000+n*300);
b0 = (2.2+m+n)/5.5;

rng(m+10*n)
x0_slx = [(-1)^n*(-m/10-rand(1)*m/5); (-1)^m*(n/20+rand(1)*n/100)];

%% Experiment setup (fixed, do not modify)
Ts = 20/a1/1e4; % fundamental step size
Tfin = 36/a1*10; % simulation duration

gain = 15;
umin = -gain; umax = gain; % input saturation
ymin = -b0*gain/1.8; ymax = b0*gain/1.8; % output saturation

whtn_pow_in = 1e-9*5*(((m-1)*5+n)/5)/2; % input white noise power and sampling time
whtn_Ts_in = Ts*3;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(9); % input quantizer (DAC)

whtn_pow_out = 1e-8*5*(((m-1)*8+n)/5)/2; % output white noise power and sampling time
whtn_Ts_out = Ts*5;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(9); % output quantizer (ADC)

u_op_region = -(m+n/5)/2; % operating point

%% Input setup (can be changed/replaced/deleted)
wf=2000;%~wosc , usor decitit din step
Fmin=wf/2/pi/10;
Fmax=wf/2/pi*3;
Ain=1.5;
%% Data acquisition (use t, u, y to perform system identification)
out = sim("circuit_electric_R2025a1.slx");

t = out.tout;
u = out.u;
y = out.y;
framelen = 151;   % trebuie IMPAR,cate puncte foloseste filtru
ord = 3;     % ordinul pol
y_raw = y;   % păstram originalul
y = sgolayfilt(y_raw, ord, framelen);

%% System identification
Ay=(-2.5+9.7)/2;
Au=1.5;
K=Ay/Au;
K= 2.2250;
%%
wr=pi/(0.0134252-0.0120043)%pulsatia de rezonanta(determinarea frecvenței dintr-un interval de timp dintre două repere (ex: perioadă/semiperioadă) în zona de rezonanță)
Ay=-1.48+9.89;%amplitudinea maximă a semnalului de ieșire în zona de rezonanță
Au=1.5;%amplitudinea semnalului de intrare de tip Chirp
Mr=Ay/Au;%factorul de amplificare la rezonanță
r=roots([4*Mr^2,0,-4*Mr^2,0,K^2])%rădăcinile ecuației caracteristice asociate metodei vârfului de rezonanță,
zeta=0.327;%factorul de amortizare al sistemului
wn = wr / sqrt(1 - 2*zeta^2)%pulsația naturală a sistemului, calculată din pulsația de rezonanță și factorul de amortizare


A=[0 1 ; -wn^2, -2*zeta*wn]
    B=[0;K*wn^2];
    C=[1,0];
    D=0;

    sys=ss(A,B,C,D)
    ysim2=lsim(sys,u,t,[y(1),0]);

    figure
    plot(t,u,t,y,t,ysim2)

    J=1/sqrt(length(t)*norm(y-ysim2))
    empn=norm(y-ysim2)/norm(y-mean(y))*100
