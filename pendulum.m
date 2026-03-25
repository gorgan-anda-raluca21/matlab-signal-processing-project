%%
% Nume si prenume: Gorgan Raluca 
%

clearvars
clc

%% Magic numbers (replace with received numbers)
m = 5;
n = 5;

%% Process data and experiment setup (fixed, do not modify)
Ts = 500e-6; % fundamental step size

u_star = 1.2+n*0.075;
delta = 0.125;
delta_spab = 0.075;

umin = -5; umax = 5; % input saturation
ymin = -100; ymax = 100; % output saturation

g = 9.81;
% pendulum parameters
M = 0.8-n/48;
l = 1.2-m/24;
b = 0.3+m/24;
% measurement
c1 = 180/pi;
c2 = 4+n/2;

% (theta0,omega0)
rng(m+10*n)
x0_slx = [(n+3)/50,(-1)^(n+1)*m/20];

% input white noise power and sampling time
whtn_pow_in = 1e-10*(Ts*1e4)/2; 
whtn_Ts_in = Ts*2;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(13); % input quantizer (DAC)

% output white noise power and sampling time
whtn_pow_out = 1e-3*Ts; 
whtn_Ts_out = Ts*2;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(13); % output quantizer (ADC)

meas_rep = round(7+n/2); % data acquisition hardware sampling limitation

%% Input setup (can be changed/replaced/deleted)
%%CALIBRARE SPAB
 Tfin = 50; % simulation duration
 t0 = Tfin/4;
 t1 = Tfin/2;
 tu=0.2;
 tpo=1.5;
%durata necesara stabilizarii dupa CI nenule
t1=16;
 %calibrarea SPAB
 N=4;
 p=round(tpo/N/Ts);

 DeltaT=1.2*(2^N-1)*p*Ts*3;

 [input_LUT_dSPACE,Tfin] = generate_input_signal(Ts,t1,DeltaT,N,p,u_star,delta,delta_spab);


%% Data acquisition (use t, u, y to perform system identification)
out = sim("pendul_R2022b1.slx");

t = out.tout;
u = out.u;
y = out.y;

subplot(211)
plot(t,u)
subplot(212)
plot(t,y)
shg

%% System identification

i1= 52614;
i2=97060;
i3=112933;
i4=164876;

N=22;

t_id=t(i1:N:i2);
u_id=u(i1:N:i2);
u_id_new=u_id-mean(u_id);
y_id=y(i1:N:i2);
y_id_new=y_id-mean(y_id);

t_vd=t(i3:N:i4);
u_vd=u(i3:N:i4);
u_vd_new=u_vd-mean(u_vd);
y_vd=y(i3:N:i4);
y_vd_new=y_vd-mean(y_vd);
%%
figure
subplot(2,2,1)
plot(t_id,u_id_new)

subplot(2,2,3)
plot(t_id,y_id_new)

subplot(2,2,2)
plot(t_vd,u_vd_new)

subplot(2,2,4)
plot(t_vd,y_vd_new)
%%
dat_id=iddata(y_id_new,u_id_new,t_id(2)-t_id(1));
dat_vd=iddata(y_vd_new,u_vd_new,t_vd(2)-t_vd(1));
 
model_armax=armax(dat_id,[2,2,7,1]) %autocorelatie
figure 
resid(model_armax,dat_vd)
figure 
compare(model_armax,dat_vd)

%%
model_bj = bj(dat_id,[2,2,1,10,1]) % intercorelatia
figure 
resid(model_bj,dat_vd)
figure 
compare(model_bj,dat_vd)

%%
model_n4sid=n4sid(dat_id,1:15)
figure
resid(model_n4sid,dat_vd)
figure 
compare(model_n4sid,dat_vd)

%%
model_ssest=ssest(dat_id,2)
figure
resid(model_ssest,dat_vd)
figure 
compare(model_ssest,dat_vd)
 