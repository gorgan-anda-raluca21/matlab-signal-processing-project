
% Nume si prenume:Gorgan Raluca 
%

clearvars
clc

%%IDNTIFICAREA SISTEMELOR DE ORDIN II folosind METODA REGRESIEI LINIARE
%% Magic numbers 
m = 5; 
n = 5; 

%% Process data (fixed, do not modify)
a1 = 2*(0.15+(m+n/20)/30)*(1000+n*300); %parametru care controlează viteza dinamicii sistemului
a2 = (1000+n*300);%parametru asociat dinamicii interne a procesului,
b0 = (2.2+m+n)/5.5;%câștigul static aproximativ al sistemului

rng(m+10*n) %generator de nr aleatoare 
x0_slx = [(-1)^n*(-m/10-rand(1)*m/5); (-1)^m*(n/20+rand(1)*n/100)]; %cond initiale 

%% Experiment setup (fixed, do not modify)
Ts = 20/a1/1e4; % fundamental step size
Tfin = 36/a1; % simulation duration

gain = 15; %amplificare maxima pe intrare 
umin = -gain; umax = gain; % input saturation
ymin = -b0*gain/1.8; ymax = b0*gain/1.8; % output saturation

whtn_pow_in = 1e-9*5*(((m-1)*5+n)/5)/2; % input white noise power and sampling time
whtn_Ts_in = Ts*3; %putere zgomot alb input
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(9); % input quantizer (DAC)

whtn_pow_out = 1e-8*5*(((m-1)*8+n)/5)/2; % output white noise power and sampling time
whtn_Ts_out = Ts*5;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(9); % output quantizer (ADC)

u_op_region = -(m+n/5)/2; % operating point

%% Input setup (can be changed/replaced/deleted)
u0 = 0; %%valoarea inițială a semnalului de intrare
ust = 5;%amplitudinea variației de tip treaptă aplicată pe intrare
t1 = 12/a1; % momentul de timp la care este aplicată treapta 

%% Data acquisition (use t, u, y to perform system identification)
out = sim("circuit_electric_R2022b1.slx");

t = out.tout;%vectorul de timp al simulării
u = out.u;
y = out.y;

plot(t,u,t,y)
shg

%% System identification
%K=? 
 i1=4811 ;
 i2=5845  ;
 i3=11095 ;
 i4=11820;
 u0=mean(u(i1:i2));
 u_st=mean(u(i3:i4));
 y0=mean(y(i1:i2));
 y_st=mean(y(i3:i4));
 K=(y_st-y0)/(u_st-u0)

 %%
 % T2vec= 0.1:0.1:3.5;
 %Y_ecuatie = T1*T2vec.log(T2vec)-T2vec(Ti*T1*log(T1))+T1*Ti 
    i5=6001;
    i6=11963;
    t_aux=t(i5:i6);
    y_aux=y(i5:i6);
    figure
    plot(t_aux,abs(y_aux-y_st))
    %%
    i7= 25;
    i8=1100;
    i9=2155;
    treg=t_aux([i7 i8 i9])
    yreg=log(abs(y_aux([i7,i8,i9])-y_st))

    figure
    plot(treg,yreg)

    Areg=[sum(treg.^2),sum(treg);
        sum(treg),length(treg)]
    breg=[sum(yreg.*treg);sum(yreg)];
    theta=inv(Areg)*breg

    Re=theta(1)
    %%imag
    i10=5991;
    i11=7091;
    Tosc=2*(t(i11)-t(i10))
    Im=2*pi/Tosc

    %%
    wn=sqrt(Re^2+Im^2)
    zeta=-Re/wn
    % Calculate
    %%
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