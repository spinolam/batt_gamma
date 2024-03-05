%% Mônica Spinola Felix Oct. 2023
%% Code for project of RUL control of Batteries


%Model parameters:



Ts=0.2; %sampling time (in seconds)
soc_min=1;
soc_max=99;
%% courbe of state of charge
figure(1)
subplot 211
I_discharge=9;
time_disc=9011-8445;
Cr=time_disc*100*I_discharge/3600/(100-soc_min); %Nominal Capacity in Ah (This decreases wrt age, from 90% to 65%)
plot(100-(1:time_disc).*100*I_discharge/3600/Cr) %3659*100*2.5/3600/(100-soc_min)
subplot 212
I_charge=-2.8;
time_disc=12847-9641;
%Nominal Capacity in Ah (This decreases wrt age, from 90% to 65%)
plot(0-(1:time_disc).*100*I_charge/3600/Cr) %3659*100*2.5/3600/(100-soc_min)
Cr=2.5;
%% voc
Voc_max_mes=8.4;
Voc_min_mes=6.14;
Eo=8.41;
K2=(Voc_max_mes-Eo)*-soc_max; %second term for modeling voltage variation wrt SoC  %
K1=(Voc_min_mes-Eo+K2/soc_min)/-log(100-soc_min); %parameter modeling the voltage variation wrt SoC %
%K1=0.1;

figure(2)

plot((soc_min:1:soc_max),(Eo-K1*log(100-(soc_min:1:soc_max))-K2./(soc_min:1:soc_max)))
grid on
title('Open-circuit voltage in function of SoC')
xlabel('SoC (%)')
ylabel('Voltage (V)')
% hold on
% Eoteste=8.35;
% K1teste=0.2038; %parameter modeling the voltage variation wrt SoC
% K2teste=1; 
% 
% 
% plot((soc_max:-1:soc_min),Eoteste-K1teste*log(100-(soc_max:-1:soc_min))-K2teste./(soc_max:-1:soc_min))
 K1teste=0.1038; %parameter modeling the voltage variation wrt SoC
%K1=K1teste;
 K2teste=6.3497; %second term for modeling voltage variation wrt SoC 
%K2=K2teste;
% plot((soc_max:-1:soc_min),Eoteste-K1teste*log(100-(soc_max:-1:soc_min))-K2teste./(soc_max:-1:soc_min))
% 
% hold off
%plot((soc_min:soc_max),Eo-K1*log(100-(soc_min:soc_max))-K2./(soc_min:soc_max))
%plot(Eo-K1*log(100-(soc_max:-0.18:soc_min))-K2./(soc_max:-0.18:soc_min))

%% resistance discharge
Vmax_mes=8.3;
Vmin_mes=5;

I_discharge=2.5;
R_soc_min=(Vmin_mes-Voc_min_mes)/-I_discharge
R_soc_max=(Vmax_mes-Voc_max_mes)/-I_discharge
R0=0.04; %Nominal internal resistance (Ohms)

K3=1*(-R0+R_soc_min)*soc_min; %additional internal resistance parameter. 
        
K4=1*(-R0+R_soc_max)*(100-soc_max);
figure(3)

plot((soc_max:-1:soc_min),R0+K3./(soc_max:-1:soc_min)+K4./(100-(soc_max:-1:soc_min)))
K3teste=6; %additional internal resistance parameter. 
%subplot 212
%plot((soc_max:-1:soc_min),I_discharge*(R0+K3./(soc_max:-1:soc_min)+K4./(100-(soc_max:-1:soc_min))))
grid on
title('Nominal resistance in function of SoC')
xlabel('SoC (%)')
ylabel('Resistance (Ohm)')


%% validation descharge
figure(4)
I_discharge=2.5;
%plot(vt_cycle-(R+K3./(soc_max:-0.18:soc_min)+K4./(100-(soc_max:-0.18:soc_min)))')
soc_vector=soc_max:-0.027:soc_min;
plot((Eo-K1*log(100-soc_vector)-K2./soc_vector)-(R0+K3./soc_vector+K4./(100-soc_vector))*I_discharge)
hold on
%plot(vt_cycle)
hold off
%% Temperature
Tmin = 24;
I_discharge=2.5;
Tmin_charge=32;
Tamb = Tmin;
Tmax=40 ;
c2 = (Tmax-Tamb)/R_soc_min/I_discharge^2
c2=5;
c2c=c2;
I_charge=-3;
(Tmin_charge-Tamb)/R_soc_max/I_charge^2
%c1 = (Tmax-c2*Tamb)/R_soc_min

%c2c = (Tmax-Tamb)/R_soc_minc/I_charge^2
%% resistance charge
Vmax_mes=8.6;
Vmin_mes=7;
I_charge=-3;
R_soc_minc=(Vmin_mes-Voc_min_mes)/-I_charge
R_soc_maxc=(Vmax_mes-Voc_max_mes)/-I_charge
Rc=0.03; %Nominal internal resistance (Ohms)


K3c=0.02*(-Rc+R_soc_minc)*soc_min; %additional internal resistance parameter. 

K4c=2*(-Rc+R_soc_maxc)*(100-soc_max);
figure(5)
soc_vector=soc_min:0.03:soc_max;

plot((soc_min:1:soc_max),Rc+K3c./(soc_min:1:soc_max)+K4c./(100-(soc_min:1:soc_max)))
figure(6)
plot(Eo-K1*log(100-soc_vector)-K2./soc_vector-((Rc+K3c./soc_vector+K4c./(100-soc_vector))*I_charge))
grid on


hold on
%plot(vt_cycle)
hold off
%%
(Eo-K1*log(100-soc_min)-K2./soc_min)-((Rc+K3c./soc_min+K4c./(100-soc_min))*I_charge);
(Eo-K1*log(100-soc_max)-K2./soc_max)-((Rc+K3c./soc_max+K4c./(100-soc_max))*I_charge);
close all

%%
figure(7)
gamma_est=out.gamma_est.signals(1).values;
mode_vector=out.gamma_est.signals(2).values;
%voltage_out=out.voltage(:,1);
%current_out=out.current(:,2);
cum_energy=out.cum_energy;
plot(gamma_est(mode_vector==-1))
%hold on
%plot(gamma_est(mode_vector==-1)-highpass(gamma_est(mode_vector==-1),0.0005,5))
plot(movmean(gamma_est(mode_vector==-1),1e5))
%plot(cum_energy(mode_vector==-1)/1e3/3600,Cr./movmean(gamma_est(mode_vector==-1),2e5))
hold off