function out = Integ_fire()

step = 1;%0.1; %ms
Max_Time = 200;%500; %ms
Time_Vector = (0:step:Max_Time);
real_time = 1e-3*Time_Vector ;
E_rest = -70;
g_rest = 0.1;%0.025; %microS
Nb_Neurons = 1;
V = ones(Nb_Neurons,length(Time_Vector))*E_rest;
g_sra = zeros(Nb_Neurons,length(Time_Vector));
I = 5 ;%nA
C = 1;%0.25;
V_th = -63;
% V_spike = 30;
tau_sra = 100;
E_K = -75;

for t =1:length(Time_Vector)-1
    g_sra(1,t+1) = g_sra(1,t) - step/tau_sra*g_sra(1,t);
    V(1,t+1) = V(t) - step/C*(g_rest*(V(1,t)-E_rest) + g_sra(1,t)*(V(1,t)-E_K) - I );% + g_e(1,t)*(V(1,t)-E_e) + g_i(1,t)*(V(1,t)-E_i)    
    if V(1,t+1) >= V_th
        V(1,t+1) = E_rest;
        g_sra(1,t+1) = g_sra(1,t+1) + 0.04;
    end
        
end

figure
plot(real_time,V)
figure
plot(real_time,g_sra)