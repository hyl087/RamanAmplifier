%% 只考虑泵浦光和信号光之间的能量转移
function dP = ramanPump(z, P)
%% 参数设定
% 单位km mW 
c=3e8;
PumpLambda=1450e-9;
PumpF=c/PumpLambda;
omegaPumpF=2*pi*PumpF;
SignalLambda=1550e-9;
SignalF=c/SignalLambda;
omegaSignalF=2*pi*omegaPumpF;
Pp=-36.3;%dbm
Pr=21.5;
alpha_p=0.19;%入射光衰减系数 %dB/Km
alpha_r=0.23;%拉曼泵浦光衰减系数
g_r=1.62;%拉曼增益
PP1=10^(Pp/10);
PR1=10^(Pr/10);
PP1=0;
PP2=0.1;
PR1=200;
PR2=0;
alpha_P=log(10)*alpha_p/10;
alpha_R=log(10)*alpha_r/10;
g_R=1000*log(10)*g_r/10;
g_R=0.75e-19;
Aeff=50e-18;
Kpeff=1;
Kreff=2;
%% 微分方程组
dP = [PP1;PP2;PR1;PR2];
dP(1) = g_R.*P(3)*P(1)/(Kpeff*Aeff)+g_R.*P(4).*P(1)/(Kreff*Aeff)-alpha_P.*P(1);
dP(2) = -(g_R.*P(3)*P(2)/(Kreff*Aeff)+g_R.*P(4).*P(2)/(Kpeff*Aeff)-alpha_P.*P(2));
dP(3) = (-omegaPumpF/omegaSignalF*g_R.*P(3).*P(1)/(Kpeff*Aeff)-omegaPumpF/omegaSignalF*g_R.*P(3).*P(2)/(Kreff*Aeff)-alpha_R.*P(3));
dP(4)=  -(-omegaPumpF/omegaSignalF*g_R.*P(4)*P(1)/(Kreff*Aeff)-omegaPumpF/omegaSignalF*g_R.*P(4).*P(2)/(Kpeff*Aeff)-alpha_R.*P(4));
end
