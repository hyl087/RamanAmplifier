%% 考虑泵浦光和信号光之间的能量转移+后向瑞利散射噪声
% 在原有拉曼泵浦上继续改进
% 假设后向瑞利散射噪声足够小，不影响原有泵浦光到信号光的转移
% 加入DWDM滤除后向瑞利散射噪声
% 瑞利散射总反射比/反射镜个数,50km打入连续光，探测后向瑞利散射光的一般信号大小，与仿真结果进行比较，获得后向瑞利散射反射镜反射系数
%% 双向双信号模拟
clear;
h=waitbar(0,"计算中，请稍候！");
c=3e8;
PumpLambda=1450e-9;
PumpF=c/PumpLambda;
omegaPumpF=2*pi*PumpF;
SignalLambda=1550e-9;
SignalF=c/SignalLambda;
omegaSignalF=2*pi*omegaPumpF;
Pp=-30;% 单位mW
Pr=31.5;% 单位mW
alpha_p=0.19;%入射光衰减系数
alpha_r=0.23;%拉曼泵浦光衰减系数
alpha_P=log(10)*alpha_p/10;
alpha_R=log(10)*alpha_r/10;
g_r=1.62;%拉曼增益
PP=10^(Pp/10);
PR=10^(Pr/10);
% %瑞利散射系数设定
% alpha_bs=-42.3;
% alpha_bs=10^(alpha_bs/10);
% S=1.5e-3;
% alpha_s=0.2;
% alpha_bs=S*alpha_s;
alpha_bs=2e-3;%反射镜的总系数
SignalPower=0.1;
SignalPower_back=0.1;
power=0.5;
snrpRecord1=zeros(1,50);
snrpRecord2=zeros(1,50);
for jj=1
ForwardPower=200;
BackPower=200;
SignalPower=power(jj);
SignalPower_back=power(jj);
distance=200;
u=BackPower*exp(-alpha_R*distance);

%% 单路信号放大情况模拟
while (1)
PP1=SignalPower;
PP2=0;
PR1=ForwardPower;
PR2=u;
opts = odeset('RelTol',1e-6);
[z,p] = ode45(@ramanPump, [0,distance], [PP1;PP2;PR1;PR2],opts);
plot(z,p)
if abs(p(end,4)-BackPower)<0.01*BackPower
break;
end
u=u+0.1*u*(BackPower-p(end,4))/p(end,4); 
end
v=0;
i=0;
%% test计算
% PP1=0.1;
% PP2=10;
% PR1=200;
% PR2=0; 
% opts = odeset('RelTol',1e-4);
% [z,p] = ode45(@ramanPump, [0,50], [PP1;PP2;PR1;PR2],opts);
%% 双路泵浦双向放大迭代开始
v=SignalPower_back*exp(-alpha_P*distance);
while(1)
    i=i+1;
    if(mod(i,1000)==0)
        disp(i)
    end
PP1=SignalPower;
PP2=v;
PR1=ForwardPower;
PR2=u; 
opts = odeset('RelTol',1e-6);
[z,p] = ode45(@ramanPump, [0,distance], [PP1;PP2;PR1;PR2],opts);
plot(z,p)
if (abs(p(end,4)-BackPower)<0.01*BackPower)&&(abs(p(end,2)-SignalPower_back)<0.01*SignalPower_back)
    break;
end

if abs(p(end,4)-BackPower)>0.01*BackPower
    u=u+u*(BackPower-p(end,4))/p(end,4);
end

if abs(p(end,2)-SignalPower_back)>0.01*SignalPower_back
    v=v+0.01*v*(SignalPower_back-p(end,2))/p(end,2);
end

end
figure(1)
plot(z,p)

%% 有泵浦正向信号瑞利散射光
% 后向泵浦光终端为z(ii)噪声不影响
% 前向泵浦光已知
% 后向散射光为正向泵浦信号光乘以瑞利散射系数
len=length(z);
noise=0;
disp("泵浦正向信号的瑞利散射光计算开始=======================================================")
noiseRecord=zeros(1,len);
for ii=2:len  
    fit=p(ii,1)*alpha_bs*(z(ii)-z(ii-1))/distance;
    PP2end=p(ii,4)*exp(-alpha_R*z(ii));
    PP1=0;
    PP2=0.01*fit;
    PR1=ForwardPower;
    PR2=PP2end;
    while 1
    opts = odeset('RelTol',1e-6);
    [zrb,prb] = ode45(@ramanPump, [0,z(ii)], [PP1;PP2;PR1;PR2],opts);
    if (abs(prb(end,4)-p(ii,4))<0.01*p(ii,4))&&(abs(prb(end,2)-fit)<0.01*fit)
        break;
    end

    if abs(prb(end,4)-p(ii,4))>0.01*p(ii,4)
        PR2=PR2+PR2*(p(ii,4)-prb(end,4))/(prb(end,4));  
    end
    if abs(prb(end,2)-fit)>0.01*fit
        PP2=PP2+0.01*PP2*(fit-prb(end,2))/prb(end,2);
    end
    end
    noise=noise+PP2;
    noiseRecord(ii)=PP2;
    disp(ii/len*100)
end
figure(2)
plot(z,noiseRecord)
disp("泵浦正向瑞利散射光计算结束=======================================================")
snrp=10*log10(p(1,2)/noise);
disp("有泵浦光的信噪比为：")
disp(snrp)
snrpRecord1(jj)=snrp;
%% 有泵浦反向信号的瑞利散射光
disp("泵浦反向信号的瑞利散射光计算开始=======================================================")
%以z(ii)为起点，前向泵浦光在z（ii）处的功率为初始光
%以z(ii)为起点，后向泵浦光在z（ii）处的功率不随散射光的引入发生改变
noiseRecord2=zeros(1,len);
noise2=0;
for ii=1:len-1    
    PR2end=p(ii,4);
    PP1=p(ii,2)*alpha_bs*(z(ii+1)-z(ii))/distance;
    PP2=0;
    PR1=p(ii,3);
    PR2=PR2end;
    while 1
    opts = odeset('RelTol',1e-6);
    [zrb,prb] = ode45(@ramanPump, [0,z(end)-z(ii)+1e-6], [PP1;PP2;PR1;PR2],opts);
    if (abs(prb(end,4)-BackPower)<0.001*BackPower)%&&(abs(prb(end,3)-p(end,3))<0.001*p(end,3))
        break;
    end

    if abs(prb(end,4)-BackPower)>0.001*BackPower
        PR2=PR2+PR2*(BackPower-prb(end,4))/prb(end,4);  
    end
%     if abs(prb(end,3)-p(end,3))>0.001*p(end,3)
%         PR1=PR1+PR1*(p(end,3)-prb(end,3))/prb(end,3);  
%     end
    end
    noise2=noise2+prb(end,1);
    noiseRecord2(ii)=prb(end,1);
end
figure(3)
plot(z,noiseRecord2)
disp("泵浦反向信号的瑞利散射光计算结束=======================================================")
snrp2=10*log10(p(end,1)/noise2);
disp("有泵浦光后端的信噪比为：")
disp(snrp2)
snrpRecord2(jj)=snrp2;

%% 双向无泵浦光信号光状态
disp("无泵浦正向信号的瑞利散射光计算开始=====================================================")
i=0;
v=1e-2*SignalPower_back;

while(1)
    i=i+1;
    if(mod(i,1000)==0)
        disp(i)
    end
 PP1=SignalPower;
PP2=v;
PR1=0;
PR2=0; 
opts = odeset('AbsTol',1e-9);
[znp,pnp] = ode45(@ramanPump, [0,distance], [PP1;PP2;PR1;PR2],opts);

if (abs(pnp(end,2)-SignalPower_back)<0.001*SignalPower_back)
    break;
end

if abs(pnp(end,2)-SignalPower_back)>0.001*SignalPower_back
    v=v+v*(SignalPower_back-pnp(end,2))/pnp(end,2);
end

end

% alpha_bs=0.0225;
%% 无泵浦正向瑞利散射光
len=length(znp);
noise3=0;
noiseRecord3=zeros(1,len);
figure(6)
for ii=2:len
    fit=pnp(ii,1)*alpha_bs*(znp(ii)-znp(ii-1))/distance;
    PP1=0;
    PP2=0.01*fit;
    PR1=0;
    PR2=0;
    while 1
    opts = odeset('RelTol',1e-5);
    [zrb,prb] = ode45(@ramanPump, [0,znp(ii)], [PP1;PP2;PR1;PR2],opts);
    if abs(prb(end,2)-fit)<0.001*fit
        break;
    end
    if abs(prb(end,2)-fit)>0.001*fit
        PP2=PP2+PP2*(fit-prb(end,2))/prb(end,2);
    end
    end
    noiseRecord3(ii)=PP2;
    noise3=noise3+PP2;
    disp(ii/len*100)
end
snrn=10*log10(pnp(1,2)/noise3);
disp("无泵浦瑞正向信号的利散射光计算结束====================================================")
disp("无泵浦光的前端信噪比为：")
disp(snrn)

%% 无泵浦反向信号的瑞利散射光
len=length(znp);
noise4=0;
noiseRecord4=zeros(1,len);
for ii=1:len-1
    PP1=pnp(ii,2)*alpha_bs*(znp(ii+1)-znp(ii))/distance;
    PP2=0;
    PR1=0;
    PR2=0;
    opts = odeset('RelTol',1e-5);
    [zrb,prb] = ode45(@ramanPump, [0,1e-6+distance-znp(ii)], [PP1;PP2;PR1;PR2],opts);
    noiseRecord4(ii)=prb(end,1);
    noise4=noise4+prb(end,1);
    disp(ii/len*100)
end
snrn2=10*log10(pnp(end,1)/noise4);
disp("无泵浦反向信号的瑞利散射光计算结束====================================================")
disp("无泵浦光的后端信噪比为：")
disp(snrn2)
waitbar(jj/50)
end
close(h)
figure(7)
plot(power,snrpRecord1)
hold on
plot(power,snrpRecord2)
xlabel("signal/mW")
ylabel("SNR/dB")
% function Ourshoot
% t0=0.5; tf=2; N=250; t=linspace(t0,tf,N);
% y0=ExaSol(t0); % initial value for first equation
% A1=3; y0(2)=A1; % initial guesses for second equation
% yf0=ExaSol(tf); % get final point value
% EPSShooting=10^-5; close all; figure; figure;
% col=str2mat('b:','g:','c:','m:','k:','b--','g--','c--');
% ya=ExaSol(t); figure(1); plot(t,ya(:,1),'r-'); hold on;
% figure(2); plot(t,ya(:,2),'r-'); hold on;
% for i=1:50, % maximum loop times
%  [y,t]=Ourrk4([t0,tf],y0,N);%y0(2) is modified each loop
%  if i==1,
%  yf1=y(end,2); %get last point value for second equation
%  flagEPS=abs((y(end,:)-yf0)./yf0); % relative error
%  A2=yf0(2).*A1./yf1; y0(2)=A2; % y0(2) is modified
%  else
%  yf2=y(end,2); %get last point value for second equation
%  flagEPS=abs((y(end,:)-yf0)./yf0); % relative error
%  A3=A1+(A2-A1)*(yf0(2)-yf1)/(yf2-yf1);
%  A1=A2; A2=A3; yf1=yf2; y0(2)=A2; % y0(2) is modified
%  end
%  figure(1); plot(t,y(:,1),col(i,:)); hold on;
%  figure(2); plot(t,y(:,2),col(i,:)); hold on;
%  if all(flagEPS<EPSShooting)==1, break; end
% end 