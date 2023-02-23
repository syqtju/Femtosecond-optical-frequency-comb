c=3*10^8;
lam=1560;
fc=c/(lam*10^(-9));
wc=2*pi*fc*10^(-12)*10^(-3);%单位rad*fs^(-1)
sigpul=127.398;%单位fs
frep=50*10^6;%单位Hz
Trep=1/frep*10^12*10^3;%单位fs
fceo=20*10^6;%单位Hz
Tceo=1/fceo*10^12*10^3;%单位fs
wceo=2*pi*fceo*10^(-12)*10^(-3);%单位同wc
%%高斯包络与载波%%%%%
tfmax=500;
t0=-tfmax;tf=tfmax;dt=0.05;fs=1/dt;%dt决定时域采样频率为0.05fs；
t1=t0:dt:tf-dt;ts=0;
D=(1/(2*pi)^0.5*sigpul)*exp(-t1.^2/(2*sigpul^2));%飞秒脉冲高斯包络,横坐标fs
E=exp(wc*t1*1i);
F=D.*E;%单个飞秒脉冲,横坐标fs
FF=real(F);
f=length(F);
% plot(t1,FF);grid;
% axis([-200,200,min(FF),max(FF)]);
%%%%%%%%%%梳状函数%%%%%%%%%%
N=5;%飞秒激光脉冲个数；
NN=0;
phi0=2*pi*Trep/Tceo;
phi=NN*phi0;%相位移动
% TTrep=3*tfmax*fs;%单位fs
TTrep=2*tfmax*fs+10000;%单位fs
B=zeros(1,N*TTrep);%横坐标，单位为fs
for tt=1:N*TTrep%tt单位为fs
    NN=fix(tt/(TTrep));
    phi=NN*phi0;
    if(rem(tt,TTrep)==1)
        B(tt)=1*exp(-phi*1i);
    end;
end;
b=length(B); 
BB=real(B);
plot(BB);
%%%%%%%卷积，得时域上5个飞秒脉冲%%%%%%%%
C=conv(F,B);%单位为fs，此时高斯包络部分单位为0.01fs，高斯包络与高斯包络之间单位为fs
CC=C(1:b);%
c=length(CC);
tc=1:b;
TTTrep=20000*10^3;%单位fs
tcdet=(TTTrep-2*tfmax)/(TTrep-f);%单位fs
    for tcci=0:N-1%扩充
        for tcct=tcci*TTrep+1:tcci*TTrep+f
            tc(tcct)=tcci*TTTrep-tfmax+(tcct-tcci*TTrep)*dt;
        end;
        for tcct=tcci*TTrep+f+1:(tcci+1)*TTrep
            tc(tcct)=tcci*TTTrep+tfmax+(tcct-tcci*TTrep-f)*tcdet;
        end;
    end;
ttc=tc/1000;%单位ps
CCB=real(C(1:b));
plot(ttc,CCB);grid;
% axis([19999,20001,min(FF),max(FF)]);