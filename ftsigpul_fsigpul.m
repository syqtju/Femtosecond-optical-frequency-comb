c=3*10^8;
lam=1560;%单位nm
fc=c/(lam*10^(-9));
wc=2*pi*fc*10^(-12);%单位rad*ps^(-1)
% fsigpul=2.746*10^12/2.355;%单位Hz
frep=50*10^6;%单位Hz
Trep=1/frep*10^12;%单位ps
fceo=100*10^6;
Tceo=1/fceo*10^12;%单位ps
wceo=2*pi*fceo*10^(-12);%单位同wc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ftmin=-450;
ftmax=450;
fdt=1;
ft=ftmin:fdt:ftmax;%单位fs

ftn=length(ft);
fE=zeros(1,ftn);
E=zeros(1,ftn);
ftfs=ft*10^(-15);

ftsigpul=zeros(1,51);
ftsi=1;
for fsigpul=0.016*10^12:0.2*10^12:10.016*10^12
    %%%\频域建模%%%
    fr=fc-3.5*10^12:frep:fc+3.5*10^12;
    frn=length(fr);
    fbao=(1/((2*pi)^0.5*fsigpul))*exp(-(fr-fc).^2/(2*fsigpul^2));
    
    %%%%%%%%%%%%%%%%%%%合成%%%%%%%%%%%%%%%%%%%%%%%
    for frr=1:frn
        fE=fbao(frr)*cos(2*pi*fr(frr)*ftfs);%此时ft单位为fs。
        E=E+fE;
    end;
    %%%%%%%%%%%求包络%%%%%%%%%%%%%
    Ebao=zeros(1,ftn);
    for ftbao=1:ftn
        Ebao(ftbao)=E(ftbao)/(cos(2*pi*fc*ftfs(ftbao)));
    end;
    
    %%%%%%%%%%%%%求半宽高%%%%%%%%%%%%%%%
    tmaxhalf=0.5*max(Ebao);
    EEbao=abs(Ebao-tmaxhalf);
    half=(ftn-1)/2;
    [ezuo,tzuo]=min(EEbao(1:half));
    [eyou,tyou]=min(EEbao(half+1:ftn));
    ftsigpul(ftsi)=ft(half+tyou)-ft(tzuo);
    ftsi=ftsi+1;
end;
fsigpul=0.016*10^12:0.2*10^12:10.016*10^12;%单位Hz
%%%%%%%%%%%%曲线%%%%%%%%%%%%%%%%%%
% fts=ft/1000;%单位ps
% subplot(2,1,1)
% plot(fts,Ebao);grid;
% xlabel('时间/ps');
% ylabel('振幅');
% hold on;
% subplot(2,1,2)
plot(fsigpul,ftsigpul);grid;
xlabel('频域半高宽/Hz');
ylabel('时域脉冲宽度/fs');