% plot(ftt,ftimemag);
% axis([9999,10001,min(ftimemag),max(ftimemag)]);
% plot(fr,famp);
%%%%%%%%%%%%测定傅里叶变换得到的时域半高宽%%%%%%%%%%%%%%%%%
% fthw=[1 2];
% ftwn=1;
% ftmax=max(ftimemag);
% fthmax=0.5*ftmax;
% ftcom=0.15*10^(-8);
% ftbig=ffts/fcn;
% for ftwcir=1:ftbig%%%测定时域图半高宽
%     if (abs(ftimemag(ftwcir)-fthmax)<=ftcom)
%     fthw(ftwn)=ftwcir*fcn/(ffts*frep)*10^12;%单位ps
%     ftwn=2;
%     end;
% end;
% ftsigpul=(fthw(2)-fthw(1))*10^3;%单位fs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 所得一些数据：
% fsigpul=1.168*10^12;ftcom=0.05*10^(-8);ftsigpul=312fs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(ftt,ftimemag);grid;
% xlabel('时间/ps');
% ylabel('振幅');
% axis([49999.9,50000.1,min(ftimemag),max(ftimemag)]);%????
%%%%%%%%%%%%测定better时域半高宽%%%%%%%%%%%%%%%%%
% thw=[1 2];
% twn=1;
% tmax=max(CC);
% thmax=0.5*tmax;
% tcom=1;
% tbig=b/N;
% for twcir=1:tbig%%%测定时域图半高宽
%     if (abs(CC(twcir)-thmax)<=tcom)
%     thw(twn)=tc(twcir);%单位ps
%     twn=2;
%     end;
% end;
% tsigpul=(thw(2)-thw(1));%单位fs
%%%%%%%%%%%%%%%频谱半高宽特殊情况%%%%%%%%%%%%%%%%%%
% fsigpul=6*10^12，ftcom=0.2*10^(-8)，ftsigpul=58

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%测定没有better时的时域半高宽%%%%%%%%%%%%%%%%%
% thw=[1 2];
% twn=1;
% tmax=max(CC);
% thmax=0.5*tmax;
% tcom=0.001;
% tbig=b/N;
% for twcir=1:tbig%%%测定时域图半高宽
%     if (abs(CC(twcir)-thmax)<=tcom)
%     thw(twn)=tc(twcir);%单位ps
%     twn=2;
%     end;
% end;
% tsigpul=(thw(2)-thw(1))*10^3;%单位fs
%%%%%%%%%%%%测定没有better时的fft得到的傅里叶变换的半高宽%%%%%%%%%%%%%%%%%
% fhw=[1 2];
% fwn=1;
% fmax=max(mag);
% fhmax=0.5*fmax;
% fcom=0.1;
% for fwcir=1:NNN%%%测定时域图半高宽
%     if (abs(mag(fwcir)-fhmax)<=fcom)
%     fhw(fwn)=ff(fwcir);%单位ps
%     fwn=2;
%     end;
% end;
% ffsigpul=(fhw(2)-fhw(1));%单位10^14Hz
%%%%%%%%%%%%%%%%%%%测定傅里叶变换的得到的时域峰值坐标%%%%%%%%%%%%
% fengzhi=zeros(2,1000);
% cciii=1;
% for ccii=1:ffts
%     if(ftimemag(ccii)==2*10^(-8))
%         fengzhi(1,cciii)=ftimemag(ccii);
%         fengzhi(2,cciii)=ftt(ccii);
%         cciii=cciii+1;
%     end;
% end;
%%%%%%%%%%%测定傅里叶变换的得到的时域分辨率%%%%%%%%%%%%%%%%
% fengzhi=zeros(3,1000);
% jianju=zeros(1,1000);
% cciii=1;
% for ccii=1:ffts
%     if(ft(ccii)>49171&&ft(ccii)<49172)
%         fengzhi(1,cciii)=ftimemag(ccii);
%         fengzhi(2,cciii)=ft(ccii);
%         fengzhi(3,cciii)=ccii;
%         cciii=cciii+1;
%     end;
% end;
% cciii=1;
% for ccii=1:23
%     jianju(cciii)=fengzhi(2,ccii+1)-fengzhi(2,ccii);
%     cciii=cciii+1;
% end;
% plot(fengzhi(2,1:23),fengzhi(1,1:23))

c=3*10^8;
lam=1560;%单位nm
fc=c/(lam*10^(-9));
wc=2*pi*fc*10^(-12);%单位rad*ps^(-1)
% fsigpul=2.746*10^12/2.355;%单位Hz
fsigpul=0.001*10^12;%单位Hz
frep=50*10^6;%单位Hz
Trep=1/frep*10^12;%单位ps
fceo=100*10^6;
Tceo=1/fceo*10^12;%单位ps
wceo=2*pi*fceo*10^(-12);%单位同wc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ftmin=-6000;
ftmax=6000;
fdt=1;
ft=ftmin:fdt:ftmax;%单位fs

ftn=length(ft);
fE=zeros(1,ftn);
E=zeros(1,ftn);
ftfs=ft*10^(-15);

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
    ftsigpul=ft(half+tyou)-ft(tzuo);
%%%%%%%%%%%%曲线%%%%%%%%%%%%%%%%%%
fts=ft/1000;%单位ps
plot(fts,Ebao);grid;
xlabel('时间/ps');
ylabel('振幅');
axis([-8,8,min(Ebao),max(Ebao)]);