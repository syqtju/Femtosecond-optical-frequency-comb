c=3*10^8;
lam=1560;%单位nm
fc=c/(lam*10^(-9));
wc=2*pi*fc*10^(-12);%单位rad*ps^(-1)
fsigpul=2.746*10^12/2.355;%单位Hz
frep=50*10^6;%单位Hz
Trep=1/frep*10^12;%单位ps
fceo=100*10^6;
Tceo=1/fceo*10^12;%单位ps
wceo=2*pi*fceo*10^(-12);%单位同wc

%%%\频域建模%%%
% fn=10*10^4;
% fmax=fn*frep+fceo;
% fcn=5;
% df2=10^10;%单位Hz
% df1=frep/fcn;%单位Hz
fr=fc-3.5*10^12:frep:fc+3.5*10^12;
frn=length(fr);
fbao=(1/((2*pi)^0.5*fsigpul))*exp(-(fr-fc).^2/(2*fsigpul^2));

%%%%%%%%%%%%%%%%%%%%%变步长%%%%%%%%%%%%%%%%%%%%%%%%%%
ftmin=-450;
ftmax=450;
fdt=0.1;
ft1=ftmin:fdt:ftmax;%单位fs
ftn1=length(ft1);
fTTrep=ftn1+2;
N=5;
ft=zeros(1,N*fTTrep);
fTTTrep=20000*10^3;%单位fs
ftcdet=(fTTTrep-2*ftmax)/(fTTrep-ftn1);%单位fs,?f
 for ftcci=0:N-1%扩充
        for ftcct=ftcci*fTTrep+1:ftcci*fTTrep+ftn1
            ft(ftcct)=ftcci*fTTTrep-ftmax+(ftcct-ftcci*fTTrep)*fdt;
        end;
        for ftcct=ftcci*fTTrep+ftn1+1:(ftcci+1)*fTTrep
            ft(ftcct)=ftcci*fTTTrep+ftmax+(ftcct-ftcci*fTTrep-ftn1)*ftcdet;
        end;
end;
ftn=length(ft);
fE=zeros(1,ftn);
E=zeros(1,ftn);
ftfs=ft*10^(-15);
fii=1;
for frr=1:frn
    fE=fbao(frr)*cos(2*pi*fr(frr)*ftfs);%此时ft单位为fs。
    E=E+fE;
    fii=fii+1;
end;
fts=ft/1000;%单位ps
plot(fts,E);grid;
xlabel('时间/ps');
ylabel('振幅');
% axis([-0.5,0.5,min(E),max(E)]);
axis([79999.5,80000.5,min(E),max(E)]);
% E1=E(1004:2006);