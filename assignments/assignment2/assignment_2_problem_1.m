close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo

%%

%Ship parameters
L = 50;
B = 18;
D = 6;
T = 2.3;
C = 0.55;
d = 570;
I = 0.8175;

% other parameters and constants
g = 9.81; %gravity acc
c = 0.5772; %Euler constant 
M_NA = 2730e-3; %m possition of neutral axis 
S_NA = 0.05; %m

% Design parameters

Ar = 0.697*d/T;
ls = Ar/B;
%interpolating for ncg from table 1 in assignment descr.
        x1 = 400;
        x2 = 1200;
        x3 = 570;
        y1 = 2;
        y2 = 1;
        y3 = (y2-y1)/(x2-x1)*(x3-x1)+y1;

ncg = y3 - 1;

k1 = 110e+3;
k2 = 190e+3;

if L >= 90
    C1 = 0.044 * L + 3.75;
    
else 
    C1 = 10.75 - ((300-L)/100)^(1.5);
end

C2 = 0.01;
C3 = 1.25;
fp = 17.5;


% Statistical properties of the basic design random variables for the composite made ship
Mms_sag = 0; 
Mms_hog = (0.375 *fp*C1*C2*L^2*B*(C+0.7))*1e+3;
                                                   %Stillwater 
Sms_sag = 1.5e+6;
Sms_hog = 1.5e+6;                               

Mmw_sag = -k1*C1*L^2*B*(C+0.7);
Mmw_hog = k1*C1*L^2*C;
                                                    %Wave 
Smw_sag = 1.5e+6;
Smw_hog = 1.5e+6;

Mms1 = C3*d*(1+ncg)*(L-ls)*1e+3;    
Sms1 = 1e+7;%Slammin    g 

MSt = 249e+6;   %lognomral tensile                  
SSt = 40e+6;    

MSc = 213e+6;   %Lognomral Compressive
SSc = 30e+6;

Meta_u = 1;         %Un..
Seta_u = 0.15; 

Meta_s = 1;
Seta_s = 0.05;

Meta_w = 0.9;
Seta_w = 0.15;

Meta_s1 = 1;
Seta_s1 = 0.15;     %Uncertainty parameters

mphi = 0.9;
SPHI = 0;

%Random values

N = 100; %number of geneterated rnd values for each.. 

        %random lognormal for tensile
a = log((MSt^2)/sqrt((SSt^2)+MSt^2)) ;
b = sqrt(log((SSt^2)/(MSt^2)+1));
ST = lognrnd(a,b,N,1);


        %random lognormal for compression
a = log((MSc^2)/sqrt((SSc^2)+MSc^2)) ;
b = sqrt(log((SSc^2)/(MSc^2)+1));
SC = lognrnd(a,b,N,1);
        %Normal distributed random data 
etaU = wnormrnd(Meta_u,Seta_u^2,N,1);

etaS = wnormrnd(Meta_s,Seta_s^2,N,1);

etaW = wnormrnd(Meta_w,Seta_w^2,N,1);

etaS1 = wnormrnd(Meta_s1,Seta_s1^2,N,1);


    %Normal distributed random value for still water hog/sag
MS_hog = wnormrnd(Mms_hog,Sms_hog^2,N,1);

MS_sag = wnormrnd(Mms_sag,Sms_sag^2,N,1);

    %gumbel distributed random value for added wave water hog/sag
B = Smw_hog*sqrt(6)/pi;
E = Mmw_hog+c*B;
MW_hog = wgumbrnd(B,E,0,N,1);

B = Smw_sag*sqrt(6)/pi;
E = Mmw_sag+c*B;
MW_sag = wgumbrnd(B,E,0,N,1);

    %gumbel distributed random value for slamming
B = Sms1*sqrt(6)/pi;
E = Mms1+c*B;
MS1 = wgumbrnd(B,E,0,N,1);

    %random number generation for normal distrebutef position of NA
    
   NA = wnormrnd(M_NA,S_NA^2,N,1);
   
 phi = mphi; %constant  
 
 
%mean and variance check
Mdiff_ST = abs(mean(ST)-MSt)
Mdiff_SC = abs(mean(SC)-MSc)
Mdiff_etaU = abs(mean(etaU)-Meta_u)
Mdiff_etaS = abs(mean(etaS)-Meta_s);
Mdiff_etaW = abs(mean(etaW)-Meta_w)
Mdiff_etaS1 = abs(mean(etaS1)-Meta_s1)
Mdiff_MS_sag = abs(mean(MS_sag)-Mms_sag)
Mdiff_MS_hog = abs(mean(MS_hog)-Mms_hog)
Mdiff_MW_sag = abs(mean(MW_sag)-Mmw_sag)
Mdiff_MW_hog = abs(mean(MW_hog)-Mmw_hog)
Mdiff_MS1 = abs(mean(MS1)-Mms1)


Sdiff_ST = std(ST)-SSt
Sdiff_SC = std(SC)-SSc
Sdiff_etaU = std(etaU)-Seta_u
Sdiff_etaS = std(etaS)-Seta_s
Sdiff_etaW = std(etaW)-Seta_w
Sdiff_etaS1 = std(etaS1)-Seta_s1
Sdiff_MS_sag = std(MS_sag)-Sms_sag
Sdiff_MS_hog = std(MS_hog)-Sms_hog
Sdiff_MW_sag = std(MW_sag)-Smw_sag
Sdiff_MW_hog = std(MW_hog)-Smw_hog
Sdiff_MS1 = std(MS1)-Sms1

% mu = [MSt ;MSc ;Meta_u ;Meta_s ;Meta_w ;Meta_s1 ;Mms_sag ;Mms_hog ;Mmw_sag ;Mmw_hog ;Mms1];
% var = [SSt ;SSc ;Seta_u ;Seta_s ;Seta_w ;Seta_s1 ;Sms_sag ;Sms_hog ;Smw_sag ;Smw_hog ;Sms1];
% 
% dis = [ST;SC;etaU;etaS;etaW;etaS1;MS_sag;MS_hog;MW_sag;MW_hog;MS1];
% 
%  
%  j=N;
% k=1;
% for i = length(mu)
%    
% mu_diff(i) = abs(mean(dis(k:j))-mu(i));
%S_diff(i)=std(k:j)-var(i); 
%         k=k+N;
%         j=j+N;
% end



%LOAD
MS_sag = abs(MS_sag);
MS_hog = abs(MS_hog);
MW_sag = abs(MW_sag);
MW_hog = abs(MW_hog);
%!Sagging!
%tension
MU = ST.*I./NA;
G_bot_sag = etaU.*MU - etaS.*MS_sag - phi.*etaW.*MW_sag;     %limitstatefunction for sagging bottom.
G_bot_sag =G_bot_sag';

%compression
MU = SC.*I./(D-NA);
G_top_sag =  etaU.*MU - etaS.*MS_sag - phi.*etaW.*MW_sag;   %limitstatefunction for sagging Top.
G_top_sag =G_top_sag';

%!Hoging!

%tension
MU = ST.*I./(D-NA);
G_top_hog = etaU.*MU - etaS.*MS_hog - phi.*etaW.*MW_hog;    %limitstatefunction for hoging top.
G_top_hog =G_top_hog';


%compression

MU = SC.*I./NA;
G_bot_hog = etaU.*MU - etaS.*MS_hog - phi.*etaW.*MW_hog; %limitstatefunction for hoging bottom.   
G_bot_hog =G_bot_hog';

%Slamming
%bottom
MU = SC.*I./NA;
G_sla_bot = etaU.*MU-etaS1.*MS1;%limitstatefunction for slamming bottom.
G_sla_bot =G_sla_bot';

%top
MU = SC.*I./(D-NA);
G_sla_top= etaU.*MU-etaS1.*MS1;%limitstatefunction for slamming top.
G_sla_top=G_sla_top';









%Probability check

%slamming

r=0;
for i=1:N
    G_sl(i) = min(G_sla_bot(i),G_sla_top(i));
    
    if G_sl(i) <=0 
        r=r+1;
        
    end 
    
end 
    P_sla = r/N;
    
%Hogging    
    r=0;
for i=1:N
    
   
    G_hog(i) = min(G_bot_hog(i),G_top_hog(i));
    
    if G_hog(i) <=0 
        r=r+1;
        
    end 
    
end 
    P_hog = r/N;
  
 %Sagging
    
    r=0;
for i=1:N 
    
   
    G_sag(i) = min(G_bot_sag(i),G_top_sag(i));    
    if G_sag(i) <=0 
        r=r+1;
        
    end 
    
end 

    P_sag = r/N;
    
    fprintf('The probability of failiure due to: \n Sagging = %.2g \n Hogging = %.2g \n Slamming = %.2g',P_sag,P_hog,P_sla);