%THIS PROGRAM DEMONSTRATES HODGKIN HUXLEY MODEL IN CURRENT CLAMP EXPERIMENTS AND SHOWS ACTION POTENTIAL PROPAGATION
%Time is in msec, voltage in mV, conductances in mho/mm^2, capacitance in microF/mm^2

%------------------------------------------------------------------------%
% Modification 1 : Running a loop to find the peaks for external current
% ranging from 0 to 0.7
%------------------------------------------------------------------------%
i=1;
delta = 0.001; 
for ImpCur=0:delta:0.7 
%ImpCur=input('enter the value of the impulse current in microamperes/mm^2: '); % 0.5
%TimeTot=input('enter the time for which stimulus is applied in milliseconds');

gkmax=.36;
vk=-77; 
gnamax=1.20;
vna=50; 
gl=0.003;
vl=-54.387; 
cm=.01; 
dt=0.01;

%------------------------------------------------------------------------%
niter=100000; % Modification 2: Changed iterations from 10000 to 100000
%------------------------------------------------------------------------%

t=(1:niter)*dt;
iapp=ImpCur*ones(1,niter);
%for i=1:100
 %   iapp(1,i)=ImpCur;
 %end;
v=-64.9964;
m=0.0530;
h=0.5960;
n=0.3177;

gnahist=zeros(1,niter);
gkhist=zeros(1,niter);
vhist=zeros(1,niter);
mhist=zeros(1,niter);
hhist=zeros(1,niter);
nhist=zeros(1,niter);


for iter=1:niter
  gna=gnamax*m^3*h; 
  gk=gkmax*n^4; 
  gtot=gna+gk+gl;
  vinf = ((gna*vna+gk*vk+gl*vl)+ iapp(iter))/gtot;
  tauv = cm/gtot;
  v=vinf+(v-vinf)*exp(-dt/tauv);
  alpham = 0.1*(v+40)/(1-exp(-(v+40)/10));
  betam = 4*exp(-0.0556*(v+65));
  alphan = 0.01*(v+55)/(1-exp(-(v+55)/10));
  betan = 0.125*exp(-(v+65)/80);
  alphah = 0.07*exp(-0.05*(v+65));
  betah = 1/(1+exp(-0.1*(v+35)));
  taum = 1/(alpham+betam);
  tauh = 1/(alphah+betah);
  taun = 1/(alphan+betan);
  minf = alpham*taum;
  hinf = alphah*tauh;
  ninf = alphan*taun;
  m=minf+(m-minf)*exp(-dt/taum);
  h=hinf+(h-hinf)*exp(-dt/tauh);
  n=ninf+(n-ninf)*exp(-dt/taun);
  vhist(iter)=v; mhist(iter)=m; hhist(iter)=h; nhist(iter)=n;
end

%------------------------------------------------------------------------%
% Modification 3 : Finding the total number of peaks
%------------------------------------------------------------------------%

j=1;
  total_peaks=zeros;      % Defining a variable for total number of peaks
  peaks=findpeaks(vhist); % Using findpeaks func to get number of peaks
for a=1:length(peaks)     % Running a loop to identify APs  
   if peaks(a) >=10       % min value for a waveform to considered an AP.
            total_peaks(j)=peaks(a);
            j=j+1;
   end;
end;
if total_peaks ~= 0       % Selecting peaks with non zero value
    num_of_peaks(i)=length(total_peaks);
else
    num_of_peaks(i)=0;
end;
i=i+1

end

%------------------------------------------------------------------------%

figure(1)
%subplot(2,1,1)
plot(t,vhist)
title('voltage vs time')

figure(2)
%subplot(2,1,2)
plot(t,mhist,'r-', t,hhist,'g.',t,nhist,'b-')
legend('m','h','n')

figure(3)
gna=gnamax*(mhist.^3).*hhist; 
  gk=gkmax*nhist.^4;
  clf
  plot(t,gna,'r');
  hold on
  plot(t,gk,'b');
  legend('gna','gk')
  hold off

%------------------------------------------------------------------------%
% Modification 4: Adding one more graph, which gives the plot of I_ext vs
% Firing rate
%------------------------------------------------------------------------%

figure(4);
x=0:delta:0.7;                                              % x coordinate
plot(x,num_of_peaks*1000/(niter/100));
xlabel('I_{Ext}');
ylabel('Firing Rate')
hold on;
for l=1:length(num_of_peaks)-1                      % Defining I1, I2 & I3
if num_of_peaks(l+1)>0 && num_of_peaks(l)==0       % One peak(AP) observed
         I1=(l)*delta  
end;
if num_of_peaks(l+1)>num_of_peaks(l)+5   % More than 5 peaks(APs) increase 
         I2=(l)*delta
end;
if num_of_peaks(l+1)<num_of_peaks(l)-5   % More than 5 peaks(APs) decrease
         I3=(l)*delta
end;
end;

%------------------------------------------------------------------------%
% Plotting I1, I2, I3
%------------------------------------------------------------------------%

I1=0*num_of_peaks*1000/(niter/100)+I1;          
plot(I1,num_of_peaks*1000/(niter/100));
text(I1(2),-3,'I1');
I2=0*num_of_peaks*1000/(niter/100)+I2;
plot(I2,num_of_peaks*1000/(niter/100));
text(I2(2),-3,'I2');
I3=0*num_of_peaks*1000/(niter/100)+I3;
plot(I3,num_of_peaks*1000/(niter/100));
text(I3(2),-3,'I3');

%------------------------------------------------------------------------%
% Ledger
%------------------------------------------------------------------------%

text(0.5,100,'0-I1 = No APs');
text(0.5,95,'I1-I2 = Finite no. of APs');
text(0.5,90,'I2-I3 = Infinite no.of APs');
text(0.5,85,'After I3 = No APs');


