% [Y,f]=ftfast(y,t)
%
%	Calculates the approximated continuous
%       fourier transform for signal y.
%	t is the time vector corresponding to y and should
%	be equally spaced.
%       ( utilizes the  function fft.m )
%		

close all;

%t = linspace(0,1,201);
%omega1 = 10.0;
%omega2 = 9.0;
%time_data = sin(omega1*2*pi*t)+sin(omega2*2*pi*t);
%
%[freq_data,f] = ftfast(time_data,t); 
%
%figure(10)
%plot(t,time_data);
%
%figure(20)
%plot(f,freq_data);
%xlim([0,20]);
%grid
%
%figure(200)
%plot(f,freq_data);
%xlim([8,12]);
%grid

%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0316/20220316/population_ramsey_01_250000.0_1000_20220316_40000_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0316/20220316/darktime_ramsey_01_250000.0_1000_20220316__40000.txt');

%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0316/20220316/population_ramsey_01_250000.0_500_20220316_40000_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0316/20220316/darktime_ramsey_01_250000.0_500_20220316__40000.txt');

%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0316/20220316/population_ramsey_01_250000.0_1000_20220316_160000_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0316/20220316/darktime_ramsey_01_250000.0_1000_20220316__160000.txt');

%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0318/20220318/population_ramsey_01_250000.0_500_20220318_160000_152_80_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0318/20220318/darktime_ramsey_01_250000.0_500_20220318_160000_152_80.txt');

%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0318/20220318/population_ramsey_01_250000.0_1000_20220318_160000_80_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0318/20220318/darktime_ramsey_01_250000.0_1000_20220318_160000_80.txt');
%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0322/20220322/population_rabi_01_1000_6000_0_20220322_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0322/20220322/time_rabi_01_1000_6000_0_20220322.txt');
%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0322/20220322/population_ramsey_01_1000000.0_1000_20220322_1000_5_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0322/20220322/darktime_ramsey_01_1000000.0_1000_20220322_1000_5.txt');
%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323/population_ramsey_01_1000000.0_500_20220323_20000_5_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323/darktime_ramsey_01_1000000.0_500_20220323_20000_5.txt');
%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323/population_ramsey_01_1000000.0_500_20220323_2000_5_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323/darktime_ramsey_01_1000000.0_500_20220323_2000_5.txt');
%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323/population_ramsey_01_2000000.0_500_20220323_2000_5_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323/darktime_ramsey_01_2000000.0_500_20220323_2000_5.txt');
%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323/population_ramsey_01_1000000.0_500_20220323_2000_5_0.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323/darktime_ramsey_01_half_period_1000000.0_2000_20220324_2000_5.txt');

%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323/population_ramsey_01_1000000.0_500_20220323_20000_5_0.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323/darktime_ramsey_01_half_period_1000000.0_2000_20220324_20000_5.txt');
%ramsey_data_01 = ramsey_data_01(1:501);
%dark_time_01   = dark_time_01(1:501);

%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323/population_ramsey_01_1000000.0_500_20220323_20000_5_0.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323/darktime_ramsey_01_half_period_1000000.0_2000_20220324_20000_5.txt');
%ramsey_data_01 = ramsey_data_01(1:501);
%dark_time_01   = dark_time_01(1:501);

%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323YanivSuggestion/population_ramsey_01_250000.0_200_20220324_4000_100_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0323/20220323YanivSuggestion/darktime_ramsey_01_250000.0_200_20220324_4000_100.txt');

%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0411/population_ramsey_01_half_period_1000000.0_1000_20220412_80000_20_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0411/darktime_ramsey_01_half_period_1000000.0_1000_20220412_80000_20.txt');

ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_1.txt');
dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/darktime_ramsey_01_half_period_1000000.0_1000_20220413_80000_20.txt');

%ramsey_data_01 = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_01_half_period_1000000.0_1000_20220413_160000_40_1.txt');
%dark_time_01   = load('/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/darktime_ramsey_01_half_period_1000000.0_1000_20220413_160000_40.txt');

dark_time_01 = dark_time_01/1e9;
ramsey_data_01 = ramsey_data_01-mean(ramsey_data_01);

length(ramsey_data_01)

[freq_data_01,f_01] = ftfast(ramsey_data_01,dark_time_01);

%ramsey_data_12 = load('../0309/data0309/population_ramsey_12_250000.0_1000_2.txt');
%dark_time_12   = load('../0309/data0309/darktime_ramsey_12_250000.0_1000_20220309.txt');
%dark_time_12 = dark_time_12/1e9;
%ramsey_data_12 = ramsey_data_12-mean(ramsey_data_12);
%dark_time_12 = dark_time_12/1000;

%[freq_data_12,f_12] = ftfast(ramsey_data_12,dark_time_12);


fake_ramsey_data_01 = 0.5*exp(-dark_time_01/10).*cos(2*pi*0.25*dark_time_01)+0.5;
[fake_freq_data_01,f_01] = ftfast(fake_ramsey_data_01,dark_time_01);


figure(1)
plot(dark_time_01*1e9,ramsey_data_01,'-o','Linewidth',1.5);
xlabel('ns')
ylabel('Ramsey01')
%plot(f_01,freq_data_01);
%xlim([0,3]*1e6)
%title('Ramsey 0-1')
%xlabel('Hz');

figure(2)
plot(f_01,freq_data_01,'-o','Linewidth',1.5);
xlim([0,25]*1e5)
grid
title('Ramsey 0-1')
xlabel('Hz');

figure(3)
plot(f_01,freq_data_01,'-o','Linewidth',1.5);
xlim([0,0.3]*1e6)
grid
title('Ramsey 0-1')
xlabel('Hz');
xticks(0:0.1e5:4e5)

figure(4)
plot(f_01,abs(freq_data_01),'-o','Linewidth',1.5);
xlim([0,2]*1e6)
grid
title('Ramsey 0-1')
xlabel('Hz');
%xticks(0e6+0:0.05e6:2e6)
%xticks(0e5+0:0.05e5:2e5)

figure(5)
plot(f_01,abs(freq_data_01),'-o','Linewidth',1.5);
%xlim([0,0.3]*1e6)
grid
title('Ramsey 0-1')
xlabel('Hz');
xlim([0,2]*1e6);
xticks(0e6+0:0.1e6:2e6)

%saveas(3,"fft_01_80.png");
%saveas(5,"fft_01_160_20220318.png");
%saveas(5,"fft_01_160_20220318_76.png");
%saveas(5,"fft_01_160_20220318_152.png");

%figure(3)
%plot(f_12,freq_data_12);
%xlim([0,10]*1e6)
%title('Ramsey 1-2')
%xlabel('Hz');
%
%figure(4)
%plot(f_12,freq_data_12,'-o','Linewidth',1.5);
%xlim([0,1]*1e6)
%grid
%title('Ramsey 1-2')
%xlabel('Hz');


function [Y,f]=ftfast(y,t)
	T=t(2)-t(1);
	Y=T*fft(y);
	ny=length(y);
	f=[(0:ny-1)./ny];
	Y=fftshift(Y);
	f=fftshift(f);
	f_index=find(f>=0.5);
	f(f_index)=f(f_index)-1;
	f=1/T*f;
end




