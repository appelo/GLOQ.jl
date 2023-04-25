ramsey_data_12_0 = load('../data-set-20220413/population_ramsey_12_1000000.0_80000_20_1000_confusion_0.txt');
ramsey_data_12_1 = load('../data-set-20220413/population_ramsey_12_1000000.0_80000_20_1000_confusion_1.txt');
ramsey_data_12_2 = load('../data-set-20220413/population_ramsey_12_1000000.0_80000_20_1000_confusion_2.txt');
dark_time_12     = load('../data-set-20220413/darktime_ramsey_12_1000000.0_20_80000_1000_20220413_confusion.txt');

dark_time_12 = dark_time_12/1e9;
ramsey_data_12 = ramsey_data_12_1-mean(ramsey_data_12_1);



[freq_data_12,f_12] = ftfast(ramsey_data_12,dark_time_12);

data_12_decay = dark_time_12-mean(dark_time_12);
[freq_data_12_decay,f_12_decay] = ftfast(data_12_decay,dark_time_12);


figure(1)
plot(dark_time_12*1e9,ramsey_data_12_0,'r','Linewidth',1.5);
hold on
plot(dark_time_12*1e9,ramsey_data_12_1,'b','Linewidth',1.5);
plot(dark_time_12*1e9,ramsey_data_12_2,'k','Linewidth',1.5);
legend('0','1','2');
xlabel('ns')
ylabel('Ramsey12')

figure(100000)
plot(dark_time_12*1e9,ramsey_data_12,'k','Linewidth',1.5);

figure(4)
plot(f_12,abs(freq_data_12),'-o','Linewidth',1.5);
xlim([0,2]*1e6)
grid
title('Ramsey 1-2')
xlabel('Hz');
%xticks(0e6+0:0.05e6:2e6)
%xticks(0e5+0:0.05e5:2e5)

figure(50)
plot(f_12,abs(freq_data_12_decay),'-o','Linewidth',1.5);
grid
title('Ramsey 1-2')
xlabel('Hz');
xlim([0,2]*1e6);
xticks(0e6+0:0.025e6:2e6)
set(gca,'FontSize',15)

figure(5)
plot(f_12,abs(freq_data_12),'-o','Linewidth',1.5);
%xlim([0,0.3]*1e6)
grid
title('Ramsey 1-2')
xlabel('Hz');
xlim([0,2]*1e6);
xticks(0e6+0:0.025e6:2e6)
set(gca,'FontSize',15)

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




