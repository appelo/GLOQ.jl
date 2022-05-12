using FFTW
using DelimitedFiles
using Plots
using LinearAlgebra
using Statistics
pyplot()
# Read data
ramsey_data_01_0 = readdlm("../data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_0.txt");
ramsey_data_01_1 = readdlm("../data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_0.txt");
ramsey_data_01_2 = readdlm("../data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_0.txt");
dark_time_01     = readdlm("../data-set-20220413/darktime_ramsey_01_half_period_1000000.0_1000_20220413_80000_20.txt");

ramsey_data_12_0 = readdlm("../data-set-20220413/population_ramsey_12_1000000.0_80000_20_1000_confusion_0.txt");
ramsey_data_12_1 = readdlm("../data-set-20220413/population_ramsey_12_1000000.0_80000_20_1000_confusion_1.txt");
ramsey_data_12_2 = readdlm("../data-set-20220413/population_ramsey_12_1000000.0_80000_20_1000_confusion_2.txt");
dark_time_12     = readdlm("../data-set-20220413/darktime_ramsey_12_1000000.0_20_80000_1000_20220413_confusion.txt");

dark_time_01 = dark_time_01/1e9;
ramsey_data_01 = ramsey_data_01_1.-mean(ramsey_data_01_1);

dark_time_12 = dark_time_12/1e9;
ramsey_data_12 = ramsey_data_12_2.-mean(ramsey_data_12_2);

function ftfast(y,t)
	T=t[2]-t[1];
	Y=T*fft(y);
	ny=length(y);
	f=collect(0:ny-1)./ny;
	Y=fftshift(Y);
	f=fftshift(f);
	f[f.>=0.5]=f[f.>=0.5].-1;
	f=1/T*f;
    return Y[:],f[:]
end

@time freq_data_01,f_01 = ftfast(ramsey_data_01,dark_time_01);
@time freq_data_12,f_12 = ftfast(ramsey_data_12,dark_time_12);


fig_fft_01 = plot(f_01,abs.(freq_data_01),
                  line=(2.5,:dash),marker=(8,:circle),
                  xlim=(0.5e6,1.5e6),
                  size=(1200,800),
                  legend=:false,
                  xtickfontsize=18,ytickfontsize=18,titlefontsize=24,
                  title = "Ramsey 0-1 with nominal 1MHz detuning\n Magnitude of population for state 1"
                  );
xticks!(fig_fft_01,
        [0.6e6,1.01e6,1.4e6],
        ["0.6MHz","1.01MHz","1.4MHz"])
plot!(fig_fft_01,[1.01e6;1.01e6],[0;3.5e-6],
    line=(:dash,2.0),c=:red)
display(fig_fft_01)


fig_fft_12 = plot(f_12,abs.(freq_data_12),
                  line=(2.5,:dash),marker=(8,:circle),
                  xlim=(0.5e6,1.5e6),
                  size=(1200,800),
                  legend=:false,
                  xtickfontsize=18,ytickfontsize=18,titlefontsize=24,
                  title = "Ramsey 1-2 with nominal 1MHz detuning\n Magnitude of population for state 2"
                  );
xticks!(fig_fft_12,
        [0.6e6,1.03e6-0.131e6,1.03e6+0.131e6,1.4e6],
        ["0.6MHz","1.03MHz-131KHz","1.03MHz+131KHz","1.4MHz"])
plot!(fig_fft_12,[1.03e6-0.131e6;1.03e6-0.131e6],[0;2e-6],
      line=(:dash,2.0),c=:red)
plot!(fig_fft_12,[1.03e6+0.131e6;1.03e6+0.131e6],[0;2.1e-6],
      line=(:dash,2.0),c=:red)
display(fig_fft_12)

savefig(fig_fft_01,"../picture/fft_01.png")
savefig(fig_fft_12,"../picture/fft_12.png")

#=
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



=#



