clear; clc;

n_subjects = 40;
fs = 1000;

EMG = cell(n_subjects+1,5);
[EMG{1,:}] = deal("Subject","Raw EMG","Dynamo","Rect. EMG","Smooth EMG");

%set(0,'DefaultFigureVisible','off');

for i = 1:n_subjects
    if i<10
        subject = sprintf("00%d",i);
    else
        subject = sprintf("0%d",i);
    end

    filename = sprintf("/Users/boramert/Desktop/Yüksek Lisans/Data_Results/Data/EMG_Data/EMG/EMG_%s.mat",subject);
    if exist(filename,'file')
        data = load(filename);

        emg = data.data(:,1);
        dyno = data.data(:,2);
        clear data;

        t = linspace(0,length(emg/1000),length(emg));

        [b,a] = butter(4,[50,300]/500,'bandpass');
        emg_filt = filtfilt(b,a,emg);
        emg_rect = abs(emg);
        emg_filt_rect = abs(emg_filt);

        sig = sqrt(movmean(emg_filt_rect.^2, 250));
        smooth = movmean(sig, 10000);

        [EMG{i+1,:}] = deal(subject,emg,dyno,emg_filt_rect,smooth);

%         f = figure;
%         fig = tiledlayout(2,1);
% 
%         nexttile;
%         plot(t,emg);
%         title("Raw Data");
%         xlabel("Time (ms)");
%         ylabel("Amplitude");
%         ylim([-2.5 2.5]);
%         nexttile;
%         plot(sig);
%         hold on;
%         plot(smooth);
%         hold on;
%         plot(t,(dyno-mean(dyno))/100);
%         hold off;
%         title("Rectified + Smooth + Dynamo");
%         xlabel("Time (ms)");
%         ylabel("Amplitude");
%         ylim([-0.25 1.25]);
%         legend("Rectified Signal","Smooth Signal","Dynamo Signal");
% 
%         title(fig, subject);
% 
%         f = gcf;
%         name = sprintf('figs/EMG_%s.png',subject);
%         exportgraphics(f,name,'Resolution',300);
    end
end

EMG(all(cellfun(@isempty, EMG),2),:) = [];
clearvars -except EMG;
%close all;
%set(0,'DefaultFigureVisible','on');
