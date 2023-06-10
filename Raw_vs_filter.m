clear; clc;

restoredefaultpath
addpath /Users/boramert/Documents/MATLAB/fieldtrip-20221223
addpath /Users/boramert/Documents/MATLAB/fieldtrip-20221223/external/artinis
addpath /Users/boramert/Documents/MATLAB/fieldtrip-20221223/external/mne
ft_defaults

sarco = ["006", "009", "014", "016", "019", "020", "023", "025", "032"];
experiment = {'Grip', 'Nback', 'Oddball'};

nirs_data = cell(1,6);

num = 1;
for subj = 1:33
    for exp = 1:3
        sname = sprintf('%03d', subj);
        ex = experiment{exp};
        eve_file = sprintf('/Users/boramert/Desktop/Yüksek Lisans/Matlab_Kod/Raw_vs_Filter/fif/%s_%s_eve.fif', sname, ex);
        raw_file = sprintf('/Users/boramert/Desktop/Yüksek Lisans/Matlab_Kod/Raw_vs_Filter/fif/%s_%s_raw_raw.fif', sname, ex);
        filt_file = sprintf('/Users/boramert/Desktop/Yüksek Lisans/Matlab_Kod/Raw_vs_Filter/fif/%s_%s_filt_raw.fif', sname, ex);

        if(isfile(raw_file))
            nirs_data{num, 1} = sprintf("%s %s",sname,ex);
            [event, ~] = fiff_read_events(eve_file);
            nirs_data{num, 4} = event;

            cfg         = [];
            cfg.dataset = raw_file;
            data_mp     = ft_preprocessing(cfg);
            ft_datatype(data_mp); % should return 'raw'

            % Make fs exactly 10 Hz
            cfg                   = [];
            cfg.resamplefs        = 10;
            data_mp             = ft_resampledata(cfg, data_mp);

            nirs_data{num,2} = [{data_mp.label}, data_mp.trial(:)'];

            cfg         = [];
            cfg.dataset = filt_file;
            data_mp     = ft_preprocessing(cfg);
            ft_datatype(data_mp); % should return 'raw'

            % Make fs exactly 10 Hz
            cfg                   = [];
            cfg.resamplefs        = 10;
            data_mp             = ft_resampledata(cfg, data_mp);

            nirs_data{num,3} = [{data_mp.label}, data_mp.trial(:)'];

            nirs_data{num,6} = 0 + 1*(ismember(sname, sarco));
            num = num+1;
        end
    end
end

for i=1:length(nirs_data)
    boxcar = zeros(1,length(nirs_data{i,2}{1,2}));
    if contains(nirs_data{i,1}, 'Grip')
        t = 35*10;
    elseif contains(nirs_data{i,1}, 'Nback')
        t = 37.5*10;
    elseif contains(nirs_data{i,1}, 'Oddball')
        t = 47.4*10;
    end

    for k=1:length(nirs_data{i,4})
        if nirs_data{i,4}(k,3) == 1
            boxcar(nirs_data{i,4}(k,1):nirs_data{i,4}(k,1)+t) = 1/950;
        end
    end
    nirs_data{i,5} = boxcar;
end

clearvars -except nirs_data;
set(0,'DefaultFigureVisible','off');

for type = 2:3
    for i=1:height(nirs_data)
        fprintf("Processiing Figure No: %d...\n", i);
        figure;
        t = linspace(0,length(nirs_data{i,2}{1,2})/10,length(nirs_data{i,2}{1,2}));
        num = 1;
        p = area(nirs_data{i,5}); hold on;
        p.FaceAlpha = 0.2;
        if type == 2
            titl = sprintf("%s %s",nirs_data{i,1}, "Raw");
        else
            titl = sprintf("%s %s",nirs_data{i,1}, "Filtered");
        end
        for ch=1:2:height(nirs_data{i,2}{1,2})-1
            plot(nirs_data{i,type}{1,2}(ch,:)+(0.00005*num),'r'); hold on;
            plot(nirs_data{i,type}{1,2}(ch+1,:)+(0.00005*num),'b'); hold on;
            text(length(nirs_data{i,type}{1,2}(ch+1,:)),mean(nirs_data{i,3}{1,2}(ch+1,:)+(0.00005*num)),nirs_data{i,type}{1,1}{ch}(1:5),'Interpreter','none');
            xlim tight;
            ylim([0 0.001052631578947]);
            title(titl);
            num = num +1;
        end
        disp("Plotting complete...\n");
        f = gcf;
        name = sprintf('/Users/boramert/Desktop/Yüksek Lisans/Matlab_Kod/Raw_vs_Filter/plt/%s %s.png',titl,nirs_data{i,1}{1}(5:end));
        exportgraphics(f,name,'Resolution',300);
        disp("Plot saved...\n");
    end
end

close all;
set(0,'DefaultFigureVisible','on');

% TODO: Add channel names - Add EMG & Dyno
