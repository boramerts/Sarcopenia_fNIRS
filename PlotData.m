clear; clc;

reject = ["005", "008", "010", "011"];
sarco = ["009", "014", "015", "016", "017"];
n_subjects = 17;

% Define the number of strings
numStrings = n_subjects;

% Define the format string
formatStr = '%03d';

% Create the matrix of strings
strMatrix = cell(numStrings, 1);
for i = 1:numStrings
    strMatrix{i} = sprintf(formatStr, i);
end

grip_epochs = cell(2,1); % 2 = rest, 1 = block
nback_epochs = cell(2,1); % 2 = Nback, 1 = Oback
oddball_epochs = cell(2,1); % 2 = odd, 1 = std

experiments = {'Grip', 'Nback', 'Oddball'};
for i = 1:n_subjects
    for k = 1:3
        if ~(ismember(strMatrix{i,1},reject))
            fiff_file = '%s_%s_epochs.fif';
            fiff_file = sprintf(fiff_file, strMatrix{i,1}, experiments{1,k});

            events_file = '%s_%s_events.fif';
            events_file = sprintf(events_file, strMatrix{i,1}, experiments{1,k});

            if (ismember(strMatrix{i,1},sarco))
                status = "Sarcopeny";
            else
                status = "Control";
            end

            if k == 1
                try
                    cfg = [];
                    cfg.dataset = fiff_file;
                    grip_epochs{1,i} = ft_preprocessing(cfg);

                    event_file = mne_read_events(events_file);
                    grip_epochs{1,i}.trialinfo = event_file(:,3);
                    grip_epochs{1,i}.cfg.trl(:,4) = event_file(:,3);
                    grip_epochs{2,i} = status;
                catch
                    continue
                end
            elseif k == 2
                try
                    cfg = [];
                    cfg.dataset = fiff_file;
                    nback_epochs{1,i} = ft_preprocessing(cfg);

                    event_file = mne_read_events(events_file);
                    nback_epochs{1,i}.trialinfo = event_file(:,3);
                    nback_epochs{1,i}.cfg.trl(:,4) = event_file(:,3);
                    nback_epochs{2,i} = status;
                catch
                    continue
                end
            elseif k == 3
                try
                    cfg = [];
                    cfg.dataset = fiff_file;
                    oddball_epochs{1,i} = ft_preprocessing(cfg);

                    event_file = mne_read_events(events_file);
                    oddball_epochs{1,i}.trialinfo = event_file(:,3);
                    oddball_epochs{1,i}.cfg.trl(:,4) = event_file(:,3);
                    oddball_epochs{2,i} = status;
                catch
                    continue
                end
            end
        end
    end
end
channels = grip_epochs{1,1}.label;

% Grip--------------------------------

grip_block = cell(1,length(grip_epochs));
grip_rest = cell(1,length(grip_epochs));
grip_block_mean = cell(1,length(grip_epochs));
grip_rest_mean = cell(1,length(grip_epochs));

for i = 1:n_subjects
    if ~isempty(grip_epochs{1,i})
        for trials = 1:length(grip_epochs{1,i}.trialinfo)
            if grip_epochs{1,i}.trialinfo(trials) == 1
                grip_block{trials,i} = grip_epochs{1,i}.trial{1,trials};
            else
                grip_rest{trials,i} = grip_epochs{1,i}.trial{1,trials};
            end
        end
        for k = 1:length(grip_block(:,1))
            if ~isempty(grip_block{k,i})
                temp_block = grip_block(:,i);
                temp_block = temp_block(~cellfun('isempty',temp_block));
                if length(temp_block) == 2
                    grip_block_mean{1,i} = mean(cat(3,temp_block{:}),3);
                else
                    grip_block_mean{1,i} = mean(cat(length(temp_block),temp_block{:}),length(temp_block));
                end
            end
        end
        for k = 1:length(grip_rest(:,1))
            if ~isempty(grip_rest{k,i})
                temp_rest = grip_rest(:,i);
                temp_rest = temp_rest(~cellfun('isempty',temp_rest));
                if length(temp_rest) == 2
                    grip_rest_mean{1,i} = mean(cat(3,temp_rest{:}),3);
                else
                    grip_rest_mean{1,i} = mean(cat(length(temp_rest),temp_rest{:}),length(temp_rest));
                end
            end
        end
    end
end

ch_means_block_control = cell(3,length(channels));
ch_means_rest_control = cell(3,length(channels));
ch_means_block_sarco = cell(3,length(channels));
ch_means_rest_sarco = cell(3,length(channels));
for i = 1:length(channels)
    ch_means_block_control{1,i} = channels{i};
    ch_means_rest_control{1,i} = channels{i};
    ch_means_block_sarco{1,i} = channels{i};
    ch_means_rest_sarco{1,i} = channels{i};
    control_num_1 = 1;
    sarco_num_1 = 1;
    control_num_2 = 1;
    sarco_num_2 = 1;
    for subj = 1:n_subjects
        if ~isempty(grip_epochs{1,subj})
            if (ismember(channels{i,1},grip_epochs{1,subj}.label))
                if strcmp(grip_epochs{2,subj},"Sarcopeny")
                    ch_index  = find(cellfun(@(x)isequal(x,channels{i,1}),grip_epochs{1,subj}.label));
                    if ~isempty(grip_block_mean{1,subj})
                        ch_means_block_sarco{2,i}(sarco_num_1,:) = grip_block_mean{1,subj}(ch_index,:);
                        sarco_num_1 = sarco_num_1 + 1;
                    end
                    if ~isempty(grip_rest_mean{1,subj})
                        ch_means_rest_sarco{2,i}(sarco_num_2,:) = grip_rest_mean{1,subj}(ch_index,:);
                        sarco_num_2 = sarco_num_2 + 1;
                    end
                end
                if strcmp(grip_epochs{2,subj},"Control")
                    ch_index  = find(cellfun(@(x)isequal(x,channels{i,1}),grip_epochs{1,subj}.label));
                    if ~isempty(grip_block_mean{1,subj})
                        ch_means_block_control{2,i}(control_num_1,:) = grip_block_mean{1,subj}(ch_index,:);
                        control_num_1 = control_num_1 + 1;
                    end
                    if ~isempty(grip_rest_mean{1,subj})
                        ch_means_rest_control{2,i}(control_num_2,:) = grip_rest_mean{1,subj}(ch_index,:);
                        control_num_2 = control_num_2 + 1;
                    end
                end
            end
        end
    end
    ch_means_block_control{3,i} = mean(ch_means_block_control{2,i});
    ch_means_rest_control{3,i} = mean(ch_means_rest_control{2,i});
    ch_means_block_sarco{3,i} = mean(ch_means_block_sarco{2,i});
    ch_means_rest_sarco{3,i} = mean(ch_means_rest_sarco{2,i});
end

% Nback--------------------------------

nback = cell(1,length(nback_epochs));
oback = cell(1,length(nback_epochs));
nback_mean = cell(1,length(nback_epochs));
oback_mean = cell(1,length(nback_epochs));

for i = 1:n_subjects
    if ~isempty(nback_epochs{1,i})
        for trials = 1:length(nback_epochs{1,i}.trialinfo)
            if nback_epochs{1,i}.trialinfo(trials) == 1
                nback{trials,i} = nback_epochs{1,i}.trial{1,trials};
            else
                oback{trials,i} = nback_epochs{1,i}.trial{1,trials};
            end
        end
        for k = 1:length(nback(:,1))
            if ~isempty(nback{k,i})
                temp_nback = nback(:,i);
                temp_nback = temp_nback(~cellfun('isempty',temp_nback));
                if length(temp_nback) == 2
                    nback_mean{1,i} = mean(cat(3,temp_nback{:}),3);
                else
                    nback_mean{1,i} = mean(cat(length(temp_nback),temp_nback{:}),length(temp_nback));
                end
            end
        end
        for k = 1:length(oback(:,1))
            if ~isempty(oback{k,i})
                temp_oback = oback(:,i);
                temp_oback = temp_oback(~cellfun('isempty',temp_oback));
                if length(temp_oback) == 2
                    oback_mean{1,i} = mean(cat(3,temp_oback{:}),3);
                else
                    oback_mean{1,i} = mean(cat(length(temp_oback),temp_oback{:}),length(temp_oback));
                end
            end
        end
    end
end

ch_means_nback_sarco = cell(3,length(channels));
ch_means_oback_sarco = cell(3,length(channels));
ch_means_nback_control = cell(3,length(channels));
ch_means_oback_control = cell(3,length(channels));
for i = 1:length(channels)
    ch_means_nback_sarco{1,i} = channels{i};
    ch_means_oback_sarco{1,i} = channels{i};
    ch_means_nback_control{1,i} = channels{i};
    ch_means_oback_control{1,i} = channels{i};
    control_num_1 = 1;
    sarco_num_1 = 1;
    control_num_2 = 1;
    sarco_num_2 = 1;
    for subj = 1:n_subjects
        if ~isempty(nback_epochs{1,subj})
            if (ismember(channels{i,1},nback_epochs{1,subj}.label))
                if strcmp(nback_epochs{2,subj},"Sarcopeny")
                    ch_index  = find(cellfun(@(x)isequal(x,channels{i,1}),nback_epochs{1,subj}.label));
                    if ~isempty(nback_mean{1,subj})
                        ch_means_nback_sarco{2,i}(sarco_num_1,:) = nback_mean{1,subj}(ch_index,:);
                        sarco_num_1 = sarco_num_1 + 1;
                    end
                    if ~isempty(oback_mean{1,subj})
                        ch_means_oback_sarco{2,i}(sarco_num_2,:) = oback_mean{1,subj}(ch_index,:);
                        sarco_num_2 = sarco_num_2 + 1;
                    end
                end
                if strcmp(nback_epochs{2,subj},"Control")
                    ch_index  = find(cellfun(@(x)isequal(x,channels{i,1}),nback_epochs{1,subj}.label));
                    if ~isempty(nback_mean{1,subj})
                        ch_means_nback_control{2,i}(control_num_1,:) = nback_mean{1,subj}(ch_index,:);
                        control_num_1 = control_num_1 + 1;
                    end
                    if ~isempty(oback_mean{1,subj})
                        ch_means_oback_control{2,i}(control_num_2,:) = oback_mean{1,subj}(ch_index,:);
                        control_num_2 = control_num_2 + 1;
                    end
                end
            end
        end
    end
    ch_means_nback_sarco{3,i} = mean(ch_means_nback_sarco{2,i});
    ch_means_oback_sarco{3,i} = mean(ch_means_oback_sarco{2,i});
    ch_means_nback_control{3,i} = mean(ch_means_nback_control{2,i});
    ch_means_oback_control{3,i} = mean(ch_means_oback_control{2,i});
end

% Oddball--------------------------------

odd = cell(1,length(oddball_epochs));
std = cell(1,length(oddball_epochs));
odd_mean = cell(1,length(oddball_epochs));
std_mean = cell(1,length(oddball_epochs));

for i = 1:n_subjects
    if ~isempty(oddball_epochs{1,i})
        for trials = 1:length(oddball_epochs{1,i}.trialinfo)
            if oddball_epochs{1,i}.trialinfo(trials) == 1
                odd{trials,i} = oddball_epochs{1,i}.trial{1,trials};
            else
                std{trials,i} = oddball_epochs{1,i}.trial{1,trials};
            end
        end
        for k = 1:length(odd(:,1))
            if ~isempty(odd{k,i})
                temp_odd = odd(:,i);
                temp_odd = temp_odd(~cellfun('isempty',temp_odd));
                if length(temp_odd) == 2
                    odd_mean{1,i} = mean(cat(3,temp_odd{:}),3);
                else
                    odd_mean{1,i} = mean(cat(length(temp_odd),temp_odd{:}),length(temp_odd));
                end
            end
        end
        for k = 1:length(std(:,1))
            if ~isempty(std{k,i})
                temp_std = std(:,i);
                temp_std = temp_std(~cellfun('isempty',temp_std));
                if length(temp_std) == 2
                    std_mean{1,i} = mean(cat(3,temp_std{:}),3);
                else
                    std_mean{1,i} = mean(cat(length(temp_std),temp_std{:}),length(temp_std));
                end
            end
        end
    end
end

ch_means_odd_sarco = cell(3,length(channels));
ch_means_std_sarco = cell(3,length(channels));
ch_means_odd_control = cell(3,length(channels));
ch_means_std_control = cell(3,length(channels));
for i = 1:length(channels)
    ch_means_odd_sarco{1,i} = channels{i};
    ch_means_std_sarco{1,i} = channels{i};
    ch_means_odd_control{1,i} = channels{i};
    ch_means_std_control{1,i} = channels{i};
    control_num_1 = 1;
    sarco_num_1 = 1;
    control_num_2 = 1;
    sarco_num_2 = 1;
    for subj = 1:n_subjects
        if ~isempty(oddball_epochs{1,subj})
            if (ismember(channels{i,1},oddball_epochs{1,subj}.label))
                if strcmp(oddball_epochs{2,subj},"Sarcopeny")
                    ch_index  = find(cellfun(@(x)isequal(x,channels{i,1}),oddball_epochs{1,subj}.label));
                    if ~isempty(odd_mean{1,subj})
                        ch_means_odd_sarco{2,i}(sarco_num_1,:) = odd_mean{1,subj}(ch_index,:);
                        sarco_num_1 = sarco_num_1 + 1;
                    end
                    if ~isempty(std_mean{1,subj})
                        ch_means_std_sarco{2,i}(sarco_num_2,:) = std_mean{1,subj}(ch_index,:);
                        sarco_num_2 = sarco_num_2 + 1;
                    end
                end
                if strcmp(oddball_epochs{2,subj},"Control")
                    ch_index  = find(cellfun(@(x)isequal(x,channels{i,1}),oddball_epochs{1,subj}.label));
                    if ~isempty(odd_mean{1,subj})
                        ch_means_odd_control{2,i}(control_num_1,:) = odd_mean{1,subj}(ch_index,:);
                        control_num_1 = control_num_1 + 1;
                    end
                    if ~isempty(std_mean{1,subj})
                        ch_means_std_control{2,i}(control_num_2,:) = std_mean{1,subj}(ch_index,:);
                        control_num_2 = control_num_2 + 1;
                    end
                end
            end
        end
    end
    ch_means_odd_sarco{3,i} = mean(ch_means_odd_sarco{2,i});
    ch_means_std_sarco{3,i} = mean(ch_means_std_sarco{2,i});
    ch_means_odd_control{3,i} = mean(ch_means_odd_control{2,i});
    ch_means_std_control{3,i} = mean(ch_means_std_control{2,i});
end

clear i ch_index subj trials sub_num cfg event_file experiments events_file fiff_file formatStr i k numStrings strMatrix sarco_num1 sarco_num2 control_num_1 control_num_2 status

% Plot Graphs



set(0,'DefaultFigureVisible','off');

for ch = 1:2:length(channels)
%Grip
    figure
    t = linspace(-2,20,224);
    plot(t,ch_means_block_control{3,ch},'r'); hold on; 
    plot(t,ch_means_block_control{3,ch+1},'b'); hold on; 
    plot(t,ch_means_block_sarco{3,ch},'--r'); hold on;
    plot(t,ch_means_block_sarco{3,ch+1},'--b'); hold off;
    legend("Control HbO","Control Hb","Sarcopeny HbO","Sarcopeny Hb");
    titl = sprintf("Block %s", channels{ch}(1:5));
    title(titl);
    xlabel("Time (s)");
    ylabel("µM");
    xlim([-2.5 20]);
    f = gcf;
    name = sprintf('Graf/%s Block.png',channels{ch}(1:5));
    exportgraphics(f,name,'Resolution',300);

    figure
    t = linspace(-2,20,224);
    plot(t,ch_means_rest_control{3,ch},'r'); hold on; 
    plot(t,ch_means_rest_control{3,ch+1},'b'); hold on; 
    plot(t,ch_means_rest_sarco{3,ch},'--r'); hold on;
    plot(t,ch_means_rest_sarco{3,ch+1},'--b'); hold off;
    legend("Control HbO","Control Hb","Sarcopeny HbO","Sarcopeny Hb");
    titl = sprintf("Rest %s", channels{ch}(1:5));
    title(titl);
    xlabel("Time (s)");
    ylabel("µM");
    xlim([-2.5 25]);
    f = gcf;
    name = sprintf('Graf/%s Rest.png',channels{ch}(1:5));
    exportgraphics(f,name,'Resolution',300);
%Nback
    figure
    t = linspace(-2,27,296);
    plot(t,ch_means_oback_control{3,ch},'r'); hold on; 
    plot(t,ch_means_oback_control{3,ch+1},'b'); hold on; 
    plot(t,ch_means_oback_sarco{3,ch},'--r'); hold on;
    plot(t,ch_means_oback_sarco{3,ch+1},'--b'); hold off;
    legend("Control HbO","Control Hb","Sarcopeny HbO","Sarcopeny Hb");
    titl = sprintf("Oback %s", channels{ch}(1:5));
    title(titl);
    xlabel("Time (s)");
    ylabel("µM");
    xlim([-2.5 30]);
    f = gcf;
    name = sprintf('Graf/%s Oback.png',channels{ch}(1:5));
    exportgraphics(f,name,'Resolution',300);

    figure
    t = linspace(-2,27,296);
    plot(t,ch_means_nback_control{3,ch},'r'); hold on; 
    plot(t,ch_means_nback_control{3,ch+1},'b'); hold on; 
    plot(t,ch_means_nback_sarco{3,ch},'--r'); hold on;
    plot(t,ch_means_nback_sarco{3,ch+1},'--b'); hold off;
    legend("Control HbO","Control Hb","Sarcopeny HbO","Sarcopeny Hb");
    titl = sprintf("nback %s", channels{ch}(1:5));
    title(titl);
    xlabel("Time (s)");
    ylabel("µM");
    xlim([-2.5 30]);
    f = gcf;
    name = sprintf('Graf/%s nback.png',channels{ch}(1:5));
    exportgraphics(f,name,'Resolution',300);
%Oddball
    figure
    t = linspace(-2,37,397);
    plot(t,ch_means_odd_control{3,ch},'r'); hold on; 
    plot(t,ch_means_odd_control{3,ch+1},'b'); hold on; 
    plot(t,ch_means_odd_sarco{3,ch},'--r'); hold on;
    plot(t,ch_means_odd_sarco{3,ch+1},'--b'); hold off;
    legend("Control HbO","Control Hb","Sarcopeny HbO","Sarcopeny Hb");
    titl = sprintf("Oddball %s", channels{ch}(1:5));
    title(titl);
    xlabel("Time (s)");
    ylabel("µM");
    xlim([-2.5 40]);
    f = gcf;
    name = sprintf('Graf/%s Odd.png',channels{ch}(1:5));
    exportgraphics(f,name,'Resolution',300);

    figure
    t = linspace(-2,37,397);
    plot(t,ch_means_std_control{3,ch},'r'); hold on; 
    plot(t,ch_means_std_control{3,ch+1},'b'); hold on; 
    plot(t,ch_means_std_sarco{3,ch},'--r'); hold on;
    plot(t,ch_means_std_sarco{3,ch+1},'--b'); hold off;
    legend("Control HbO","Control Hb","Sarcopeny HbO","Sarcopeny Hb");
    titl = sprintf("Std %s", channels{ch}(1:5));
    title(titl);
    xlabel("Time (s)");
    ylabel("µM");
    xlim([-2.5 40]);
    f = gcf;
    name = sprintf('Graf/%s std.png',channels{ch}(1:5));
    exportgraphics(f,name,'Resolution',300);
end

close all;
set(0,'DefaultFigureVisible','on');













