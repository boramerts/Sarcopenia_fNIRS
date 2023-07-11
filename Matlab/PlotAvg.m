% clear; clc;
% 
% restoredefaultpath
% addpath /Users/boramert/Documents/MATLAB/fieldtrip-20221223
% addpath /Users/boramert/Documents/MATLAB/fieldtrip-20221223/external/artinis
% addpath /Users/boramert/Documents/MATLAB/fieldtrip-20221223/external/mne
% ft_defaults
% 
% % List of rejected subjects & Sarco subjects
% reject = ["003","004","005", "008", "010", "011","032","034"];
% reject_oddball = "019";
% reject_nback = "018";
% reject_grip = "018";
% sarco = ["006","009", "014", "016", "019","020","023","025","037","038","039","040"];
% 
% n_subjects = 42; % Update with total number of subjects
% 
% % Define the number of strings
% numStrings = n_subjects;
% 
% % Define the format string
% formatStr = '%03d';
% 
% % Create the matrix of strings
% strMatrix = cell(numStrings, 1);
% for i = 1:numStrings
%     strMatrix{i} = sprintf(formatStr, i);
% end
% 
% grip_epochs = cell(3,1); % 2 = rest, 1 = block
% nback_epochs = cell(3,1); % 2 = Oback, 1 = Nback
% oddball_epochs = cell(3,1); % 2 = Std, 1 = Odd
% % Row 1 = Raw data
% % Row 2 = Epoch data
% % Row 3 = Group
% % Row 4 = Subject no
% 
% % Import & clear data
% experiments = {'Grip', 'Nback', 'Oddball'};
% for i = 1:n_subjects
%     for k = 1:3
%         fprintf("Subject: %s ************************* \n", strMatrix{i,1})
%         if ~(ismember(strMatrix{i,1},reject))
%             fiff_file = '/Users/boramert/Desktop/Yüksek Lisans/Exports/fNIRS_Data/%s_%s_epochs.fif'; % Path of epochs data
%             fiff_file = sprintf(fiff_file, strMatrix{i,1}, experiments{1,k});
%             data_file = '/Users/boramert/Desktop/Yüksek Lisans/Exports/fNIRS_Data/%s_%s_data.fif'; % path for raw data
%             data_file = sprintf(data_file, strMatrix{i,1}, experiments{1,k});
%             event_file = '/Users/boramert/Desktop/Yüksek Lisans/Exports/fNIRS_Data/%s_%s_events.fif'; % path for events
%             event_file = sprintf(event_file, strMatrix{i,1}, experiments{1,k});
%             cfg_events = [];
%             cfg_events.dataset = fiff_file;
%             cfg_data = [];
%             cfg_data.dataset = data_file;
% 
%             if (ismember(strMatrix{i,1},sarco)) % Add subject group info
%                 status = "Sarcopeny";
%             else
%                 status = "Control";
%             end
% 
%             if k == 1
%                 disp("Found Grip ***************")
%                 try
% 
%                     grip_epochs{1,i} = ft_preprocessing(cfg_data);
%                     grip_epochs{2,i} = ft_preprocessing(cfg_events);
% 
%                     event_file = mne_read_events(event_file);
%                     grip_epochs{2,i}.trialinfo = event_file;
%                     grip_epochs{3,i} = status;
%                     grip_epochs{4,i} = strMatrix{i,1};
%                     disp("=============================================")
%                     fprintf("Added Subject: %s Experiment: Grip", strMatrix{i,1})
%                     disp("=============================================")
%                 catch
%                     disp("Error \n")
%                     continue
%                 end
%             elseif k == 2
%                 disp("Found Nback ***************")
%                 try
%                     nback_epochs{1,i} = ft_preprocessing(cfg_data);
%                     nback_epochs{2,i} = ft_preprocessing(cfg_events);
% 
%                     event_file = mne_read_events(event_file);
%                     nback_epochs{2,i}.trialinfo = event_file;
%                     nback_epochs{3,i} = status;
%                     nback_epochs{4,i} = strMatrix{i,1};
%                     disp("=============================================")
%                     fprintf("Added Subject: %s Experiment: Nback", strMatrix{i,1})
%                     disp("=============================================")
%                 catch
%                     disp("Error \n")
%                     continue
%                 end
%             elseif k == 3
%                 disp("Found Oddball ***************")
%                 try
%                     oddball_epochs{1,i} = ft_preprocessing(cfg_data);
%                     oddball_epochs{2,i} = ft_preprocessing(cfg_events);
% 
%                     event_file = mne_read_events(event_file);
%                     oddball_epochs{2,i}.trialinfo = event_file;
%                     oddball_epochs{3,i} = status;
%                     oddball_epochs{4,i} = strMatrix{i,1};
%                     disp("=============================================")
%                     fprintf("Added Subject: %s Experiment: Oddball \n", strMatrix{i,1})
%                     disp("=============================================")
%                 catch
%                     disp("Error \n")
%                     continue
%                 end
%             end
%         end
%     end
% end
% 
% % Remove empty columns (rejects)
% grip_epochs(:, any(cellfun(@isempty, grip_epochs), 1)) = [];
% nback_epochs(:, any(cellfun(@isempty, nback_epochs), 1)) = [];
% oddball_epochs(:, any(cellfun(@isempty, oddball_epochs), 1)) = [];
% 
% % Get the list of channels from first subject
% channels = grip_epochs{1,1}.label;
% 
% clc;
% 
% data = {grip_epochs; nback_epochs; oddball_epochs};
% 
% clearvars -except data channels sarco n_subjects

[ch_means_control,ch_means_sarco,times] = getEpochsCell(data, channels);
ch_means_control = getEpochAvgs(ch_means_control);
ch_means_sarco = getEpochAvgs(ch_means_sarco);
plotAvgs(ch_means_control, ch_means_sarco, channels, times);



function plotAvgs(ch_means_control, ch_means_sarco, channels, times)
set(0,'DefaultFigureVisible','off');
path = "/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/fNIRS_analysis/Power/Epochs/";
trials = ["Rest","Grip","0-Back","2-Back","Standard","Oddball"];
for col = 1:6
    for row = 1:2:length(channels)
        figure;
        plot(times{col},ch_means_control{row,col},'r'); hold on; % HbO of control group (First Trial)
        plot(times{col},ch_means_control{row+1,col},'b'); hold on; % HbR of control group (Second Trial)
        plot(times{col},ch_means_sarco{row,col},'--r'); hold on; % HbO of sarco group (First Trial)
        plot(times{col},ch_means_sarco{row+1,col},'--b'); hold on; % HbR of sarco group (Second Trial)
        legend("Control HbO","Control HbR","Sarcopeny HbO","Sarcopeny HbR");
        title(sprintf("%s %s", trials(col),channels{row}(1:end-4)),'Interpreter', 'none');
        xlabel("Time (s)");
        ylabel("µM");
        xlim([times{col}(1) times{col}(end)]);
        f = gcf;
        filename = sprintf("%s%s_%s.png",path, trials(col), channels{row}(1:end-4));
        exportgraphics(f,filename,"Resolution",300);
        close all;
    end
end
set(0,'DefaultFigureVisible','on');
end

function ch_means = getEpochAvgs(ch_means) % Get the mean of all trials for all subjects for each channel
for row = 1:size(ch_means,1)
    for col = 1:size(ch_means,2)
        ch_means{row,col} = mean(ch_means{row,col});
    end
end
end

function [times, ch_means] = addToCell(ch_index, ch, subj, exp, trial, data, ch_means, times) % Extra function to make getEpochsCell shorter
if exp == 1
    if data{exp}{2,subj}.trialinfo(trial,3) == 2
        %fprintf("ch: %d - trial: %d - exp: %d - subj: %d \n", ch, trial, exp, subj);
        ch_means{ch,1}(end+1,:) = data{exp}{2,subj}.trial{trial}(ch_index,:);
        times{1} = data{exp}{2,subj}.time{trial};
    else
        %fprintf("ch: %d - trial: %d - exp: %d - subj: %d \n", ch, trial, exp, subj);
        ch_means{ch,2}(end+1,:) = data{exp}{2,subj}.trial{trial}(ch_index,:);
        times{2} = data{exp}{2,subj}.time{trial};
    end
elseif exp == 2
    if data{exp}{2,subj}.trialinfo(trial,3) == 2
        %fprintf("ch: %d - trial: %d - exp: %d - subj: %d \n", ch, trial, exp, subj);
        ch_means{ch,3}(end+1,:) = data{exp}{2,subj}.trial{trial}(ch_index,:);
        times{3} = data{exp}{2,subj}.time{trial};
    else
        %fprintf("ch: %d - trial: %d - exp: %d - subj: %d \n", ch, trial, exp, subj);
        ch_means{ch,4}(end+1,:) = data{exp}{2,subj}.trial{trial}(ch_index,:);
        times{4} = data{exp}{2,subj}.time{trial};
    end
else
    if data{exp}{2,subj}.trialinfo(trial,3) == 2
        %fprintf("ch: %d - trial: %d - exp: %d - subj: %d \n", ch, trial, exp, subj);
        ch_means{ch,5}(end+1,:) = data{exp}{2,subj}.trial{trial}(ch_index,:);
        times{5} = data{exp}{2,subj}.time{trial};
    else
        %fprintf("ch: %d - trial: %d - exp: %d - subj: %d \n", ch, trial, exp, subj);
        ch_means{ch,6}(end+1,:) = data{exp}{2,subj}.trial{trial}(ch_index,:);
        times{6} = data{exp}{2,subj}.time{trial};
    end
end
end

function [ch_means_control,ch_means_sarco,times] = getEpochsCell(data, channels) % Get a cell that contains all epochs of a channel for a trial
ch_means_control = cell(40,6); % 1: Rest 2: Block 3: 0-Back 4: 2-Back 5: Std 6: Oddball (For each channel)
ch_means_sarco = cell(40,6); % 1: Rest 2: Block 3: 0-Back 4: 2-Back 5: Std 6: Oddball (For each channel)

times = cell(1,6); % Get time matrices of trials

for exp = 1:3
    for ch = 1:length(channels)
        for subj = 1:length(data{exp})
            ch_index = getChannelIndex(channels{ch},data{exp}{2,subj}.label);
            if ch_index > 0 % Returns 0 if subject does not have the channel
                for trial = 1:size(data{exp}{2,subj}.trialinfo,1) % Run throught trials of subject
                    if strcmp(data{exp}{3,subj},"Control")
                        [times, ch_means_control] = addToCell(ch_index, ch, subj, exp, trial, data, ch_means_control, times);
                    else
                        [times, ch_means_sarco] = addToCell(ch_index, ch, subj, exp, trial, data, ch_means_sarco, times);
                    end
                end
            end
        end
    end
end
end

function index = getChannelIndex(channel, subj_channels) % Get index of channel in the subject struct
for ch = 1:length(subj_channels)
    if any(strcmp(channel,subj_channels{ch}))
        index = ch;
        break
    else
        index = 0;
    end
end
end









