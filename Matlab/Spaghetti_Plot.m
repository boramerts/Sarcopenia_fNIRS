clear; clc;

restoredefaultpath
addpath /Users/boramert/Documents/MATLAB/fieldtrip-20221223
addpath /Users/boramert/Documents/MATLAB/fieldtrip-20221223/external/artinis
addpath /Users/boramert/Documents/MATLAB/fieldtrip-20221223/external/mne
ft_defaults

reject = ["003","004","005", "008", "010", "011","032","034"];
% reject = ["001","003","004","005","007", "008", "010", "011","012","018","019","021","026","028","029","032","034","033","035","036","038","039"];
reject_oddball = "019";
reject_nback = "018";
reject_grip = "018";
sarco = ["006","009", "014", "016", "019","020","023","025","037","038","039","040"];

n_subjects = 40;

% Define the number of strings
numStrings = n_subjects;

% Define the format string
formatStr = '%03d';

% Create the matrix of strings
strMatrix = cell(numStrings, 1);
for i = 1:numStrings
    strMatrix{i} = sprintf(formatStr, i);
end

grip_epochs = cell(3,1); % 2 = rest, 1 = block
nback_epochs = cell(3,1); % 2 = Oback, 1 = Nback
oddball_epochs = cell(3,1); % 2 = Std, 1 = Odd

% Import & clear data
experiments = {'Grip', 'Nback', 'Oddball'};
for i = 1:n_subjects
    for k = 1:3
        fprintf("Subject: %s ************************* \n", strMatrix{i,1})
        if ~(ismember(strMatrix{i,1},reject))
            fiff_file = '/Users/boramert/Desktop/Yüksek Lisans/Exports/fNIRS_Data/%s_%s_epochs.fif';
            fiff_file = sprintf(fiff_file, strMatrix{i,1}, experiments{1,k});
            cfg = [];
            cfg.dataset = fiff_file;

            if (ismember(strMatrix{i,1},sarco))
                status = "Sarcopeny";
            else
                status = "Control";
            end

            if k == 1
                disp("Found Grip ***************")
                try

                    grip_epochs{1,i} = ft_preprocessing(cfg);

                    event_file = mne_read_events(fiff_file);
                    grip_epochs{1,i}.trialinfo = event_file;
                    grip_epochs{2,i} = status;
                    grip_epochs{3,i} = strMatrix{i,1};
                    disp("=============================================")
                    fprintf("Added Subject: %s Experiment: Grip", strMatrix{i,1})
                    disp("=============================================")
                catch
                    disp("Error \n")
                    continue
                end
            elseif k == 2
                disp("Found Nback ***************")
                try
                    nback_epochs{1,i} = ft_preprocessing(cfg);

                    event_file = mne_read_events(fiff_file);
                    nback_epochs{1,i}.trialinfo = event_file;
                    nback_epochs{2,i} = status;
                    nback_epochs{3,i} = strMatrix{i,1};
                    disp("=============================================")
                    fprintf("Added Subject: %s Experiment: Nback", strMatrix{i,1})
                    disp("=============================================")
                catch
                    disp("Error \n")
                    continue
                end
            elseif k == 3
                disp("Found Oddball ***************")
                try
                    oddball_epochs{1,i} = ft_preprocessing(cfg);

                    event_file = mne_read_events(fiff_file);
                    oddball_epochs{1,i}.trialinfo = event_file;
                    oddball_epochs{2,i} = status;
                    oddball_epochs{3,i} = strMatrix{i,1};
                    disp("=============================================")
                    fprintf("Added Subject: %s Experiment: Oddball \n", strMatrix{i,1})
                    disp("=============================================")
                catch
                    disp("Error \n")
                    continue
                end
            end
        end
    end
end

grip_epochs(:, any(cellfun(@isempty, grip_epochs), 1)) = [];
nback_epochs(:, any(cellfun(@isempty, nback_epochs), 1)) = [];
oddball_epochs(:, any(cellfun(@isempty, oddball_epochs), 1)) = [];

channels = grip_epochs{1,1}.label;

clc;

data = {grip_epochs; nback_epochs; oddball_epochs};

% Get epoch averages
for exp = 1:3
    for subj = 1:length(data{exp})
        meanOddCellColumns = zeros(size(data{exp}{1,subj}.trial{1}));
        meanEvenCellColumns = zeros(size(data{exp}{1,subj}.trial{1}));
        for i = 1:numel(data{exp}{1, subj}.trial)
            if(mod(i,2)==1)
                meanOddCellColumns = meanOddCellColumns + data{exp}{1,subj}.trial{i};
            else
                meanEvenCellColumns = meanEvenCellColumns + data{exp}{1,subj}.trial{i};
            end
        end
        numOddCellColumns = sum(mod(1:numel(data{exp}{1, subj}.trial), 2) == 1);
        meanOddCellColumns = meanOddCellColumns / numOddCellColumns;

        numEvenCellColumns = sum(mod(1:numel(data{exp}{1, subj}.trial), 2) == 0);
        meanEvenCellColumns = meanEvenCellColumns / numEvenCellColumns;

        data{exp}{1,subj}.trial = [{meanOddCellColumns},{meanEvenCellColumns}];
    end
end

% Get ROI averages
for exp = 1:3
    for subj = 1:length(data{exp})
        for cond = 1:2
            signal = data{exp}{1,subj}.trial{1,cond};
            ROI = setROI(signal,data{exp}{1,subj}.label);
            data{exp}{4,subj}{1,cond} = ROI;
%             for roi = 1:length(ROI)
%                 if cond == 1
%                     if (max(ROI{roi,1}) >= 1.5*10^-7) |  (min(ROI{roi,1}) <= -1.5*10^-7)
%                         reject(end+1) = data{exp}{3,subj};
%                         break;
%                     end
%                 end
%                 data{exp}{4,subj}{1,cond}{roi,1} = ROI{roi,1};
%                 data{exp}{4,subj}{1,cond}{roi,2} = ROI{roi,2};
%             end
        end
    end
end

clc;

% Plot Data
set(0,'DefaultFigureVisible','off');
titl = ["Rest-Control","Block-Control","Rest-Sarco","Block-Sarco";...
        "O-Back-Control","2-Back-Control","O-Back-Sarco","2-Back-Sarco";...
        "Std-Control","Odd-Control","Std-Sarco","Odd-Sarco"];
exps = ["Grip","N-Back","Oddball"];
for exp = 1:3
    for roi = 1:13
        figure("Name",sprintf("%s %s",exps(exp),data{exp}{4,1}{1,1}{roi,2}));
        t = tiledlayout(2,2);
        ax1 = nexttile;
        ax2 = nexttile;
        ax3 = nexttile;
        ax4 = nexttile;
        axs = [ax1,ax3,ax2,ax4];
        for subj = 1:length(data{exp})
            if strcmp(data{exp}{4,1}{1,1}{roi,2},data{exp}{4,subj}{1,1}{roi,2})
                for cond = 1:2
                    if strcmp(data{exp}{2,subj},"Control")
                        plot(axs(cond),data{exp}{4,subj}{1,cond}{roi,1},"r");
                        ylim(axs(cond),[(-2*10^-7) (2*10^-7)]);
                        title(axs(cond),titl(exp,cond));
                        hold(axs,"on");
                    else
                        plot(axs(cond+2),data{exp}{4,subj}{1,cond}{roi,1},"r");
                        ylim(axs(cond+2),[(-2*10^-7) (2*10^-7)]);
                        title(axs(cond+2),titl(exp,cond+2));
                        hold(axs,"on");
                    end
                end
            end
        end
        f = gcf;
        name = sprintf('/Users/boramert/Desktop/Yüksek Lisans/Matlab_Kod/SpaghettiPlot/%s %s.png',exps(exp),data{exp}{4,1}{1,1}{roi,2});
        exportgraphics(f,name,'Resolution',300);
    end
end
close all;
set(0,'DefaultFigureVisible','on');



function ROI = addToROI(data,ROI,index)
disp(inputname(1));
if(~isempty(data))
    ROI{index,1} = data;
    ROI{index,2} = inputname(1);
end
end

function ROI = setROI(signal,label)
% Primary Motor Cortex (BA4)
% Primary_Motor_Cortex_L = [1,1]
% Primary_Motor_Cortex_R = [5,5]
%
% Premotor and Supplementary Motor Cortex (BA6)
% Pre_Supplementary_Motor_Cortex_L = [2,2],[2,4],[4,2],[4,4],[3,2]
% Pre_Supplementary_Motor_Cortex_R = [5,5],[6,5],[5,6],[6,6],[7,6]
%
% Dorselateral Prefrontal Cortex (BA46 & BA9)
% Drosolateral_Prefrontal_Cortex_L = [4,4],[4,3],[2,2],[4,2],[2,4]
% Dorsolateral_Prefrontal_Cortex_R = [8,8],[8,7],[8,6],[6,8]

ROI = cell(13,2);

PMC_L_hbo = signal(contains(label,"S1_D1 hbo"),:);
PMC_L_hbr = signal(contains(label,"S1_D1 hbr"),:);
PMC_R_hbo = signal(contains(label,"S5_D5 hbo"),:);
PMC_R_hbr = signal(contains(label,"S1_D1 hbr"),:);
PSMC_L_hbo = mean(signal(contains(label,["S2_D2 hbo","S2_D4 hbo","S4_D2 hbo","S4_D4 hbo","S3_D2 hbo"]),:));
PSMC_L_hbr = mean(signal(contains(label,["S2_D2 hbr","S2_D4 hbr","S4_D2 hbr","S4_D4 hbr","S3_D2 hbr"]),:));
PSMC_R_hbo = mean(signal(contains(label,["S5_D5 hbo","S6_D5 hbo","S5_D6 hbo","S6_D6 hbo","S7_D6 hbo"]),:));
PSMC_R_hbr = mean(signal(contains(label,["S5_D5 hbr","S6_D5 hbr","S5_D6 hbr","S6_D6 hbr","S7_D6 hbr"]),:));
DPC_L_hbo = mean(signal(contains(label,["S4_D4 hbo","S4_D3 hbo","S2_D2 hbo","S4_D2 hbo","S2_D4 hbo"]),:));
DPC_L_hbr = mean(signal(contains(label,["S4_D4 hbr","S4_D3 hbr","S2_D2 hbr","S4_D2 hbr","S2_D4 hbr"]),:));
DPC_R_hbo = mean(signal(contains(label,["S8_D8 hbo","S8_D7 hbo","S8_D6 hbo","S6_D8 hbo"]),:));
DPC_R_hbr = mean(signal(contains(label,["S8_D8 hbr","S8_D7 hbr","S8_D6 hbr","S6_D8 hbr"]),:));

ROI = addToROI(PMC_L_hbo,ROI,1); %#ok<*STRSCALR>
ROI = addToROI(PMC_L_hbr,ROI,2);
ROI = addToROI(PMC_R_hbo,ROI,3);
ROI = addToROI(PMC_R_hbr,ROI,4);
ROI = addToROI(PSMC_L_hbo,ROI,5);
ROI = addToROI(PSMC_L_hbr,ROI,6);
ROI = addToROI(PSMC_R_hbo,ROI,7);
ROI = addToROI(PSMC_R_hbr,ROI,8);
ROI = addToROI(DPC_L_hbo,ROI,9);
ROI = addToROI(DPC_L_hbr,ROI,10);
ROI = addToROI(DPC_R_hbo,ROI,11);
ROI = addToROI(DPC_R_hbr,ROI,12);
end
