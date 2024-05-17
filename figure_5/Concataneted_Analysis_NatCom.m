%% Load information about the raw data for rat37

filename = 'combined.dat'; % Concatenated Raw Data filename
ratID='rTBY37';
num_files = 1;
FileInfo = dir(filename);
nChansInRawFile = 256;

if ~isempty(FileInfo)  % If the raw data is accessible
    num_samples = FileInfo.bytes / (nChansInRawFile * 2);
else  % If the raw data file is missing, use this default value to estimate the length of raw data in samples.
    num_samples = 595698304;
end

%% Load preprocessed variables
% Loading the results of spike sorting from concatenated sessions.

try
    load('Github_Data.mat');
catch
disp('Download row spikes from Drobpox')
% Define the URL for direct download
url = 'https://www.dropbox.com/scl/fi/moxem3u5df546wofs9jxz/Github_Data.mat?rlkey=oadoony6j03izk1koasauozjt&dl=1';
outputFileName = 'downloaded_Github_Data.mat';
% Download the file
websave(outputFileName, url);
% Load all variables from the .mat file into the workspace
load(outputFileName);
end

% Activities - [number of assemblies x assembly activation strength per bin (bin = 25 ms)]
% ActivityMatrix - [all individual units x z-scored firing rate]
% AssemblyTemplates - [all individual neurons x assembly number (weight)]
% r37 - [detected units per assembly x number of sessions]
% Ripple.TimeStamp - [number of ripples x (start, stop) times] in seconds; .eventID - 
% label for the location where the ripple was detected; .centerTimepoints - center
% of the Sharp Wave Ripple determined at the maximum power of ripples (150-200 Hz).
% Spike - Results of spike sorting from concatenated data
% spikeClusters - Cluster assignment for each spike
% SpikeSites - Site with the peak spike amplitude (central site)
% spikeTimes - Timing of spikes in ADC sample units
% waveform - Spike waveforms for each spike

%% Constants

bin = 25; % to construct binned z-scored firing rate matrix
sample_rate = 20000; %Hz
real_time_scale = [3 7 9 14 17 21 24 28 31 35 37 42 44 49 52 56 59 63 66 70 76 80 84 87 90]; % Days for the sessions

%Several time scales are necessary for binning the concatenated data while
%keeping track of the beginning and the end of the session. The first
%session is shorter that is the resion of extraction from the total sample
% The total length of the 1+24 (1st session is shorter then the others) sessions is calculated as sample = (24 sessions * 1200 sec * 20000 sample rate). 
% This value needs to be extracted from the total size of the concatenated sessions to determine the length 
% of the first session, which was only 984.9152 sec instead of the 1200 sec.

time_bins=linspace(0,(num_samples./sample_rate),ceil((num_samples./sample_rate)/(bin/1000)));
time_absolute=(num_samples-(24*1200*20000):(1200*20000):num_samples);
center_of_relative_times=([0, (0:1200:(24*1200))+((num_samples-(24*1200*20000))./20000)]);
center_of_relative_times_hist=([0, (0:1200:(24*1200))+((num_samples-(24*1200*20000))./20000)]);
Recording_CellID = unique(spikeClusters(spikeClusters > 0)); % Cluster ID<=0 represent noise.
Assembly_cellID = intersect(find(ismember(Spike.clusterNotes(num_files).clusterNotes, {'single'})), Recording_CellID); % Selecting single units; each cluster is labeled: "multi" or "single"
binnend_activity_strength=[linspace(0,((num_samples-(24*1200*20000)-1)./20000),ceil((num_samples-(24*1200*20000)-1)./20000)) linspace((num_samples-(24*1200*20000))./20000,(num_samples./sample_rate),24*1200)];

%% This is only for rat 37 based on anatomical confirmation!

PutativeStructure = {};
PutativeStructure(spikeSites(Assembly_cellID) <= 64) = {'iHP'};
PutativeStructure(spikeSites(Assembly_cellID) > 64 & spikeSites(Assembly_cellID) <= 128) = {'dHP'};
PutativeStructure(spikeSites(Assembly_cellID) > 128 & spikeSites(Assembly_cellID) <= 192) = {'RSC'};
PutativeStructure(spikeSites(Assembly_cellID) > 192 & spikeSites(Assembly_cellID) <= 202) = {'IL'};
PutativeStructure(spikeSites(Assembly_cellID) > 202 & spikeSites(Assembly_cellID) <= 220) = {'PrL'};
PutativeStructure(spikeSites(Assembly_cellID) > 220) = {'Cg1'};

%% Construction matrix for assembly analysis -> binning data

opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 1000;

Activitymatrix = zeros(length(Assembly_cellID), length(time_bins));
for Cell_ID = 1:length(Assembly_cellID)
    Activitymatrix(Cell_ID, :) = histc(spikeTimes(spikeClusters == (Assembly_cellID(Cell_ID))) ./ sample_rate, time_bins);
    %    Activitymatrix_times(Cell_ID, :) = {spikeTimes(spikeClusters == (Assembly_cellID(Cell_ID)))};
end


%% Selection for the paper (Inter areal and Intra areal ensemble)

%selected_assembly_plot=[];% plot all assemblies
selected_assembly_plot=[3 20 12 18 32 28]; 

%% Detecting cell assemblies in large neuronal populations. J Neurosci Methods 220(2):149-66. 10.1016/j.jneumeth.2013.04.010
% https://github.com/tortlab/Cell-Assembly-Detection/tree/master

if isempty(AssemblyTemplates) % Check if assembly detection results are already provided. If you recalculate the assembly patterns, be aware that the assembly IDs may be reshuffled and the signs of the weights can also be flipped.
    AssemblyTemplates = assembly_patterns(Activitymatrix);
    Activities = assembly_activity(AssemblyTemplates, Activitymatrix);
    selected_assembly_plot=[]; % In this case, the selection is not accurate because the assembly IDs are shuffled.
end

%% Otsu's Methode for threshold of the weight of the neurons/assembly (soring out the memebrs vs. non-members)

figure(1)
unique_PutativeStructure = unique(PutativeStructure);
color_for_stem = colormap("lines");

disp ('Plot Figure 5a')
subplot_num=1;

if ~isempty(selected_assembly_plot)
    loopRange = selected_assembly_plot;
else
    loopRange = 1:size(AssemblyTemplates, 2);
end


for assembly_num=loopRange 

    significant(:,assembly_num)=graythresh(abs(AssemblyTemplates(:,assembly_num)));
    Significant_Neurons_Assemblies(:,assembly_num)=abs(AssemblyTemplates(:, assembly_num))>=significant(:,assembly_num);

    for num_structure=1:length(unique_PutativeStructure)
        sorted_neurons_for_Assembly_sig=intersect(find(Significant_Neurons_Assemblies(:,assembly_num)),find(ismember(PutativeStructure,unique_PutativeStructure(num_structure))));
        sorted_neurons_for_Assembly_non_sig=intersect(find(~Significant_Neurons_Assemblies(:,assembly_num)),find(ismember(PutativeStructure,unique_PutativeStructure(num_structure))));

        if ~isempty(selected_assembly_plot)
            subplot(2,size(AssemblyTemplates(:,selected_assembly_plot),2),subplot_num);
        else
            subplot(2,size(AssemblyTemplates,2),assembly_num);
        end

        if   ~isempty(sorted_neurons_for_Assembly_sig)
            stem(sorted_neurons_for_Assembly_sig,AssemblyTemplates(sorted_neurons_for_Assembly_sig,assembly_num),'filled','Color', color_for_stem(num_structure,:));
        end

        subtitle(num2str(assembly_num));
        hold on;
        stem(sorted_neurons_for_Assembly_non_sig,AssemblyTemplates(sorted_neurons_for_Assembly_non_sig,assembly_num),'Color', color_for_stem(num_structure,:));
        ylim([-0.6 0.6])
        xlim([0 size(AssemblyTemplates,1)+1])
        view(90,-90)
        box on;

    end
    if assembly_num==1
        xticks(1:size(AssemblyTemplates,1))
        xticklabels({Assembly_cellID})
        xlabel('Single Neurons Cell ID#')
    else
        xticks(1:size(AssemblyTemplates,1))
        xticklabels(PutativeStructure)
    end


    if ~isempty(selected_assembly_plot)
        subplot(2,size(AssemblyTemplates(:,selected_assembly_plot),2),subplot_num);
            subplot(2,size(AssemblyTemplates(:,selected_assembly_plot),2),subplot_num+size(AssemblyTemplates(:,selected_assembly_plot),2));
    else
        subplot(2,size(AssemblyTemplates,2),assembly_num+size(AssemblyTemplates,2));
    end



    PieChart_Properties=piechart(categorical(PutativeStructure(find(Significant_Neurons_Assemblies(:,assembly_num)))));
    PieChart_Properties.ColorOrder=(color_for_stem(find(ismember(unique_PutativeStructure,categorical(PutativeStructure(find(Significant_Neurons_Assemblies(:,assembly_num)))))),:));

    subplot_num=subplot_num+1;
end
sgtitle ("Weight of inter- and intra- areal (hippocampal and cortical) neuronal ensembles")

%% Plotting Figure 5c

figure(2)
disp('plot Figure 5c')

%% Assembly activation > 95th across concataneted sessions

Percentiles=2*std(Activities');
clear Activities_Assembly
for assembly_num=1:size(AssemblyTemplates,2)
    Activities_Assembly(assembly_num,:)=Activities(assembly_num,:)>Percentiles(1,assembly_num);
    %Activities_Assembly_l(assembly_num,:)=Activities(assembly_num,:)<Percentiles(1,assembly_num);
    %Activities_Assembly_h(assembly_num,:)=Activities(assembly_num,:)>Percentiles(2,assembly_num);
end


color_assemblies=colormap("parula");
color_assemblies=color_assemblies(1:floor(length(color_assemblies)./size(AssemblyTemplates,2)):size(color_assemblies,1),:);
bin_activities=binnend_activity_strength(1:60:length(binnend_activity_strength));
activity_strength=[];
subplot_num=1;

if ~isempty(selected_assembly_plot)
    loopRange = selected_assembly_plot;
else
    loopRange = 1:size(AssemblyTemplates, 2);
end


for assembly_num=loopRange 
    disp(['Assembly:' num2str(assembly_num)])

    if ~isempty(selected_assembly_plot);
        ax(subplot_num)=subplot(length(selected_assembly_plot),1,subplot_num);

    else isempty(selected_assembly_plot);
        ax(assembly_num)=subplot(size(AssemblyTemplates,2),1,assembly_num);
    end

    hold on;
    activity_strength=[];
    activity_strength_temp=[];
    group_time_vector_activty_strength = discretize(time_bins,bin_activities);
    activity_strength= grpstats(Activities(assembly_num,:),group_time_vector_activty_strength,"mean");
    Activity_strength(assembly_num,:)= activity_strength;

    yi = smooth(bin_activities(1:end-1),  Activity_strength(assembly_num,:),4);
    yyaxis left


    plot(bin_activities(1:end-1),yi,'LineWidth',2)
    hold on;  plot(bin_activities(yi>2*std(yi)), yi(yi>2*std(yi)),'.r','LineWidth',2)
    hold on; yline(2*std(yi),'--','LineWidth',2)
    ylabel(['Assembly ' num2str(assembly_num)])
    xlim([0 max(time_bins)])
    Activity_strength_smoothed(assembly_num,:)=yi;
    yyaxis right

    u = repelem(r37(assembly_num,:),2);
    pseudo_center_of_relative_times=zeros(1,length(repelem(r37(assembly_num,:),2)));
    pseudo_center_of_relative_times(1:2:length(repelem(r37(assembly_num,:),2)))=center_of_relative_times(1:end-1);
    pseudo_center_of_relative_times(2:2:length(repelem(r37(assembly_num,:),2)))=center_of_relative_times(2:end);
    plot(pseudo_center_of_relative_times,repelem(r37(assembly_num,:),2),'LineWidth',2) 
    ylim([0 1.1])

    xticks(center_of_relative_times(1:end-1)+([(diff(center_of_relative_times)./2)]))
    xticklabels(real_time_scale);hold on;
    xline(center_of_relative_times)

    subplot_num=subplot_num+1;
end

axes(ax(floor(numel(ax)/2)));
ylabel('Detectable neurons in ensemble (%)')

axes(ax(numel(ax)));
xlabel('Days')
sgtitle('Ensemble activation strength [a.u.]')


%% (figure 5b) Calculate the mean firing rate of members during periods when the
% assembly activation strength exceeds 2 standard deviations. 

for Cell_ID=1:length(Assembly_cellID)
    Activitymatrix_for_firing(Cell_ID,:)={spikeTimes(spikeClusters==(Assembly_cellID(Cell_ID)))./sample_rate};
end

%% Comparison of the the z-scored firing rates of neuronal members and non-members within the same ensemble triggered by ensemble activations. 
disp('Comparison of the the z-scored firing rates between members and non-members')
Collect_Cells_Matrix_trig_assembly_activation;

%% Plotting and claculate Statistic.
disp('Plot Figure 5b')
PlotFig5B;


%% Grabbing activation strength around Sharp Wave Ripples for all ensembles.
disp('Plot Figure 5e')
figure()
Relative_Time=linspace(0,num_samples./sample_rate,length(Activities));
time_on=Ripple.centerTimepoint-0.5;
time_off=Ripple.centerTimepoint+0.5;%(ripples.timestamps(:,2));
center= Ripple.centerTimepoint;

for darab_ripples=1:length(center);
    Assembly_Activation_During_Ripple_long(:,:,darab_ripples)=Activities(:,find(time_on(darab_ripples) <Relative_Time & Relative_Time < time_off(darab_ripples)));
end

if ~isempty(selected_assembly_plot);
    loopRange = selected_assembly_plot;
else
    loopRange = 1:size(Assembly_Activation_During_Ripple_long,1);
end
assembly_time=linspace(-500,500,40);
counter=1;
for assembly_num=loopRange %[3 20 12 18 32 28]%size(Assembly_Activation_During_Ripple,1)

    if ~isempty(selected_assembly_plot);
          subplot(length(selected_assembly_plot),1,counter);

    else isempty(selected_assembly_plot);
         subplot(size(Assembly_Activation_During_Ripple_long,1),1,counter);
    end

    patch_data_sig=squeeze(mean(Assembly_Activation_During_Ripple_long(assembly_num,:,:),3));
    errBar=(std(squeeze(Assembly_Activation_During_Ripple_long(assembly_num,:,:))'))./sqrt(size(squeeze(Assembly_Activation_During_Ripple_long(assembly_num,:,:)),2));
    plot(assembly_time,patch_data_sig,'LineWidth',2,'Color',color_for_stem( counter,:));

    p1 =   patch([assembly_time flip(assembly_time)], [patch_data_sig+errBar flip(patch_data_sig-errBar)],color_for_stem( counter,:));
    p1.FaceAlpha=0.3;
    p1.LineStyle='--';
    title(['Assembly ID:' num2str(assembly_num)]);
    axis square
    Assembly_activation_during_SWR_mean(counter,:)= patch_data_sig; % mean activation strengh/ensemble for fig5d
    Assembly_activation_during_SWR_SEM(counter,:)= errBar; % SEM activation strengh/ensemble for fig5d
    counter=counter+1;


end

%% Collect data for for python script to calculate lifetime of ensembly vs activation strengh of a given ensembly during shapr wave ripples fig5d-f
% We only focus on the activation strengh during Sharp Wave Ripple -25+25 around center. 

time_on=Ripple.centerTimepoint-0.025;
time_off=Ripple.centerTimepoint+0.025;%(ripples.timestamps(:,2));
center= Ripple.centerTimepoint;

%   data in Assembly_Activation_During_Ripple by a python script fig5df.py

for darab_ripples=1:length(center);
    Assembly_Activation_During_Ripple(:,darab_ripples)=mean(Activities(:,find(time_on(darab_ripples) <Relative_Time & Relative_Time < time_off(darab_ripples))),2);
end

disp('Additionally, we repeated this analysis and collected these values ( Assembly_Activation_During_Ripple_XX) from other rats for Figures 5D and 5F.')
