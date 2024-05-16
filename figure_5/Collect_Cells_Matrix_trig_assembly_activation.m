%% Kernels for smothing the binned z-scored firing rate
sigma_for_trial=0.005;
bin_trial=0.001;
edges_trial=(-0.5:0.015:0.5);
edges_norm_for_trial=(-0.5*sigma_for_trial:bin_trial:0.5*sigma_for_trial);
kernel_trial = normpdf(edges_norm_for_trial,0,sigma_for_trial);
kernel_trial=kernel_trial*bin_trial;

sigma = .005;                            %Standard deviation of the kernel 5-15 ms
edges_norm=[-0.5*sigma:.001:0.5*sigma];     %Time ranges form -X0*st. dev. to X1*st. dev.
kernel = normpdf(edges_norm,0,sigma);
kernel = kernel*.001;

edges=(-0.5:0.005:0.5); %5ms bin;
edges_large=(-0.5:0.025:0.5); %25ms bin

if ~isempty(selected_assembly_plot)
    loopRange = selected_assembly_plot;
else
    loopRange = 1:size(AssemblyTemplates,2);
end

CellID_counter=1;
for assembly_num=loopRange
    disp(['Assembly: ' num2str(assembly_num)]);
    for num_structure=1:length(unique(PutativeStructure))
        sorted_neurons_for_Assembly_sig=intersect(find(Significant_Neurons_Assemblies(:,assembly_num)),find(ismember(PutativeStructure,unique_PutativeStructure(num_structure))));
        sorted_neurons_for_Assembly_non_sig=intersect(find(~Significant_Neurons_Assemblies(:,assembly_num)),find(ismember(PutativeStructure,unique_PutativeStructure(num_structure)))');

        center= time_bins(find(Activities_Assembly(assembly_num,:)));
        time_on=  center-0.5; 
        time_off=  center+0.5;
        Activitymatrix_firing_Total_sig=[];
        Activitymatrix_firing_Total_non_sig=[];


        for num_sig=sorted_neurons_for_Assembly_sig'

            trial_ripple_ON=1;
            SpikeSecund=[];
            SpikeSecund=Activitymatrix_for_firing{num_sig,:}(:,:);

            if length(SpikeSecund)>100
                spikes_trials=zeros(length(center),length(edges));
                spikes_trials_assembly=zeros(length(center),length(edges_trial));
                spikes_trials_assembly_random=zeros(length(center),length(edges_trial));
                spikes_trials_assembly_Gauss=zeros(length(center),length(edges_trial));
                spikes_trials_assembly_Gauss_rand=zeros(length(center),length(edges_trial));

                for trial=1:length(center)

                    if ~isempty(SpikeSecund (find(time_on(trial) <= SpikeSecund   & SpikeSecund   <= time_off(trial))))
                        spike(trial_ripple_ON,:).time=(SpikeSecund(find(time_on(trial) <= SpikeSecund & SpikeSecund <=time_off(trial))))-(center(trial));
                        spike(trial_ripple_ON,:).assembly_ID=[trial, assembly_num];
                        spikes_gauss(trial_ripple_ON,:)=(histc(spike(trial_ripple_ON,:).time,edges)');
                        spikes_trials(trial,:)=(histc(spike(trial_ripple_ON,:).time,edges)');
                        spikes_trials_assembly(trial,:)=(histc(spike(trial_ripple_ON,:).time,edges_trial)');
                        s=conv(spikes_trials_assembly(trial,:),kernel_trial);
                        center_gauss= ceil(length(edges_trial)/2);
                        spikes_trials_assembly_Gauss(trial,:)=s(ceil(length(s)/2)-( center_gauss-1):ceil(length(s)/2)+( center_gauss-1));
                        trial_ripple_ON=trial_ripple_ON+1;

                    end

                end

                Cells_assembly(CellID_counter).individual_trials_count_assembly=  spikes_trials_assembly;
                Cells_assembly(CellID_counter).individual_trials_gauss_assembly=  spikes_trials_assembly_Gauss;
                Cells_assembly(CellID_counter).individual_trials_count_assembly_rand=    spikes_trials_assembly_random;
                Cells_assembly(CellID_counter).individual_trials_gauss_assembly_rand=   spikes_trials_assembly_Gauss_rand;
                Cells_assembly(CellID_counter).AssemblyID= assembly_num;
                Cells_assembly(CellID_counter).Significance= 1;
                Cells_assembly(CellID_counter).Structure= unique_PutativeStructure(num_structure);
                CellID_counter=CellID_counter+1;

            end

        end


        for num_sig= sorted_neurons_for_Assembly_non_sig'

            trial_ripple_ON=1;
            SpikeSecund=[];
            SpikeSecund=Activitymatrix_for_firing{num_sig,:}(:,:);

            if length(SpikeSecund)>100
                spikes_trials=zeros(length(center),length(edges));
                spikes_trials_assembly=zeros(length(center),length(edges_trial));
                spikes_trials_assembly_random=zeros(length(center),length(edges_trial));
                spikes_trials_assembly_Gauss=zeros(length(center),length(edges_trial));
                spikes_trials_assembly_Gauss_rand=zeros(length(center),length(edges_trial));


                for trial=1:length(center)

                    if ~isempty(SpikeSecund (find(time_on(trial) <= SpikeSecund   & SpikeSecund   <= time_off(trial))))
                        spike(trial_ripple_ON,:).time=(SpikeSecund(find(time_on(trial) <= SpikeSecund & SpikeSecund <=time_off(trial))))-(center(trial));
                        spike(trial_ripple_ON,:).assembly_ID=[trial, assembly_num];
                        spikes_gauss(trial_ripple_ON,:)=(histc(spike(trial_ripple_ON,:).time,edges)');
                        spikes_trials(trial,:)=(histc(spike(trial_ripple_ON,:).time,edges)');
                        spikes_trials_assembly(trial,:)=(histc(spike(trial_ripple_ON,:).time,edges_trial)');
                        s=conv(spikes_trials_assembly(trial,:),kernel_trial);
                        center_gauss= ceil(length(edges_trial)/2);
                        spikes_trials_assembly_Gauss(trial,:)=s(ceil(length(s)/2)-( center_gauss-1):ceil(length(s)/2)+( center_gauss-1));
                        trial_ripple_ON=trial_ripple_ON+1;

                    end

                end

                Cells_assembly(CellID_counter).individual_trials_count_assembly=  spikes_trials_assembly;
                Cells_assembly(CellID_counter).individual_trials_gauss_assembly=  spikes_trials_assembly_Gauss;
                Cells_assembly(CellID_counter).individual_trials_count_assembly_rand= spikes_trials_assembly_random;
                Cells_assembly(CellID_counter).individual_trials_gauss_assembly_rand= spikes_trials_assembly_Gauss_rand;
                Cells_assembly(CellID_counter).AssemblyID= assembly_num;
                Cells_assembly(CellID_counter).Significance= 0;
                Cells_assembly(CellID_counter).Structure= unique_PutativeStructure(num_structure);
                CellID_counter=CellID_counter+1;
              
            end

        end
    end
end
