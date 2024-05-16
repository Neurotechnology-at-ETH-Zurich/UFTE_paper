
fig_id=3;
statistic_counter=1;
for assembly_num=unique([Cells_assembly.AssemblyID])
    disp(['Assembly: ' num2str(assembly_num)]);
    figure(fig_id);
    sub_counter=1;
    for num_structure=1:length(unique(PutativeStructure))
        sub_bar1=subplot(length(unique(PutativeStructure)),2,sub_counter);

        significant=find(and([Cells_assembly.AssemblyID]==assembly_num, ([Cells_assembly.Significance]))); %Significant means "member of an assembly"
        significant=significant(ismember([Cells_assembly(significant).Structure],unique_PutativeStructure(num_structure)));

        if ~isempty(significant)

            non_significant=find(and([Cells_assembly.AssemblyID]==assembly_num, ([Cells_assembly.Significance]==0)));
            non_significant=non_significant(ismember([Cells_assembly(non_significant).Structure],unique_PutativeStructure(num_structure)));

            Reshaped_plotting_sig=reshape([Cells_assembly(significant).individual_trials_gauss_assembly],[size(Cells_assembly(significant(1)).individual_trials_gauss_assembly,1) size(Cells_assembly(significant(1)).individual_trials_gauss_assembly,2)*length(significant)]);
            Reshaped_plotting_non_sig=reshape([Cells_assembly(non_significant).individual_trials_gauss_assembly],[size(Cells_assembly(non_significant(1)).individual_trials_gauss_assembly,1) size(Cells_assembly(non_significant(1)).individual_trials_gauss_assembly,2)*length(non_significant)]);

            errBar_SIG=(reshape(std(Reshaped_plotting_sig),length(mean(Reshaped_plotting_sig))/length(significant),length(significant)))';

            Reshaped_plotting_sig_mean=reshape(mean(Reshaped_plotting_sig),length(mean(Reshaped_plotting_sig))/length(significant),length(significant))';
            Reshaped_plotting_sig_mean=(Reshaped_plotting_sig_mean-mean(Reshaped_plotting_sig_mean(:,4:16),2))./errBar_SIG;
            errBar_non_SIG=(reshape(std(Reshaped_plotting_non_sig),length(mean(Reshaped_plotting_non_sig))/length(non_significant),length(non_significant)))';


            Reshaped_plotting_non_sig_mean=reshape(mean(Reshaped_plotting_non_sig),length(mean(Reshaped_plotting_non_sig))/length(non_significant),length(non_significant))';
            Reshaped_plotting_non_sig_mean=(Reshaped_plotting_non_sig_mean-mean(Reshaped_plotting_non_sig_mean(:,4:16),2))./errBar_non_SIG;

            if  ~isempty(find(sum(isnan(Reshaped_plotting_non_sig_mean)')))
                Reshaped_plotting_non_sig_mean=Reshaped_plotting_non_sig_mean(find(sum((isnan(Reshaped_plotting_non_sig_mean))')==0),:);
                errBar_non_SIG= errBar_non_SIG(find(sum((isnan(Reshaped_plotting_non_sig_mean))')==0),:);
            end

            if  ~isempty(find(sum(isinf(Reshaped_plotting_non_sig_mean)')))
                Reshaped_plotting_non_sig_mean=Reshaped_plotting_non_sig_mean(find(sum((isinf(Reshaped_plotting_non_sig_mean))')==0),:);
                errBar_non_SIG= errBar_non_SIG(find(sum((isinf(Reshaped_plotting_non_sig_mean))')==0),:);
            end

            if size(Reshaped_plotting_sig_mean,1)~=1

                patch_data_sig=(mean(Reshaped_plotting_sig_mean(:,4:64))); %((mean(mean(Reshaped_plotting_sig(:,4:64,:),3),1))-mean(mean(mean(Reshaped_plotting_sig(:,4:16,:),3))))./mean(std(mean(Reshaped_plotting_sig(:,4:16,:),3)))+mean(std(mean(Reshaped_plotting_sig(:,4:16,:),3)));
                errBar_sig=std(errBar_SIG(:,4:64)); %mean(std(mean(Reshaped_plotting_sig(:,4:16,:),3)));
                hold on;
                plot(edges_trial(4:64),patch_data_sig,'Color', color_for_stem(num_structure,:));
                p1 =   patch([edges_trial(4:64) flip(edges_trial(4:64))], [patch_data_sig+errBar_sig flip(patch_data_sig-errBar_sig)],color_for_stem(num_structure,:));
                p1.EdgeColor=color_for_stem(num_structure,:);
                p1.FaceColor=color_for_stem(num_structure,:);
                p1.FaceAlpha=0.3;
                p1.LineStyle='--';

                ylim([-0.05 max(patch_data_sig+errBar_sig)+0.01]);
                % if  assembly_num==8
                %     ylim([-0.05 0.7]);
                % else
                %     ylim([-0.05 0.2]);
                % end


                patch_data_non_sig=mean(Reshaped_plotting_non_sig_mean(:,4:64)); %((mean(mean(Reshaped_plotting_non_sig(:,4:64,:),3),1))-mean(mean(mean(Reshaped_plotting_non_sig(:,4:16,:),3))))./mean(std(mean(Reshaped_plotting_non_sig(:,4:16,:),3)));
                errBar_non_sig=std(errBar_non_SIG(:,4:64));%std(Reshaped_plotting_non_sig_mean(:,4:64)); %mean(std(mean(Reshaped_plotting_non_sig(:,4:16,:),3)));

                plot(edges_trial(4:64),patch_data_non_sig,'Color',[17 17 17]/255);
                hold on;
                p =   patch([edges_trial(4:64) flip(edges_trial(4:64))], [patch_data_non_sig+  errBar_non_sig flip(patch_data_non_sig-  errBar_non_sig)], [0.800000011920929 0.800000011920929 0.800000011920929]);
                p.EdgeColor=[0.800000011920929 0.800000011920929 0.800000011920929];
                p.FaceColor= [0.800000011920929 0.800000011920929 0.800000011920929];
                p.FaceAlpha=0.1;
                p.LineStyle='--';
                xlim([-0.4 0.4]);

                axis square
                % Plot

                sub_counter=sub_counter+1;
                sub_bar=subplot(length(unique(PutativeStructure)),2,sub_counter);

                % Data
                y = [mean(mean(Reshaped_plotting_non_sig_mean(:,32:38))); mean(mean(Reshaped_plotting_sig_mean(:,32:38)))]'; % first 3 #s are pre-test, second 3 #s are post-test
                err = [mean(std(errBar_non_SIG(:,32:38))); mean(std(errBar_SIG(:,32:38)))]';


                [h,p,ci,stats] =ttest2((mean(Reshaped_plotting_non_sig_mean(:,32:38))),(mean(Reshaped_plotting_sig_mean(:,32:38))));

                STAT(statistic_counter).AssemblyID=assembly_num;
                STAT(statistic_counter).Structura=unique_PutativeStructure(num_structure);
                STAT(statistic_counter).h=h;
                STAT(statistic_counter).p=p;
                STAT(statistic_counter).tstat=stats.tstat;
                STAT(statistic_counter).df=stats.df;
                STAT(statistic_counter).sd=stats.sd;
                STAT(statistic_counter).SD=[err];
                STAT(statistic_counter).mean=[mean(Reshaped_plotting_non_sig_mean(:,32:38),'all'),mean(Reshaped_plotting_sig_mean(:,32:38),'all')];

                clear hb
                hb(1) = bar(1,y(1));
                % get the bar handles
                hold on;
                hb(2)  = bar(2,y(2));
                hb(1).BarWidth=0.5;
                hb(2).BarWidth=0.5;

                hold on;


                for k =1:size(y,2)
                    % get x positions per group
                    xpos = hb(k).XData + hb(k).XOffset;
                    % draw errorbar
                    errorbar(xpos, y(:,k), err(:,k), 'LineStyle', 'none', ...
                        'Color', 'k', 'LineWidth', 1);
                end

                % Set Axis properties
                sub_bar.XTick=[1 2];
                sub_bar.XTickLabel={'non-members'; 'members'};
                ylabel('z-scored firing rate');
                xlim([0.5 2.5]);
                if  assembly_num==8
                    ylim([-0.05 0.7]);
                else
                    ylim([-0.05 0.2]);
                end

                axis square

                clear    sub_bar
            else

                plot(edges_trial(4:64),((Reshaped_plotting_sig_mean(:,4:64))),'Color',color_for_stem(num_structure,:));
                hold on; plot(edges_trial(4:64),(mean( Reshaped_plotting_non_sig_mean(:,4:64))),'Color',[0.800000011920929 0.800000011920929 0.800000011920929]);
                patch_data_sig=((Reshaped_plotting_sig_mean(:,4:64))); %((mean(mean(Reshaped_plotting_sig(:,4:64,:),3),1))-mean(mean(mean(Reshaped_plotting_sig(:,4:16,:),3))))./mean(std(mean(Reshaped_plotting_sig(:,4:16,:),3)))+mean(std(mean(Reshaped_plotting_sig(:,4:16,:),3)));
                errBar_sig=(errBar_SIG(:,4:64)); %std(Reshaped_plotting_sig_mean(:,4:64)) %mean(std(mean(Reshaped_plotting_sig(:,4:16,:),3)));
                hold on;
                plot(edges_trial(4:64),patch_data_sig,'Color', color_for_stem(num_structure,:));
                p1 =   patch([edges_trial(4:64) flip(edges_trial(4:64))], [patch_data_sig+errBar_sig flip(patch_data_sig-errBar_sig)],color_for_stem(num_structure,:));
                p1.EdgeColor=color_for_stem(num_structure,:);
                p1.FaceColor=color_for_stem(num_structure,:);
                p1.FaceAlpha=0.3;
                p1.LineStyle='--';

                ylim([-0.05 max(patch_data_sig+errBar_sig)+0.01]);
                % if  assembly_num==8
                %     ylim([-0.05 0.7])
                % else
                %     ylim([-0.05 0.2])
                % end


                patch_data_non_sig=mean(Reshaped_plotting_non_sig_mean(:,4:64)); %((mean(mean(Reshaped_plotting_non_sig(:,4:64,:),3),1))-mean(mean(mean(Reshaped_plotting_non_sig(:,4:16,:),3))))./mean(std(mean(Reshaped_plotting_non_sig(:,4:16,:),3)));
                errBar_non_sig=std(errBar_non_SIG(:,4:64)); %std(Reshaped_plotting_non_sig_mean(:,4:64)); %mean(std(mean(Reshaped_plotting_non_sig(:,4:16,:),3)));

                plot(edges_trial(4:64), patch_data_non_sig,'Color',[17 17 17]/255);
                hold on;

                p =   patch([edges_trial(4:64) flip(edges_trial(4:64))], [patch_data_non_sig+ errBar_non_sig flip(patch_data_non_sig-errBar_non_sig)], [0.800000011920929 0.800000011920929 0.800000011920929]);
                p.EdgeColor=[0.800000011920929 0.800000011920929 0.800000011920929];
                p.FaceColor= [0.800000011920929 0.800000011920929 0.800000011920929];
                p.FaceAlpha=0.1;
                p.LineStyle='--';

                xlim([-0.4 0.4]);

                axis square;

                sub_counter=sub_counter+1;
                sub_bar=subplot(length(unique(PutativeStructure)),2,sub_counter);

                % Data
                y = [mean(mean(Reshaped_plotting_non_sig_mean(:,32:38))); mean(mean(Reshaped_plotting_sig_mean(:,32:38)))]'; % first 3 #s are pre-test, second 3 #s are post-test
                err = [mean(std(errBar_non_SIG(:,32:38))); mean(std(errBar_SIG(:,32:38)))]';
                [h,p,ci,stats] =ttest2((mean(Reshaped_plotting_non_sig_mean(:,32:38))),(mean(Reshaped_plotting_sig_mean(:,32:38))));

                STAT(statistic_counter).AssemblyID=assembly_num;
                STAT(statistic_counter).Structura=unique_PutativeStructure(num_structure);
                STAT(statistic_counter).h=h;
                STAT(statistic_counter).p=p;
                STAT(statistic_counter).tstat=stats.tstat;
                STAT(statistic_counter).df=stats.df;
                STAT(statistic_counter).sd=stats.sd;
                STAT(statistic_counter).SD=[err];
                STAT(statistic_counter).mean=[mean(Reshaped_plotting_non_sig_mean(:,32:38),'all'),mean(Reshaped_plotting_sig_mean(:,32:38),'all')];

                clear hb;
                hb(1) = bar(1,y(1));
                hold on;
                hb(2)  = bar(2,y(2));
                hold on;
                axis square
                hb(1).BarWidth=0.5;
                hb(2).BarWidth=0.5;


                for k =1:size(y,2)
                    % get x positions per group
                    xpos = hb(k).XData + hb(k).XOffset;
                    % draw errorbar
                    errorbar(xpos, y(:,k), err(:,k), 'LineStyle', 'none', ...
                        'Color', 'k', 'LineWidth', 1);
                end

                % Set Axis properties
                sub_bar.XTick=[1 2];
                sub_bar.XTickLabel={'non-members'; 'members'};
                ylabel('z-scored firing rate')
                xlim([0.5 2.5])
                ylim([-0.05 max(patch_data_sig+errBar_sig)+0.01])
                % if  assembly_num==8
                %     ylim([-0.05 0.7]);
                % else
                %     ylim([-0.05 0.2]);
                % end
                axis square
                clear    sub_bar

            end

            try
                title(['Assembly: ' num2str(assembly_num) ' Structure:' char(unique([Cells_assembly(significant).Structure])) ' n=' num2str(size(Cells_assembly(significant(1)).individual_trials_gauss_assembly,1))]);
            catch
            end
            sub_counter=sub_counter+1;

            clear p p1

            statistic_counter=statistic_counter+1;
        end

    end


    fig=gcf;
    fig.Renderer = 'painters';
    fig.PaperUnits = 'points';
    fig.PaperPosition = [0 0 600 800];
    fig.PaperSize = [600 800];
    % saveas(fig,['Assembly_Firing_Rate_' num2str( fig_id)],'svg')
    % saveas(fig,['Assembly_Firing_Rate_' num2str( fig_id)],'tif')
    % saveas(fig,['Assembly_Firing_Rate_' num2str( fig_id)],'fig')
    fig_id=fig_id+1;
end
