%% Generate multiple batches for a given Production phase. 
% IndPenSim_V2.01 Main file
% Last editted 9th of July 2019
% Objective: Generates batches for analysis
%% Copyright
% Stephen Goldrick Apr 2019 
% Contact: s.goldrick@ucl.ac.uk
% University College London, The University of Manchester, Newcastle University and Perceptive Engineering
% Please reference: "The Development of an Industrial Scale Fed-Batch
% Fermentation Simulation", Stepen Goldrick, Andrei Stefen, David Lovett,
% Gary Montague, Barry Lennox, Journal of Biotechnology 2015
%(https://www.sciencedirect.com/science/article/pii/S0168165614009377)
% And also please reference:  
% and "Modern day control challenges for industrial-scale fermentation
% processes" Goldrick et al. 2018 Computers and Chemical Engineering 
% All rights reserved. Copyright (c) The University of Manchester, University College London, Newcastle University and Perceptive Engineering.
% Clear workspace
 close all 
 clear all 
 clc

Data_generation_flag=1; % 1-  Generate data for a fixed time interval
                         % 2- Generate Fixed number of batches (generate
                         % Normal and Fault batches
Operational_days = 336; 
Bioreactor_turn_around_time =3;
Production_Phase_in_years = 1; % Selecting a Production phase in year (1 year =  365 days)
Save_batch_flag = 1; % 1- save data 0 - don't save data
Batches_File_Name  = 'IndPenSim_V2_export_V7_CampaignALL'; % File name for saving generated batches
Generate_batch_flag = 1;% 1- Generate Batches and Plot figures ,  0 -  Load data from workspace and plot figures

% Setup campaign configurations - 5 different campaigns
% Each campaign has different control strategies
num_campaigns = 5;

% Initialize the campaign data structure
All_campaign_data = struct();

for campaign_no = 1:num_campaigns
    % Reset batch data for this campaign
    Raw_Batch_data = struct();
    
    % Configure campaign parameters based on campaign number
    % Configure different operation modes for each campaign
    switch campaign_no
        case 1
            % Campaign 1: Operator controlled, Fixed batch length, No Raman
            control_strategy_val = 1; % Operator controlled 
            batch_length_val = 0; % Fixed batch length
            raman_spec_val = 0; % No Raman
            campaign_name = 'Campaign1';
        case 2
            % Campaign 2: Recipe driven, Variable batch length, No Raman
            control_strategy_val = 0; % Recipe driven
            batch_length_val = 1; % Variable batch length
            raman_spec_val = 0; % No Raman
            campaign_name = 'Campaign2';
        case 3
            % Campaign 3: Operator controlled, Variable batch length, No Raman
            control_strategy_val = 1; % Operator controlled
            batch_length_val = 1; % Variable batch length
            raman_spec_val = 0; % No Raman
            campaign_name = 'Campaign3';
        case 4
            % Campaign 4: Recipe driven, Fixed batch length, No Raman
            control_strategy_val = 0; % Recipe driven
            batch_length_val = 0; % Fixed batch length
            raman_spec_val = 0; % No Raman
            campaign_name = 'Campaign4';
        case 5
            % Campaign 5: Operator controlled, Fixed batch length, No Raman
            control_strategy_val = 1; % Operator controlled 
            batch_length_val = 0; % Fixed batch length
            raman_spec_val = 0; % No Raman
            campaign_name = 'Campaign5';
    end
    
    fprintf('Generating %s\n', campaign_name);
    
    % Automatic calculation of batches for time period
    % Start with batch 1 and continue until we exceed the time limit
    Batch_no = 0;
    Batch_start = 1;
    Summary_of_batch_lengths = 0;
    Summmary_of_campaign = [];
    
    Max_theoritical_number_of_batches =  round((Production_Phase_in_years*Operational_days)/11);
    
    % Set up Batch_run_flags with arrays of the right size
    Batch_run_flags.Batch_fault_order_reference = zeros(Max_theoritical_number_of_batches,1);
    Batch_run_flags.Control_strategy = control_strategy_val*ones(1,Max_theoritical_number_of_batches);
    Batch_run_flags.Batch_length = batch_length_val*ones(1,Max_theoritical_number_of_batches);
    Batch_run_flags.Raman_spec = raman_spec_val*ones(1,Max_theoritical_number_of_batches);
    
    % Store campaign flags for reference
    campaign_flags.Control_strategy = control_strategy_val;
    campaign_flags.Batch_length = batch_length_val;
    campaign_flags.Raman_spec = raman_spec_val;

    % Generate batches until we exceed the campaign duration
    while Summary_of_batch_lengths < Production_Phase_in_years*Operational_days && Batch_no < Max_theoritical_number_of_batches
        Batch_no = Batch_no + 1;
        matFileName = sprintf('Batch_%02d', Batch_no);
        
        % Run the IndPenSim simulation for this batch
        fprintf('Running Batch %d...\n', Batch_no);
        Xref = indpensim_run(Batch_no, Batch_run_flags);
        Raw_Batch_data.(matFileName) = Xref;
        
        % Collect batch summary statistics
        Summmary_of_campaign(Batch_no,:) = [
            Raw_Batch_data.(matFileName).Stats.Penicllin_harvested_during_batch, 
            Raw_Batch_data.(matFileName).Stats.Penicllin_harvested_end_of_batch,
            Raw_Batch_data.(matFileName).Stats.Penicllin_yield_total,  
            Raw_Batch_data.(matFileName).Fg.t(end)/24
        ];
        
        % Calculate cumulative batch lengths including turnaround time
        Summary_of_batch_lengths = sum(ceil(Summmary_of_campaign(:,4))) + Batch_no*Bioreactor_turn_around_time;
        
        fprintf('Batch %d completed. Total campaign time: %d days\n', Batch_no, Summary_of_batch_lengths);
    end
    
    % Final number of batches in this campaign
    Num_of_Batches = Batch_no;
    
    % Generate batch records for this campaign
    campaign_filename = sprintf('%s_%s', Batches_File_Name, campaign_name);
    Batch_Records = Generate_Batch_records(Raw_Batch_data,campaign_filename,Batch_run_flags);
    
    % Save campaign data
    if Save_batch_flag == 1
        savename = campaign_filename;
        save(savename, 'Batch_Records', 'Num_of_Batches', 'Batch_start', 'Summmary_of_campaign', 'Raw_Batch_data', 'campaign_flags');
    end
    
    % Store campaign data in overall structure
    All_campaign_data.(campaign_name).Batch_Records = Batch_Records;
    All_campaign_data.(campaign_name).Num_of_Batches = Num_of_Batches;
    All_campaign_data.(campaign_name).Summmary_of_campaign = Summmary_of_campaign;
    All_campaign_data.(campaign_name).campaign_flags = campaign_flags;
    
    % Create visualization of campaign data
    if Generate_batch_flag == 1
        % Generate plots for this campaign
        All_variables_names = fieldnames(Batch_Records.Batch_01);
        Var_all = size(All_variables_names, 1);
        
        IndPenSim_QbD_Figure_properties;
        set(0, 'DefaultFigureWindowStyle', 'docked');
        
        for plot_no = 1:1:Var_all(end)
            % Defining the figure properties
            figure;
            hold on;
            set(gcf, 'color', [1 1 1]);
            set(gca, 'color', [1 1 1]);
            grid on;
            
            Temp_Var = All_variables_names(plot_no);
            Temp_Var_name = regexprep(Temp_Var,'_',' ','emptymatch');
            Temp_title = (Temp_Var_name);
            set(gcf,'name',[char(Temp_title),' - ',campaign_name]);
            
            % Plotting all batches
            for batch_idx = Batch_start:1:Num_of_Batches
                matFileName = sprintf('Batch_%02d', batch_idx);
                Var_name_temp = char(Temp_Var);
                if isfield(Batch_Records.(matFileName).(Var_name_temp), 't') == 1
                    Temp_T_value = Batch_Records.(matFileName).(Var_name_temp).t;
                    Temp_Y_value = Batch_Records.(matFileName).(Var_name_temp).y;
                    plot(Temp_T_value,Temp_Y_value,'Color',cmap(batch_idx,:),'lineStyle','-','linewidth',Line_width_fig);
                    if diff(Temp_T_value) == 0
                        plot(Temp_T_value(1),Temp_Y_value(1),'d','MarkerSize',20,'MarkerFaceColor','r');
                    end
                    xlabel(['Time [', Batch_Records.(matFileName).(Var_name_temp).tUnit, ']'], 'fontsize', Font_size_fig1);
                    ylabel([Batch_Records.(matFileName).(Var_name_temp).name, ' [', Batch_Records.(matFileName).(Var_name_temp).yUnit, ']'], 'fontsize', Font_size_fig1);
                    title([Batch_Records.(matFileName).(Var_name_temp).name, ' plot - ', campaign_name], 'fontsize', Font_size_fig1);
                else
                    close(gcf);
                end
            end
        end
    end
    
    % Clear campaign-specific variables for next iteration
    clear Raw_Batch_data Batch_Records Summmary_of_campaign;
end

% Save the complete simulation data with all campaigns
Generate_CSV_file_flag = 1;
if Generate_CSV_file_flag == 1
    All_variables = 1;
end
fprintf('Simulation complete. Generated %d campaigns.\n', num_campaigns);

%% Copyright
% Stephen Goldrick 
% University of Manchester, Newcastle University and Perceptive Engineering and University College London
%
% All rights reserved. Copyright (c) Newcastle University, University of Manchester and Perceptive Engineering.

% Please reference: DOI: 10.1016/j.jbiotec.2014.10.029 and https://doi.org/10.1016/j.compchemeng.2019.05.037




