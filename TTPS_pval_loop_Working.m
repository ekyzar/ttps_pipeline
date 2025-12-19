% This script must be run after the counts of CFOS, TDT, Dapi, and Overlap
% (cells double-labeled for CFOS and TDT) have already been loaded using
% the Batch_Load_Data_BrainJ_Both_Hemi script.
%
% The script below is the analysis of 5 brains in an 'AA' group
% compared to a 'Homecage-Homecage' group

num_samples = 300;

%% AA group

A1_o = A1c.Dapi_total;
A2_o = A2c.Dapi_total;
A3_o = A3c.Dapi_total;
A4_o = A4c.Dapi_total;
A5_o = A5c.Dapi_total;
Dapi_total_1 = ([A1_o A2_o A3_o A4_o A5_o]);

A1_o = A1c.CFOS;
A2_o = A2c.CFOS;
A3_o = A3c.CFOS;
A4_o = A4c.CFOS;
A5_o = A5c.CFOS;
CFOS_1 = ([A1_o A2_o A3_o A4_o A5_o]);

A1_o = A1c.TDT;
A2_o = A2c.TDT;
A3_o = A3c.TDT;
A4_o = A4c.TDT;
A5_o = A5c.TDT;
TDT_1 = ([A1_o A2_o A3_o A4_o A5_o]);

A1_o = A1c.Overlap;
A2_o = A2c.Overlap;
A3_o = A3c.Overlap;
A4_o = A4c.Overlap;
A5_o = A5c.Overlap;
Overlap_1 = ([A1_o A2_o A3_o A4_o A5_o]);
%% Homecage-Homecage

H1_o = H1c.Dapi_total;
H2_o = H2c.Dapi_total;
H3_o = H3c.Dapi_total;
H4_o = H4c.Dapi_total;
Dapi_total_2 = ([H1_o H2_o H3_o H4_o]);

H1_o = H1c.CFOS;
H2_o = H2c.CFOS;
H3_o = H3c.CFOS;
H4_o = H4c.CFOS;
CFOS_2 = ([H1_o H2_o H3_o H4_o]);

H1_o = H1c.TDT;
H2_o = H2c.TDT;
H3_o = H3c.TDT;
H4_o = H4c.TDT;
TDT_2 = ([H1_o H2_o H3_o H4_o]);

H1_o = H1c.Overlap;
H2_o = H2c.Overlap;
H3_o = H3c.Overlap;
H4_o = H4c.Overlap;
Overlap_2 = ([H1_o H2_o H3_o H4_o]);

columns1 = size(CFOS_1,2);
columns2 = size(CFOS_2,2);

p_values = nan(length(CFOS_1),columns1*columns2);
column = 1;

columns1 = size(CFOS_1,2);
columns2 = size(CFOS_2,2);
Areas = length(CFOS_1);


p_values = nan(length(CFOS_1),columns1*columns2);
pseudo_p = nan(length(CFOS_1),1);
column = 1;

figure;
hold on
value_kdiff = [];

for kk = 1:Areas
    
    A.Dapi_total = Dapi_total_1(kk,:);
    A.CFOS = CFOS_1(kk,:);
    A.TDT = TDT_1(kk,:);
    A.Overlap = Overlap_1(kk,:);
   
    
    B.Dapi_total = Dapi_total_2(kk,:);
    B.CFOS = CFOS_2(kk,:);
    B.TDT = TDT_2(kk,:);
    B.Overlap = Overlap_2(kk,:);
    
   
    
    %%
    
    Samples_Area1 = nan(num_samples^2,columns1);
    k_Area1 = nan(columns1,1);
    for i = 1:columns1;
        if A.Dapi_total(i)~=0 && A.Dapi_total(i)<150000
            k1 = floor(A.Overlap(i));
            k_Area1(i) = k1;
            
            N = floor(A.Dapi_total(i));
            k = floor(A.Overlap(i));
            K = floor(A.CFOS(i));
            n = floor(A.TDT(i));
            
            a = K-k;
            b = n-k;
            d = N-n;
            c = d-a;
            
            cont_t = [k b; a c];
            cont_t = floor(cont_t);
            if ~(cont_t<0)
                
                %[h,p,stats] = fishertest(cont_t);
                [h1,p1,stats] = fishertest(cont_t,'Tail', 'right','Alpha', 0.05);
                %
                % Sampling from the hypergeometric distribution
                samples_1 = hygernd(N, K, n, num_samples);
            end
            Samples_Area1(:,i) = samples_1(:);
        end
    end
    
    Samples_Area2 = nan(num_samples^2,columns2);
    k_Area2 = nan(columns2,1);
    for i = 1:columns2;
        if B.Dapi_total(i)~=0 && B.Dapi_total(i)<150000
            
            k2 = floor(B.Overlap(i));
            k_Area2(i) = k2;
            %        if k1>k2
            
            % Parameters
            N = floor(B.Dapi_total(i));
            k = floor(B.Overlap(i));
            K = floor(B.CFOS(i));
            n = floor(B.TDT(i));
            
            
            k2 = k;
            
            a = K-k;
            b = n-k;
            d = N-n;
            c = d-a;
            
            cont_t = [k b; a c];
            cont_t = floor(cont_t);
            if ~(cont_t<0)
                
                [h,p,stats] = fishertest(cont_t);
                
                
                % Sampling from the hypergeometric distribution
                samples_2 = hygernd(N, K, n, num_samples);
            end
            Samples_Area2(:,i) = samples_2(:);
        end
    end
    
    for i = 1:columns1
        samples_1 = Samples_Area1(:,i);
        k1 = k_Area1(i);
        if ~isnan(k1)
            for j = 1:columns2
                samples_2 = Samples_Area2(:,j);
                k2 = k_Area2(j);
                if ~isnan(k2)
                    % Making a new distribution by subtracting random sample 2 from random
                    % sample 1.
                    samples_new = (samples_1-samples_2);
                    
                    
                    % Subtracting the true observed sample overlap (k) from dataset 1 and 2.
                    kdiff = k1-k2;
                    
                    % %                     %Plotting the new sample distribution in a histogram
                    %                                              % Creates a new figure window
                    %                                             histogram(samples_new, 'BinMethod', 'integers', 'FaceColor', 'blue');
                    %                                             xlabel('Sample Value');
                    %                                             ylabel('Frequency');
                    %                                             title('Histogram of Samples from Hypergeometric Distribution');
                    %                                             grid on; % Adds a grid to the plot for better readability
                    %                                             hold on
                    %                                             plot([kdiff kdiff],[0 700],'r')
                    
                    % Calculating the p-value for the observed overlap given the new
                    % distribution and the new overlap, kdiff:
                    
                    % Flatten the matrix to a single vector for easier processing
                    flattened_samples = samples_new(:);
                    
                    %                     % Counting how many samples in flattened_samples are as extreme as or more extreme than kdiff
                    %                     extreme_count = sum(abs(flattened_samples) >= abs(kdiff));
                    
                    % RIGHT SIDED ONLY Counting how many samples in flattened_samples are as extreme as or more extreme than kdiff
                    % For a 'Lack-of-overlap' or LEFT-SIDED analysis, the
                    % sign below would be <=
                    extreme_count = sum((flattened_samples) >= (kdiff));
                    
                    % Calculating the p-value, using the total number of values in the flattened matrix
                    total_samples = numel(flattened_samples);
                    p_value_kdiff = extreme_count / total_samples;
                    value_kdiff = [value_kdiff kdiff];
                    p_values(kk,column) = p_value_kdiff;
                    % Displaying the p-value
                    % disp(['P-value for kdiff = ' num2str(p_value_kdiff)]);
                    column = column + 1;
                end
            end
        end
    end
    
    % Now calculate the Monte-Carlo on the Pseudo-Populations:
    N = floor(nanmean(A.Dapi_total));
    k = floor(nanmean(A.Overlap));
    K = floor(nanmean(A.CFOS));
    n = floor(nanmean(A.TDT));
    k1=k;
    
    if ~isnan(k1)
        
        a = K-k;
        b = n-k;
        d = N-n;
        c = d-a;
        
        cont_t = [k b; a c];
        cont_t = floor(cont_t);
        if ~(cont_t<0)
            
            %[h,p,stats] = fishertest(cont_t);
            [h1,p1,stats] = fishertest(cont_t,'Tail', 'right','Alpha', 0.05);
            %
            % Sampling from the hypergeometric distribution
            samples_1 = hygernd(N, K, n, num_samples);
        end
        
        N = floor(nanmean(B.Dapi_total));
        k = floor(nanmean(B.Overlap));
        K = floor(nanmean(B.CFOS));
        n = floor(nanmean(B.TDT));
        k2 = k;
        
        if ~isnan(k2)
            a = K-k;
            b = n-k;
            d = N-n;
            c = d-a;
            
            cont_t = [k b; a c];
            cont_t = floor(cont_t);
            if ~(cont_t<0)
                
                [h,p,stats] = fishertest(cont_t);
                
                samples_2 = hygernd(N, K, n, num_samples);
            end
            
            samples_new = (samples_1-samples_2);
            
            
            % Subtracting the true observed sample overlap (k) from dataset 1 and 2.
            kdiff = k1-k2;
            
            % Flatten the matrix to a single vector for easier processing
            flattened_samples = samples_new(:);
            
            %                     % Counting how many samples in flattened_samples are as extreme as or more extreme than kdiff
            %                     extreme_count = sum(abs(flattened_samples) >= abs(kdiff));
            
            % RIGHT SIDED ONLY Counting how many samples in flattened_samples are as extreme as or more extreme than kdiff
            extreme_count = sum((flattened_samples) >= (kdiff));
            
            % Calculating the p-value, using the total number of values in the flattened matrix
            total_samples = numel(flattened_samples);
            p_value_kdiff = extreme_count / total_samples;
            pseudo_p(kk) = p_value_kdiff;
        end
    end
    
    column = 1;
    disp(['Brain Area: ',num2str(kk)]);
end

p_values2 = p_values;
p_values2(isnan(p_values2))=1;
p_values2(p_values2==0)=(1/90000);
IndPvalues= sum((p_values==0)')./size(p_values,2);
IndPvalues = 1- IndPvalues;
[harmonic_p2 combined_p2 stouff_p2 exp_p2] = comb_p_val(p_values2);
Ranked = find(harmonic_p2 < 0.05/620); % change this based on multiple comparisons in your data
Ranked_p = harmonic_p2(Ranked)
[sortedPValues, originalIndices] = sort(Ranked_p, 'ascend');
Ranked = Ranked(originalIndices);

cd('F:\TRAP_Fast_ReFed_Brains\Processed_V5_2023\Evan');
save('EL_HH_pval_loop_min_thresh');

Brain_Area = A1c.Brain_Area;
Brain_Acr = A1c.Brain_Acr;
Brain_ID = A1c.Brain_ID;
Harmonic = harmonic_p2';
Fischer = combined_p2';
Stouffer = stouff_p2';
Exp = exp_p2';
DataTable_Subtraction = table(Brain_Area, Brain_Acr, Brain_ID, Harmonic, Fischer, Stouffer, Exp);
writetable(DataTable_Subtraction, 'AA_HH_Subtraction.csv');



