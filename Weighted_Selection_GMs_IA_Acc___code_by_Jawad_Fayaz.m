clear all; clc; close all; fclose all; direc = pwd;
%% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author : JAWAD FAYAZ (email: j.fayaz@ucl.ac.uk)
% visit: (https://jfayaz.github.io)
%
% ------------------------------ Instructions ------------------------------------- 
% This code selects ground motions based on weighted multi-IM normalized objective 
% function using different intensity measures including Sa, PGV, CAV, Ia, and D5-95
%
% Citation: 
%         % Jawad Fayaz, Miguel Medalla, Pablo Torres, and Carmine Galasso (2023). 
%         % "Recurrent Neural Networks based Generalized Ground Motion Model for Chilean Subduction Seismic Environment". 
%         % Structural Safety, Vol. 100, Pages 102282
%
%         % Jawad Fayaz and Carmine Galasso (2022). 
%         % "A Generalized Ground Motion Model for Consistent Mainshock-Aftershock Intensity Measures using Successive Recurrent Neural Networks". 
%         % Bulletin of Earthquake Engineering, Vol. 20, Pages 6467-6486
%     
% INPUTS:
% This code uses the target IMs provided in 'Target_IMs.xlsx' file to
% select best matching GMs using the IMs provided in 'IM_Data.mat' based on
% on weighted multi-IM normalized objective function. The example files are
% provided with the code.
% The 'IM_Data.mat' file which must be in the current folder and contain 2 variables
%      'Spectra'       --> n x 1  Cell structure containing the spectra of 'n' GMs
%      'Other_IMs'     --> (n+1) x 4  Cell structure containing the 4 IMs of 'n' GMs with 4 titles in the same order as the example file
% The 'Target_IMs.xlsx' file must contain 3 columns in the same format as
% the example file. The column 1 can be empty, the code uses the columns 2
% and 3 of the excel
%
% OUTPUT:
% The outputs will be provided in 'SEL_IMs_GMs.mat' file which contain the scale and unscaled versions of the selected IMs. 
% The output file will contain
%     'Sel_CAV_Unscaled','Sel_Spectra_Unscaled','Sel_Spectra_Scaled','Sel_PGV_Unscaled','Sel_Ia_Unscaled','Sel_CAV_Scaled','Sel_PGV_Scaled','Sel_Ia_Scaled','Sel_D595','Sel_Error','Sel_GMs_ids'
% The variables are self-explanatory 
% 'Sel_GMs_ids' --> n x 1  array containing indices of the selected GM IMs
% with the same order as the 'IM_Data.mat' file
%%-------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ======================== USER INPUTS =============================== %%
%%
% Number of required GMs
No_of_GMs = 10;

% Primary spectral period
T         = 0.8;

% Allowable scaling factors
Scale_Factors = [0.5:0.1:2]; 

% Tolerance level for mismatch of Duration (values between: 0 and 0.99)
D595_Tol = 0.25;

% Sa spectral period range to match between the target and GMs spectra 
T1       = 0.5*T;
T2       = 2*T;
T_range  = [ [T1:0.05:T-0.05], T, [T+0.05:0.05:T2] ];

% Weights for different IMs (though the weights are normalized in the code,
% try to keep the sum of the weights equal to 1) (values between: 0 and 1)
Wt_Sa   = 0.6;
Wt_Ia   = 0.1;
Wt_PGV  = 0.1;
Wt_CAV  = 0.1;

% Plot finally selected IMs  (options: 'Yes' or 'No')
Plot_IMs  = 'Yes';

%%%%%%================= END OF USER INPUT ========================%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =============== Prelim Computations ================
load('IM_Data.mat')
target   = xlsread('Target_IMs.xlsx');
tar_CAV  = target(1,2);
tar_PGV  = target(2,2);
tar_D595 = target(3,2);
tar_Ia   = target(4,2);
tar_Sa   = target(5:end,:);

GMs_D595 = cell2mat(Other_IMs(2:end,1));
GMs_Ia   = cell2mat(Other_IMs(2:end,2));
GMs_PGV  = cell2mat(Other_IMs(2:end,3));
GMs_CAV  = cell2mat(Other_IMs(2:end,4));
tmp      = horzcat(Spectra{:});
GMs_Sa   = tmp(:,[2:2:end])';
Periods  = tmp(:,1);

r_D595       = find( GMs_D595 > (1-D595_Tol)*tar_D595  &  GMs_D595 < (1+D595_Tol)*tar_D595 );
SelGMs_D595  = GMs_D595(r_D595,:);
SelGMs_Ia    = GMs_Ia(r_D595,:);
SelGMs_PGV   = GMs_PGV(r_D595,:);
SelGMs_CAV   = GMs_CAV(r_D595,:);
SelGMs_Sa    = GMs_Sa(r_D595,:);

[Sel_Spectra_Unscaled,Sel_Spectra_Scaled,Sel_CAV_Unscaled,Sel_PGV_Unscaled,Sel_Ia_Unscaled,Sel_Sa_Unscaled,Sel_CAV_Scaled,Sel_PGV_Scaled,Sel_Ia_Scaled,Sel_Sa_Scaled,Sel_D595,Sel_Error,Sel_GMs_ids] = IM_match(No_of_GMs,tar_CAV,tar_PGV,tar_Ia,tar_Sa,SelGMs_D595,SelGMs_Ia,SelGMs_PGV,SelGMs_CAV,SelGMs_Sa,Wt_Sa,Wt_Ia,Wt_PGV,Wt_CAV,T_range,Periods,Scale_Factors); 
Sel_GMs_ids  = r_D595(Sel_GMs_ids);

save('SEL_IMs_GMs.mat','Sel_CAV_Unscaled','Sel_Spectra_Unscaled','Sel_Spectra_Scaled','Sel_PGV_Unscaled','Sel_Ia_Unscaled','Sel_CAV_Scaled','Sel_PGV_Scaled','Sel_Ia_Scaled','Sel_D595','Sel_Error','Sel_GMs_ids','-v7.3')

%% =============== Plotting ================
if strcmpi(Plot_IMs,'yes')==1
    figure()
    subplot(3,2,1);hold on
    plot([T1 T1],[0.001 5],'--','linewidth',2,'Color','b')
    plot([T2 T2],[0.001 5],'--','linewidth',2,'Color','b')
    plot([T T],[0.001 5],'-','linewidth',3,'Color','b')
    plot(repmat(Periods',No_of_GMs,1)',Sel_Sa_Unscaled','linewidth',1.5,'Color',[0.3 0.3 0.3],'DisplayName','Unscaled')
    plot(tar_Sa(:,1),tar_Sa(:,2),'linewidth',4,'Color','r','DisplayName','Target')
    set(gca,'fontsize',16,'FontName', 'Times New Roman','LineWidth', 1.25,'TickDir','out','TickLength', [0.01 0.01])
    ylabel ('Sa (g)','fontsize',18)
    xlabel('T (s)','fontsize',18)
    set(gca,'Yscale','log')
    title('Unscaled Spectra','fontsize',18)

    subplot(3,2,2);hold on
    plot([T1 T1],[0.001 5],'--','linewidth',2,'Color','b')
    plot([T2 T2],[0.001 5],'--','linewidth',2,'Color','b')
    plot([T T],[0.001 5],'-','linewidth',3,'Color','b')
    plot(repmat(Periods',No_of_GMs,1)',Sel_Sa_Scaled','linewidth',1.5,'Color',[0.3 0.3 0.3],'DisplayName','Scaled')
    plot(tar_Sa(:,1),tar_Sa(:,2),'linewidth',4,'Color','r','DisplayName','Target')
    set(gca,'fontsize',16,'FontName', 'Times New Roman','LineWidth', 1.25,'TickDir','out','TickLength', [0.01 0.01])
    ylabel ('Sa (g)','fontsize',18)
    xlabel('T (s)','fontsize',18)
    set(gca,'Yscale','log')
    title('Scaled Spectra','fontsize',18)

    subplot(3,2,3);hold on
    histogram(Sel_D595,'DisplayName','GMs')
    plot([tar_D595,tar_D595],[0,No_of_GMs/2],'linewidth',6,'Color','r','DisplayName','Target')
    set(gca,'fontsize',16,'FontName', 'Times New Roman','LineWidth', 1.25,'TickDir','out','TickLength', [0.01 0.01])
    ylabel ('Histogram','fontsize',18)
    xlabel('D5-95','fontsize',18)
    legend('fontsize',14)

    subplot(3,2,4);hold on
    histogram(Sel_Ia_Unscaled,'DisplayName','Unscaled')
    histogram(Sel_Ia_Scaled,'DisplayName','Scaled')
    plot([tar_Ia,tar_Ia],[0,No_of_GMs/2],'linewidth',6,'Color','r','DisplayName','Target')
    set(gca,'fontsize',16,'FontName', 'Times New Roman','LineWidth', 1.25,'TickDir','out','TickLength', [0.01 0.01])
    ylabel ('Histogram','fontsize',18)
    xlabel('Ia','fontsize',18)
    legend('fontsize',14)

    subplot(3,2,5);hold on
    histogram(Sel_CAV_Unscaled,'DisplayName','Unscaled')
    histogram(Sel_CAV_Scaled,'DisplayName','Scaled')
    plot([tar_CAV,tar_CAV],[0,No_of_GMs/2],'linewidth',6,'Color','r','DisplayName','Target')
    set(gca,'fontsize',16,'FontName', 'Times New Roman','LineWidth', 1.25,'TickDir','out','TickLength', [0.01 0.01])
    ylabel ('Histogram','fontsize',18)
    xlabel('CAV','fontsize',18)
    legend('fontsize',14)

    subplot(3,2,6);hold on
    histogram(Sel_PGV_Unscaled,'DisplayName','Unscaled')
    histogram(Sel_PGV_Scaled,'DisplayName','Scaled')
    plot([tar_PGV,tar_PGV],[0,No_of_GMs/2],'linewidth',6,'Color','r','DisplayName','Target')
    set(gca,'fontsize',16,'FontName', 'Times New Roman','LineWidth', 1.25,'TickDir','out','TickLength', [0.01 0.01])
    ylabel ('Histogram','fontsize',18)
    xlabel('PGV','fontsize',18)
    legend('fontsize',14)
end

%% =============== Computing Wieghted Error and Selecting GMs ================
function [Sel_Spectra_Unscaled,Sel_Spectra_Scaled,Sel_CAV_Unscaled,Sel_PGV_Unscaled,Sel_Ia_Unscaled,Sel_Sa_Unscaled,Sel_CAV_Scaled,Sel_PGV_Scaled,Sel_Ia_Scaled,Sel_Sa_Scaled,Sel_D595,Sel_Error,Sel_GMs_ids] = IM_match(No_of_GMs,tar_CAV,tar_PGV,tar_Ia,tar_Sa,SelGMs_D595,SelGMs_Ia,SelGMs_PGV,SelGMs_CAV,SelGMs_Sa,Wt_Sa,Wt_Ia,Wt_PGV,Wt_CAV,T_range,Periods,Scale_Factors) 
    SFs = Scale_Factors';
    c1  = find(abs(Periods-T_range(1)) == min(abs(Periods-T_range(1))));
    c2  = find(abs(Periods-T_range(end)) == min(abs(Periods-T_range(end))));
    r1  = find(abs(tar_Sa(:,1)-T_range(1)) == min(abs(tar_Sa(:,1)-T_range(1))));
    r2  = find(abs(tar_Sa(:,1)-T_range(end)) == min(abs(tar_Sa(:,1)-T_range(end))));

    GMsSa_int   = [];
    for i = 1:size(SelGMs_Sa,1)
        GMsSa_int(i,:) = interp1(Periods',SelGMs_Sa(i,:),T_range);   
    end
        
    for sf = 1:length(SFs)
        GMsCAV    = SelGMs_CAV.*SFs(sf);
        GMsPGV    = SelGMs_PGV.*SFs(sf);
        GMsIa     = SelGMs_Ia.*(SFs(sf)^2);
        GMsSaint  = GMsSa_int.*(SFs(sf)^2);
       
        tarSa_int    = interp1(tar_Sa(:,1),tar_Sa(:,2),T_range);        
       
        IA_CAV(:,sf) = compute_acc(tar_CAV,GMsCAV);
        IA_PGV(:,sf) = compute_acc(tar_PGV,GMsPGV);
        IA_Ia(:,sf)  = compute_acc(tar_Ia,GMsIa);
        IA_Sa(:,sf)  = compute_IA(tarSa_int,GMsSaint);
              
        error(:,sf) =  -1 .* (Wt_Sa.*IA_Sa(:,sf) + Wt_Ia.*IA_Ia(:,sf) + Wt_PGV.*IA_PGV(:,sf) + Wt_CAV.*IA_CAV(:,sf))./(Wt_Sa+Wt_Ia+Wt_CAV+Wt_PGV);   %% Weighted Error
    end
    
    [err_min, err_sf]        = min(error');            %% Minimum error for each GM among all SFs
    [~,idxs_sorted]          = sort(err_min');         %% Sorting the minimum errors for the GMs
    err_sf_sorted            = err_sf(idxs_sorted')';
    
    
    Sel_D595                = SelGMs_D595(idxs_sorted(1:No_of_GMs),:);
            
    Sel_CAV_Unscaled        = SelGMs_CAV(idxs_sorted(1:No_of_GMs),:);
    Sel_PGV_Unscaled        = SelGMs_PGV(idxs_sorted(1:No_of_GMs),:);
    Sel_Ia_Unscaled         = SelGMs_Ia(idxs_sorted(1:No_of_GMs),:);
    Sel_Sa_Unscaled         = SelGMs_Sa(idxs_sorted(1:No_of_GMs),:);

    Sel_CAV_Scaled          = SelGMs_CAV(idxs_sorted(1:No_of_GMs),:).*SFs(err_sf_sorted(1:No_of_GMs));
    Sel_PGV_Scaled          = SelGMs_PGV(idxs_sorted(1:No_of_GMs),:).*SFs(err_sf_sorted(1:No_of_GMs));
    Sel_Ia_Scaled           = SelGMs_Ia(idxs_sorted(1:No_of_GMs),:).*(SFs(err_sf_sorted(1:No_of_GMs)).^2);
    Sel_Sa_Scaled           = SelGMs_Sa(idxs_sorted(1:No_of_GMs),:).*SFs(err_sf_sorted(1:No_of_GMs));   

    Sel_SFs                 = SFs(err_sf_sorted(1:No_of_GMs));
    Sel_GMs_ids             = idxs_sorted(1:No_of_GMs);
    
    for i = 1:No_of_GMs
        Sel_Error(i,1)       = error(idxs_sorted(i,1),err_sf_sorted(i,1));
        Sel_Spectra_Unscaled{i,1} = [Periods, Sel_Sa_Unscaled(i,:)'];
        Sel_Spectra_Scaled{i,1}   = [Periods, Sel_Sa_Scaled(i,:)'];
    end
    
end

function acc = compute_acc(X,Xhat)
    acc = 1 - abs((X - Xhat)./(X+Xhat));
end

function IA = compute_IA(X,Xhat)
    IA  = 1 - sum((Xhat-X).^2,2)./sum((abs(Xhat-mean(X))+abs(X-mean(X))).^2,2);
end

