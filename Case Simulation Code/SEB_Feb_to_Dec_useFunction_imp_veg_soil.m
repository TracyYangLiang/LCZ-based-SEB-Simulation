%% Feb–Dec IMP LST/PET simulation for Paris (parfor, row-wise). 
clc,clear;

addpath 'D:\Paris Case Study\Case Simulation Code'; 	% Add the folder containing updated scripts/functions to the MATLAB search path
load('Model_input_Feb_to_Dec_paris.mat');	
load('SEB_imp_Jan_results_model_input_paris.mat');	

LST_result_previous_month=LST_result_summary;
PET_result_previous_month=PET_result_summary;
beta_imp=0;

row_max=size(Height,1);
column_max=size(Height,2);
LW_input1a=LW_input;
SW_input1a=SW_input;
t2m_input1a=t2m;
wind_speed_10m_input1a=wind_speed_10m;
qa_1000hpa_10m_input1a=qa_1000hpa;
P_input1a=P_input;

Height=Height;
z0m_input(z0m_input<0.001) = 0;
z0m_input=z0m_input;

parfor row_no=1:row_max
[LST_result_summary_Feb_to_Dec(:,:,row_no),PET_result_summary_Feb_to_Dec(:,:,row_no),beta_imp_summary_Feb_to_Dec(:,:,row_no),qsat_summary_Feb_to_Dec(:,:,row_no)]=one_row_LST_imp_Feb_to_Dec(row_no,column_max,...
beta_imp,LST_result_previous_month,PET_result_previous_month,...
LW_input1a,SW_input1a,t2m_input1a,...
wind_speed_10m_input1a,qa_1000hpa_10m_input1a,P_input1a,...
Height,z0m_input);
row_no
end
disp("Job compeletely done (impervious)");

%imagesc(permute(LST_result_summary_Feb_to_Dec(8042,:,:), [3 2 1]));colorbar


save('SEB_imp_results_Feb_to_Dec_pairs.mat','LST_result_summary_Feb_to_Dec','-v7.3');

%%   Feb–Dec veg LST/PET simulation for Paris (parfor, row-wise).
clc,clear;

addpath 'D:\Paris Cas Study\11'; 	% Add the folder containing updated scripts/functions to the MATLAB search path
load('Model_input_Feb_to_Dec_paris.mat');	
load('SEB_veg_Jan_results_model_input_paris.mat');	

LST_result_previous_month=LST_result_summary; % copy Jan LST obtained from April calculation
PET_result_previous_month=PET_result_summary;
beta_veg=0.1;

row_max=size(Height,1);
column_max=size(Height,2);
LW_input1a=LW_input;
SW_input1a=SW_input;
t2m_input1a=t2m;
wind_speed_10m_input1a=wind_speed_10m;
qa_1000hpa_input1a=qa_1000hpa;

P_input1a=P_input;

Height=Height;
z0m_input=z0m_input;


parfor row_no=1:row_max
[LST_result_summary_Feb_to_Dec(:,:,row_no),PET_result_summary_Feb_to_Dec(:,:,row_no),beta_veg_summary_Feb_to_Dec(:,:,row_no),qsat_summary_Feb_to_Dec(:,:,row_no)]=one_row_LST_veg_Feb_to_Dec(row_no,column_max,...
beta_veg,LST_result_previous_month,PET_result_previous_month,...
LW_input1a,SW_input1a,t2m_input1a,...
wind_speed_10m_input1a,qa_1000hpa_input1a,P_input1a,...
Height,z0m_input);
row_no
end
disp("Job compeletely done");
save('SEB_veg_results_Feb_to_Dec_Paris.mat','LST_result_summary_Feb_to_Dec','PET_result_summary_Feb_to_Dec','beta_veg_summary_Feb_to_Dec','qsat_summary_Feb_to_Dec','-v7.3');

%%  Feb–Dec SOIL LST/PET simulation for Paris (parfor, row-wise).
clc,clear;

addpath 'D:\Paris Cas Study\11'; 	% Add the folder containing updated scripts/functions to the MATLAB search path
load('Model_input_Feb_to_Dec_paris.mat');	
load('SEB_soil_Jan_results_model_input_pairs.mat');	

LST_result_previous_month=LST_result_summary; % copy Jan LST obtained from April calculation
PET_result_previous_month=PET_result_summary;
beta_soil=0;

row_max=size(Height,1);
column_max=size(Height,2);
LW_input1a=LW_input;
SW_input1a=SW_input;
t2m_input1a=t2m;
wind_speed_10m_input1a=wind_speed_10m;
qa_1000hpa_input1a=qa_1000hpa;
P_input1a=P_input;

Height=Height;
z0m_input=z0m_input;

parfor row_no=1:row_max
[LST_result_summary_Feb_to_Dec(:,:,row_no),PET_result_summary_Feb_to_Dec(:,:,row_no),beta_soil_summary_Feb_to_Dec(:,:,row_no),qsat_summary_Feb_to_Dec(:,:,row_no)]=one_row_LST_soil_Feb_to_Dec(row_no,column_max,...
beta_soil,LST_result_previous_month,PET_result_previous_month,...
LW_input1a,SW_input1a,t2m_input1a,...
wind_speed_10m_input1a,qa_1000hpa_input1a,P_input1a,...
Height,z0m_input);
row_no
end
disp("Job compeletely done");

save('SEB_soil_results_Feb_to_Dec_paris.mat','LST_result_summary_Feb_to_Dec','PET_result_summary_Feb_to_Dec','beta_soil_summary_Feb_to_Dec','qsat_summary_Feb_to_Dec','-v7.3');