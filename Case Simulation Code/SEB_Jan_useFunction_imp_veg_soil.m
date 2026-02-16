clc,clear
%%  Jan IMP LST/PET simulation for Paris (parfor, row-wise). 
addpath 'D:\Paris Case Study\Case Simulation Code'; 	% Add the folder containing updated scripts/functions to the MATLAB search path
load('Model_input_Jan_paris.mat');	
row_max=size(Height,1);
column_max=size(Height,2);
LW_input1a=LW_input1;
LW_input2a=LW_input2;
SW_input1a=SW_input1;
SW_input2a=SW_input2;
t2m_input1a=t2m1;
t2m_input2a=t2m2;
wind_speed_10m_input1a=wind_speed_10m1;
wind_speed_10m_input2a=wind_speed_10m2;
qa_1000hpa_10m_input1a=qa_1000hpa1;
qa_1000hpa_10m_input2a=qa_1000hpa2;
Height=Height;

for row_no=1:row_max
    [LST_result_summary(:,:,row_no),PET_result_summary(:,:,row_no)]=one_row_LST_imp(row_no,column_max,...
        LW_input1a,LW_input2a,SW_input1a,SW_input2a,t2m_input1a,t2m_input2a,...
        wind_speed_10m_input1a,wind_speed_10m_input2a,qa_1000hpa_10m_input1a,qa_1000hpa_10m_input2a,...
        Height,z0m_input);
    row_no
end
disp("Job compeletely done (impervious)");
a = squeeze(LST_result_summary(2,:,:));
b = squeeze(LST_result_summary(3,:,:));

save('SEB_imp_Jan_results_model_input_paris.mat','LST_result_summary','PET_result_summary');

%%  Jan veg LST/PET simulation for Paris (parfor, row-wise)
clc,clear
addpath 'D:\Paris Cas Study\CAS Simulation Code'; 	% Add the folder containing updated scripts/functions to the MATLAB search path
load('Model_input_Jan_paris.mat');
row_max=size(Height,1);
column_max=size(Height,2);
LW_input1a=LW_input1;
LW_input2a=LW_input2;
SW_input1a=SW_input1;
SW_input2a=SW_input2;
t2m_input1a=t2m1;
t2m_input2a=t2m2;
wind_speed_10m_input1a=wind_speed_10m1;
wind_speed_10m_input2a=wind_speed_10m2;
qa_1000hpa_10m_input1a=qa_1000hpa1;
qa_1000hpa_10m_input2a=qa_1000hpa2;
Height=Height;
z0m_input=z0m_input;


for row_no=1:row_max
[LST_result_summary(:,:,row_no),PET_result_summary(:,:,row_no),ra]=one_row_LST_veg(row_no,column_max,...
LW_input1a,LW_input2a,SW_input1a,SW_input2a,t2m_input1a,t2m_input2a,...
wind_speed_10m_input1a,wind_speed_10m_input2a,qa_1000hpa_10m_input1a,qa_1000hpa_10m_input2a,...
Height,z0m_input);
row_no
end
disp("Job compeletely done");
a = squeeze(LST_result_summary(2,:,:));
b = squeeze(LST_result_summary(3,:,:));

save('SEB_veg_Jan_results_model_input_paris.mat','LST_result_summary','PET_result_summary');
%%  Jan SOIL LST/PET simulation for Paris (parfor, row-wise).
clc,clear

addpath 'D:\Paris Cas Study\CAS Simulation Code'; 	% Add the folder containing updated scripts/functions to the MATLAB search path
load('Model_input_Jan_paris.mat');	
row_max=size(Height,1);
column_max=size(Height,2);
LW_input1a=LW_input1;
LW_input2a=LW_input2;
SW_input1a=SW_input1;
SW_input2a=SW_input2;
t2m_input1a=t2m1;
t2m_input2a=t2m2;
wind_speed_10m_input1a=wind_speed_10m1;
wind_speed_10m_input2a=wind_speed_10m2;
qa_1000hpa_10m_input1a=qa_1000hpa1;
qa_1000hpa_10m_input2a=qa_1000hpa2;
Height=Height;
z0m_input=z0m_input;


for row_no=1:row_max
[LST_result_summary(:,:,row_no),PET_result_summary(:,:,row_no)]=one_row_LST_soil(row_no,column_max,...
LW_input1a,LW_input2a,SW_input1a,SW_input2a,t2m_input1a,t2m_input2a,...
wind_speed_10m_input1a,wind_speed_10m_input2a,qa_1000hpa_10m_input1a,qa_1000hpa_10m_input2a,...
Height,z0m_input);
row_no
end
disp("Job completed done");
a = squeeze(LST_result_summary(2,:,:));
b = squeeze(LST_result_summary(3,:,:));

save('SEB_soil_Jan_results_model_input_pairs.mat','LST_result_summary','PET_result_summary');