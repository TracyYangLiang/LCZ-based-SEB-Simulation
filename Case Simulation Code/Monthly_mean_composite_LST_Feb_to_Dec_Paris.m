clc,clear;
%% 
addpath 'D:\Paris Case Study\Case Simulation Code'; 	% Add the folder containing updated scripts/functions to the MATLAB search path
load('Pv_Ps_Pimp_paris.mat');
%% load LST_Pv_Ps_Pimp

load('SEB_imp_results_Feb_to_Dec_pairs.mat','LST_result_summary_Feb_to_Dec');

LST_imp=permute(LST_result_summary_Feb_to_Dec,[3 2 1]);

load('SEB_veg_results_Feb_to_Dec_Paris.mat','LST_result_summary_Feb_to_Dec');

LST_veg=permute(LST_result_summary_Feb_to_Dec,[3 2 1]);

load('SEB_soil_results_Feb_to_Dec_paris.mat','LST_result_summary_Feb_to_Dec');

LST_soil=permute(LST_result_summary_Feb_to_Dec,[3 2 1]);
%%

Pimp_all_year = zeros(21,21,8042);
h=[2, 29*24,31*24,30*24,31*24,30*24,31*24,31*24,30*24,31*24,30*24,31*24];
t = 1;
for i = 1:length(h)
    hs = h(i);
    Pimp_all_year(:,:,t:t+hs-1) = repmat(Pi2(:,:,i), [1, 1, hs]);
    t = t + hs;
end
Pveg_all_year = zeros(21,21,8042);
t = 1;
for i = 1:length(h)
    hs = h(i);
    Pveg_all_year(:,:,t:t+hs-1) = repmat(Pv2(:,:,i), [1, 1, hs]);
    t = t + hs;
end
Psoil_all_year = zeros(21,21,8042);
t = 1;
for i = 1:length(h)
    hs = h(i);
    Psoil_all_year(:,:,t:t+hs-1) = repmat(Ps2(:,:,i), [1, 1, hs]);
    t = t + hs;
end

LST_composite=LST_imp.*Pimp_all_year+LST_veg.*Pveg_all_year+LST_soil.*Psoil_all_year;

    %%

start_index = 3;
LST_0000_all = [];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          nnbhnnbLST_0000_all = [];
LST_0100_all = [];
LST_0200_all = [];
LST_0300_all = [];
LST_0400_all = [];
LST_0500_all = [];
LST_0600_all = [];
LST_0700_all = [];
LST_0800_all = [];
LST_0900_all = [];
LST_1000_all = [];
LST_1100_all = [];
LST_1200_all = [];
LST_1300_all = [];
LST_1400_all = [];
LST_1500_all = [];
LST_1600_all = [];
LST_1700_all = [];
LST_1800_all = [];
LST_1900_all = [];
LST_2000_all = [];
LST_2100_all = [];
LST_2200_all = [];
LST_2300_all = [];
LST_2400_all = [];


for k = 2:length(h)
    end_index = start_index + h(k) - 1;
    current_month = LST_composite(:,:,start_index:end_index);
    n = h(k)/24;
    nn(k) = n;

   LST_0100 = [];
    for day=1:n
        LST_0100(:,:,day)=current_month(:,:,1+24*(day-1));
    end
    LST_0100_all = cat(3,LST_0100_all,LST_0100);
    LST_0100_final(:,:,k) = mean(LST_0100,3); 

    LST_0200 = [];
    for day=1:n
        LST_0200(:,:,day)=current_month(:,:,2+24*(day-1));
    end
    LST_0200_all = cat(3,LST_0200_all,LST_0200);
    LST_0200_final(:,:,k) = mean(LST_0200,3);    

    LST_0300 = [];
    for day=1:n
        LST_0300(:,:,day)=current_month(:,:,3+24*(day-1));
    end
    
    LST_0300_all = cat(3,LST_0300_all,LST_0300);
    LST_0300_final(:,:,k) = mean(LST_0300,3);

    LST_0400 = [];
    for day=1:n
        LST_0400(:,:,day)=current_month(:,:,4+24*(day-1));
    end
    
    LST_0400_all = cat(3,LST_0400_all,LST_0400);
    LST_0400_final(:,:,k) = mean(LST_0400,3); 

    LST_0500 = [];
    for day=1:n
    LST_0500(:,:,day)=current_month(:,:,5+24*(day-1));
    end
    
    LST_0500_all = cat(3,LST_0500_all,LST_0500);
    LST_0500_final(:,:,k) = mean(LST_0500,3);

    LST_0600 = [];
    for day=1:n
        LST_0600(:,:,day)=current_month(:,:,6+24*(day-1));
    end
    
    LST_0600_all = cat(3,LST_0600_all,LST_0600);
    LST_0600_final(:,:,k) = mean(LST_0600,3);

    LST_0700 = [];
    for day=1:n
    LST_0700(:,:,day)=current_month(:,:,7+24*(day-1));
    end
    
    LST_0700_all = cat(3,LST_0700_all,LST_0700);
    LST_0700_final(:,:,k) = mean(LST_0700,3);

    LST_0800 = [];
    for day=1:n
    LST_0800(:,:,day)=current_month(:,:,8+24*(day-1));
    end
    
    LST_0800_all = cat(3,LST_0800_all,LST_0800);
    LST_0800_final(:,:,k) = mean(LST_0800,3);

    LST_0900 = [];
    for day=1:n
    LST_0900(:,:,day)=current_month(:,:,9+24*(day-1));
    end
    
    LST_0900_all = cat(3,LST_0900_all,LST_0900);
    LST_0900_final(:,:,k) = mean(LST_0900,3);

    LST_1000 = [];
    for day=1:n
    LST_1000(:,:,day)=current_month(:,:,10+24*(day-1));
    end
    
    LST_1000_all = cat(3,LST_1000_all,LST_1000);
    LST_1000_final(:,:,k) = mean(LST_1000,3);

    LST_1100 = [];
    for day=1:n
    LST_1100(:,:,day)=current_month(:,:,11+24*(day-1));
    end
    
    LST_1100_all = cat(3,LST_1100_all,LST_1100);
    LST_1100_final(:,:,k) = mean(LST_1100,3);

    LST_1200 = [];
    for day=1:n
    LST_1200(:,:,day)=current_month(:,:,12+24*(day-1));
    end
    
    LST_1200_all = cat(3,LST_1200_all,LST_1200);
    LST_1200_final(:,:,k) = mean(LST_1200,3);


    LST_1300 = [];
    for day=1:n
    LST_1300(:,:,day)=current_month(:,:,13+24*(day-1));
    end
    
    LST_1300_all = cat(3,LST_1300_all,LST_1300);
    LST_1300_final(:,:,k) = mean(LST_1300,3);


    LST_1400 = [];
    for day=1:n
    LST_1400(:,:,day)=current_month(:,:,14+24*(day-1));
    end
    
    LST_1400_all = cat(3,LST_1400_all,LST_1400);
    LST_1400_final(:,:,k) = mean(LST_1400,3);

    LST_1500 = [];
    for day=1:n
    LST_1500(:,:,day)=current_month(:,:,15+24*(day-1));
    end
    
    LST_1500_all = cat(3,LST_1500_all,LST_1500);
    LST_1500_final(:,:,k) = mean(LST_1500,3);

    LST_1600 = [];
    for day=1:n
    LST_1600(:,:,day)=current_month(:,:,16+24*(day-1));
    end
    
    LST_1600_all = cat(3,LST_1600_all,LST_1600);
    LST_1600_final(:,:,k) = mean(LST_1600,3);

    LST_1700 = [];
    for day=1:n
    LST_1700(:,:,day)=current_month(:,:,17+24*(day-1));
    end
    
    LST_1700_all = cat(3,LST_1700_all,LST_1700);
    LST_1700_final(:,:,k) = mean(LST_1700,3);

    LST_1800 = [];
    for day=1:n
    LST_1800(:,:,day)=current_month(:,:,18+24*(day-1));
    end
    
    LST_1800_all = cat(3,LST_1800_all,LST_1800);
    LST_1800_final(:,:,k) = mean(LST_1800,3);
    
    LST_1900 = [];
    for day=1:n
    LST_1900(:,:,day)=current_month(:,:,19+24*(day-1));
    end
    
    LST_1900_all = cat(3,LST_1900_all,LST_1900);
    LST_1900_final(:,:,k) = mean(LST_1900,3);

    LST_2000 = [];
    for day=1:n
    LST_2000(:,:,day)=current_month(:,:,20+24*(day-1));
    end
    
    LST_2000_all = cat(3,LST_2000_all,LST_2000);
    LST_2000_final(:,:,k) = mean(LST_2000,3);

    LST_2100 = [];
    for day=1:n
    LST_2100(:,:,day)=current_month(:,:,21+24*(day-1));
    end
    
    LST_2100_all = cat(3,LST_2100_all,LST_2100);
    LST_2100_final(:,:,k) = mean(LST_2100,3);


    LST_2200 = [];
    for day=1:n
    LST_2200(:,:,day)=current_month(:,:,22+24*(day-1));
    end
    
    LST_2200_all = cat(3,LST_2200_all,LST_2200);
    LST_2200_final(:,:,k) = mean(LST_2200,3);

    LST_2300 = [];
    for day=1:n
    LST_2300(:,:,day)=current_month(:,:,23+24*(day-1));
    end
    
    LST_2300_all = cat(3,LST_2300_all,LST_2300);
    LST_2300_final(:,:,k) = mean(LST_2300,3);
    
    LST_2400 = [];
    for day=1:n
    LST_2400(:,:,day)=current_month(:,:,24*(day));
    end
    
    LST_2400_all = cat(3,LST_2400_all,LST_2400);
    LST_2400_final(:,:,k) = mean(LST_2400,3);

    start_index = end_index + 1;
end
   save('LST_composite.mat','LST_0800_final','LST_0900_final','LST_1000_final','LST_1100_final','LST_1200_final','LST_1300_final','LST_1400_final','LST_1500_final','LST_1600_final','LST_1700_final','LST_1800_final','LST_1900_final','LST_2000_final','LST_2100_final','LST_2200_final','LST_2300_final','LST_2400_final');
  