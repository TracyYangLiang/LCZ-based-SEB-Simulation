clc,clear;

addpath 'D:\Paris Case Study\Case Simulation Code'; 	% Add the folder containing updated scripts/functions to the MATLAB search path
load('Pv_Ps_Pimp_paris.mat');

%% load LST_Pv_Ps_Pimp

load('SEB_imp_Jan_results_model_input_paris.mat','LST_result_summary');
LST_imp=permute(LST_result_summary,[3 2 1]);
load('SEB_veg_Jan_results_model_input_paris.mat','LST_result_summary');
LST_veg=permute(LST_result_summary,[3 2 1]);
load('SEB_soil_Jan_results_model_input_pairs.mat','LST_result_summary');
LST_soil=permute(LST_result_summary,[3 2 1]);

h=[752];

Pimp_Jan = zeros(21,21,752);
Pimp_Jan(:,:,:) = repmat(Pi2(:,:,1), [1, 1, h]);

Pveg_Jan = zeros(21,21,752);
Pveg_Jan(:,:,:) = repmat(Pv2(:,:,1), [1, 1, h]);

Psoil_Jan = zeros(21,21,752);
Psoil_Jan(:,:,:) = repmat(Ps2(:,:,1), [1, 1, h]);

LST_composite=LST_imp.*Pimp_Jan+LST_veg.*Pveg_Jan+LST_soil.*Psoil_Jan;

n = 31;
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

current_month = LST_composite(:,:,9:752);
LST_0800 = [];

    for day=1:n
        LST_0800(:,:,day)=current_month(:,:,8+24*(day-1));
    end
    LST_0800_all = cat(3,LST_0800_all,LST_0800);
    LST_0800_final(:,:,1) = mean(LST_0800,3); 


    LST_0900 = [];
    for day=1:n
        LST_0900(:,:,day)=current_month(:,:,9+24*(day-1));
    end
    LST_0900_all = cat(3,LST_0900_all,LST_0900);
    LST_0900_final(:,:,1) = mean(LST_0900,3);    

    LST_1000 = [];
    for day=1:n
        LST_1000(:,:,day)=current_month(:,:,10+24*(day-1));
    end
    
    LST_1000_all = cat(3,LST_1000_all,LST_1000);
    LST_1000_final(:,:,1) = mean(LST_1000,3);



    LST_1100 = [];
    for day=1:n
        LST_1100(:,:,day)=current_month(:,:,11+24*(day-1));
    end
    
    LST_1100_all = cat(3,LST_1100_all,LST_1100);
    LST_1100_final(:,:,1) = mean(LST_1100,3); 

    LST_1200 = [];
    for day=1:n
    LST_1200(:,:,day)=current_month(:,:,12+24*(day-1));
    end
    
    LST_1200_all = cat(3,LST_1200_all,LST_1200);
    LST_1200_final(:,:,1) = mean(LST_1200,3);

    LST_1300 = [];
    for day=1:n
        LST_1300(:,:,day)=current_month(:,:,13+24*(day-1));
    end
    
    LST_1300_all = cat(3,LST_1300_all,LST_1300);
    LST_1300_final(:,:,1) = mean(LST_1300,3);

    LST_1400 = [];
    for day=1:n
    LST_1400(:,:,day)=current_month(:,:,14+24*(day-1));
    end
    
    LST_1400_all = cat(3,LST_1400_all,LST_1400);
    LST_1400_final(:,:,1) = mean(LST_1400,3);

    LST_1500 = [];
    for day=1:n
    LST_1500(:,:,day)=current_month(:,:,15+24*(day-1));
    end
    
    LST_1500_all = cat(3,LST_1500_all,LST_1500);
    LST_1500_final(:,:,1) = mean(LST_1500,3);

    LST_1600 = [];
    for day=1:n
    LST_1600(:,:,day)=current_month(:,:,16+24*(day-1));
    end
    
    LST_1600_all = cat(3,LST_1600_all,LST_1600);
    LST_1600_final(:,:,1) = mean(LST_1600,3);

    LST_1700 = [];
    for day=1:n
    LST_1700(:,:,day)=current_month(:,:,17+24*(day-1));
    end
    
    LST_1700_all = cat(3,LST_1700_all,LST_1700);
    LST_1700_final(:,:,1) = mean(LST_1700,3);

    LST_1800 = [];
    for day=1:n
    LST_1800(:,:,day)=current_month(:,:,18+24*(day-1));
    end
    
    LST_1800_all = cat(3,LST_1800_all,LST_1800);
    LST_1800_final(:,:,1) = mean(LST_1800,3);

    LST_1900 = [];
    for day=1:n
    LST_1900(:,:,day)=current_month(:,:,19+24*(day-1));
    end
    
    LST_1900_all = cat(3,LST_1900_all,LST_1900);
    LST_1900_final(:,:,1) = mean(LST_1900,3);


    LST_2000 = [];
    for day=1:n
    LST_2000(:,:,day)=current_month(:,:,20+24*(day-1));
    end
    
    LST_2000_all = cat(3,LST_2000_all,LST_2000);
    LST_2000_final(:,:,1) = mean(LST_2000,3);


    LST_2100 = [];
    for day=1:n
    LST_2100(:,:,day)=current_month(:,:,21+24*(day-1));
    end
    
    LST_2100_all = cat(3,LST_2100_all,LST_2100);
    LST_2100_final(:,:,1) = mean(LST_2100,3);

    LST_2200 = [];
    for day=1:n
    LST_2200(:,:,day)=current_month(:,:,22+24*(day-1));
    end
    
    LST_2200_all = cat(3,LST_2200_all,LST_2200);
    LST_2200_final(:,:,1) = mean(LST_2200,3);

    LST_2300 = [];
    for day=1:n
    LST_2300(:,:,day)=current_month(:,:,23+24*(day-1));
    end
    
    LST_2300_all = cat(3,LST_2300_all,LST_2300);
    LST_2300_final(:,:,1) = mean(LST_2300,3);

    LST_2400 = [];
    for day=1:n
    LST_2400(:,:,day)=current_month(:,:,24+24*(day-1));
    end
    
    LST_2400_all = cat(3,LST_2400_all,LST_2400);
    LST_2400_final(:,:,1) = mean(LST_2400,3);

save('LST_composite_jan.mat','LST_0800_final','LST_0900_final','LST_1000_final','LST_1100_final','LST_1200_final','LST_1300_final','LST_1400_final','LST_1500_final','LST_1600_final','LST_1700_final','LST_1800_final','LST_1900_final','LST_2000_final','LST_2100_final','LST_2200_final','LST_2300_final','LST_2400_final');

