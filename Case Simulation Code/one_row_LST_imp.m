function [LST_result,PET_result] = one_row_LST_imp(row_no,column_max,...
    LW_input1a,LW_input2a,SW_input1a,SW_input2a,t2m_input1a,t2m_input2a,...
    wind_speed_10m_input1a,wind_speed_10m_input2a,qa_1000hpa_10m_input1a,qa_1000hpa_10m_input2a,...
    Height,z0m_input)
LW1_a=LW_input1a(row_no,1:column_max,:);
LW1=squeeze(LW1_a);
LW_input1=LW1';
LW2_a=LW_input2a(row_no,1:column_max,:);
LW2=squeeze(LW2_a);
LW_input2=LW2';

SW1_a=SW_input1a(row_no,1:column_max,:);
SW1=squeeze(SW1_a);
SW_input1=SW1';
SW2_a=SW_input2a(row_no,1:column_max,:);
SW2=squeeze(SW2_a);
SW_input2=SW2';

t2m_a=t2m_input1a(row_no,1:column_max,:);
t2m1_a2=squeeze(t2m_a);
t2m1=t2m1_a2';
t2m_b=t2m_input2a(row_no,1:column_max,:);
t2m_b2=squeeze(t2m_b);
t2m2=t2m_b2';

wind_speed_a=wind_speed_10m_input1a(row_no,1:column_max,:);
wind_speed_a2=squeeze(wind_speed_a);
wind_speed_10m1=wind_speed_a2';
wind_speed_b=wind_speed_10m_input2a(row_no,1:column_max,:);
wind_speed_b2=squeeze(wind_speed_b);
wind_speed_10m2=wind_speed_b2';

qa_1000hpa_a=qa_1000hpa_10m_input1a(row_no,1:column_max,:);
qa_1000hpa_a2=squeeze(qa_1000hpa_a);
qa_1000hpa1=qa_1000hpa_a2';
qa_1000hpa_b=qa_1000hpa_10m_input2a(row_no,1:column_max,:);
qa_1000hpa_b2=squeeze(qa_1000hpa_b);
qa_1000hpa2=qa_1000hpa_b2';

Height=Height(row_no,1:column_max);
z0m_input=z0m_input(row_no,1:column_max);
beta_imp_result=repmat(0,752,column_max);

Stefan_Boltzmann_constant=5.67*10^(-8); % W m-2 K-4
surface_emissivity=0.9; % impervious
kv=0.4; %% Von Karman's constant
rho=1.21; % mean air density
Cp=1013; % specific heat of air at constant pressure
H_humidity_meansure=2; % Height of airT/humidity
H_wind_measure=10; % Height of wind speed measurement
patm=101.325; %% atmospheric pressure (kPa)
rho_w=1000; % water density (kg/m3)
a1=0.238; % paramter for heat storage (impervious),averaged value of paved and buildings
a2=0.427; % paramter for heat storage (impervious), averaged value of paved and buildings
a3=-16.7; % paramter for heat storage (impervious),averaged value of paved and buildings

t2m=cat(1,t2m1,t2m2);
Ta=t2m(737:738,:); % start at 0am Jan 2011 local time
Ts_now_2=Ta(1,:);
Ts_now_1=Ta(2,:);
LST_first2hours(1,:)=Ta(1,:); % LST of the first hour, Ta 0:00am local time
LST_first2hours(2,:)=Ta(2,:); % LST of the second hour, Ta 1:00am local time
LST_result=LST_first2hours; % LST of the first 2 hours
conversion_factor=3600000; % for PET calculation, convert PET unit from m/s to mm/h, 1000mm/(1/3600h)
PET_first2hours(1,:)=repmat(0,1,column_max); % PET of the first hour (0am local time), assigned 0
PET_first2hours(2,:)=repmat(0,1,column_max); % PET of the second hour (1am local time), assigned 0
PET_result=PET_first2hours; % PET of the first 2 hours

if mean(mean(LST_result))>273.15
    esat_result=0.61078.*exp(17.27.*(LST_result-273.15)./((LST_result)));
else
    esat_result=0.61078.*exp(21.875.*(LST_result-273.15)./((LST_result-273.15)+265));
end
qsat_result=0.622.*esat_result./(patm-0.378.*esat_result);
%% Part 1 Net shortwave radiation
SW_input=cat(1,SW_input1,SW_input2);
SW=SW_input(737:1488,:); % start at 0am Jan 2011 local time
[row,column]=size(SW);
SW_aftercorrection=zeros(row,column);

for Row=1:row
    for Col=1:column
        if SW(Row,Col)>0
            SW_aftercorrection(Row,Col)=SW(Row,Col);
        else
            SW_aftercorrection(Row,Col)=0;
        end
    end
end
Albedo=0.177; % albedo of impevious177
SW_net=SW_aftercorrection.*(1-Albedo);
%% Part 2 Net longwave radiation
LW_input=cat(1,LW_input1,LW_input2);
LW_downwards=LW_input(737:1488,:); % start at 0am Jan 2011 local time
LW_downwards_corrected=LW_downwards;
%% Part 3 Sensible heat
z0m=repmat(z0m_input,752,1);
Ta=t2m(737:1488,:); % start at 0am Jan 2011 local time
wind_speed_10m=cat(1,wind_speed_10m1,wind_speed_10m2);
wind_speed=wind_speed_10m(737:1488,:); % start at 0am Jan 2011 local time
wind_speed(wind_speed<1)=1; % wind speed correction
wind_speed(wind_speed>10)=10; % wind speed correction

Height=repmat(Height,752,1);
d=0; %% zero plane displacement height
kinematic_viscosity_air=5*10^(-11)*Ta.^2+7*10^(-8)*Ta-10^(-5); % retrived from statistics relationship, Ta in Kevin
friction_velocity=kv.*wind_speed./log(10./z0m);
Re=z0m.*friction_velocity./kinematic_viscosity_air; % Reynolds number
z0h=z0m.*exp((-kv)*(0.6/5*sqrt(Re)-3.5));
z0h(z0h<0.001)=0.001;

ra_part1=(H_wind_measure+z0m-d)./z0m;
ra_part2=(H_humidity_meansure+z0h-d)./z0h;
ra_part3=log(ra_part1);
ra_part4=log(ra_part2);
ra_part5=kv.^2.*wind_speed;
ra=ra_part3.*ra_part4./ra_part5;

ra(ra==inf)=50;
ra(ra==-inf)=50;
ra(ra==0)=50;

Lamda=1000*(2501.3-2.361*(Ta-273.15)); % Latent heat vaporization/condensation [J/kg]
qa_1000hpa=cat(1,qa_1000hpa1,qa_1000hpa2);
qa=qa_1000hpa(737:1488,:); % start at 0am Jan 2011 local time
beta_imp=0; % No precipitation in Jan, then beta_imp=0
%%
for i=3:752
    %% Part 5 Ground heat
    SW_net_now=SW_net(i,:); % the current hour
    SW_net_now_1=SW_net(i-1,:); % 1 hour ago, now-1
    LW_downwards_corrected_now=LW_downwards_corrected(i,:);
    LW_downwards_corrected_now_1=LW_downwards_corrected(i-1,:);
    Qnow_1=SW_net_now_1+LW_downwards_corrected_now_1-(1-surface_emissivity).*LW_downwards_corrected_now_1-surface_emissivity.*Stefan_Boltzmann_constant.*(Ts_now_1.^4); % net radiation of 2 hour ago
    Qnow_part1=SW_net_now+LW_downwards_corrected_now-(1-surface_emissivity).*LW_downwards_corrected_now; % SW_net + LW_downwards-LW_upwards_part1 of current hour
    %% Part 6 Solve for LST
    Ta_now=Ta(i,:); % the current hour
    ra_now=ra(i,:);
    qa_now=qa(i,:);
    Lamda_now=Lamda(i,:);

    y=@(Ts,SW_net_now,LW_downwards_corrected_now,Ta_now,ra_now,Lamda_now,qa_now,Qnow_part1,Qnow_1)SW_net_now...
        +LW_downwards_corrected_now-(1-surface_emissivity).*LW_downwards_corrected_now-surface_emissivity.*Stefan_Boltzmann_constant.*(Ts.^4)...
        -rho*Cp.*(Ts-Ta_now)./ra_now-Lamda_now*beta_imp.*rho*(0.622*0.6108*exp(17.27*(Ts-273.15)./((Ts-273.15)+237.3))./(patm-0.378*0.6108*exp(17.27*(Ts-273.15)./((Ts-273.15)+237.3)))-qa_now)./ra_now...
        -(a1*(Qnow_part1-surface_emissivity.*Stefan_Boltzmann_constant.*(Ts.^4))+a2*((Qnow_part1-surface_emissivity.*Stefan_Boltzmann_constant.*(Ts.^4))-Qnow_1)+a3);

    LST=zeros(1,column_max);
    for j=1:column_max
        LST(1,j)=fzero(@(Ts)y(Ts,SW_net_now(1,j),LW_downwards_corrected_now(1,j),Ta_now(1,j),ra_now(1,j),Lamda_now(1,j),qa_now(1,j),Qnow_part1(1,j),Qnow_1(1,j)),LST_result(i-1,j)); % fzero initial value is assigned the previous hour's LST
    end

    if mean(mean(LST))>273.15
        esat_now=0.61078.*exp(17.27.*(LST-273.15)./((LST)));
    else
        esat_now=0.61078.*exp(21.875.*(LST-273.15)./((LST-273.15)+265));
    end
    qsat_now=0.622.*esat_now./(patm-0.378.*esat_now);
    qsat_result=[qsat_result;qsat_now(:,:)];
    PET=conversion_factor.*rho.*(qsat_now-qa_now)./ra_now./rho_w; % PET of current hour
    PET(PET<0)=0;
    PET_result=[PET_result;PET(:,:)];
    LST_result=[LST_result;LST(:,:)];
    Ts_now_1=LST_result(i,:); % LST of 1 hour ago
end
end