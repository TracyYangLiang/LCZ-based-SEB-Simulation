function [LST_result,PET_result,beta_imp_result,qsat_result] = one_row_LST_imp_Feb_to_Dec(row_no,column_max,...
    beta_imp,LST_result_previous_month,PET_result_previous_month,...
    LW_input1a,SW_input1a,t2m_input1a,...
    wind_speed_10m_input1a,qa_1000hpa_10m_input1a,P_input1a,...
    Height,z0m_input)
beta_imp_result=repmat(beta_imp,752,column_max);

LST_initial=LST_result_previous_month(:,:,row_no); % copy Jan LST obtained from April calculation
PET_initial=PET_result_previous_month(:,:,row_no); % copy Jan PET obtained from April calculation

LW1_a=LW_input1a(row_no,1:column_max,:);
LW1=squeeze(LW1_a);
LW_input=LW1';

SW1_a=SW_input1a(row_no,1:column_max,:);
SW1=squeeze(SW1_a);
SW_input=SW1';

t2m_a=t2m_input1a(row_no,1:column_max,:);
t2m1_a2=squeeze(t2m_a);
t2m=t2m1_a2';

wind_speed_a=wind_speed_10m_input1a(row_no,1:column_max,:);
wind_speed_a2=squeeze(wind_speed_a);
wind_speed_10m=wind_speed_a2';

qa_1000hpa_a=qa_1000hpa_10m_input1a(row_no,1:column_max,:);
qa_1000hpa_a2=squeeze(qa_1000hpa_a);
qa_1000hpa=qa_1000hpa_a2';

P_a=P_input1a(row_no,1:column_max,:);
P_a2=squeeze(P_a);
P_input=P_a2';

Height=Height(row_no,1:column_max);
z0m_input=z0m_input(row_no,1:column_max);
beta_imp_result=repmat(0,752,column_max);

Stefan_Boltzmann_constant=5.67*10^(-8); % W m-2 K-4
surface_emissivity=0.95; % impervious
kv=0.4; %% Von Karman's constant
rho=1.21; % mean air density
Cp=1013; % specific heat of air at constant pressure
H_humidity_meansure=2; % Height of airT/humidity
H_wind_measure=10; % Height of wind speed measurement
patm=101.325; %% atmospheric pressure (kPa)
rho_w=1000; % water density (kg/m3)


a1_spring=0.238;
a1_summer=0.238;
a1_autumn=0.238;
a1_winter=0.238;
a3_spring=-16.7;
a3_summer=-16.7;
a3_autumn=-16.7;
a3_winter=-16.7;
a1=[repmat(a1_winter,2,column_max);repmat(a1_winter,696,column_max);repmat(a1_spring,744,column_max);repmat(a1_spring,720,column_max);repmat(a1_spring,744,column_max);repmat(a1_summer,720,column_max);repmat(a1_summer,744,column_max);repmat(a1_summer,744,column_max);repmat(a1_autumn,720,column_max);repmat(a1_autumn,744,column_max);repmat(a1_autumn,720,column_max);repmat(a1_winter,744,column_max)];
a3=[repmat(a3_winter,2,column_max);repmat(a3_winter,696,column_max);repmat(a3_spring,744,column_max);repmat(a3_spring,720,column_max);repmat(a3_spring,744,column_max);repmat(a3_summer,720,column_max);repmat(a3_summer,744,column_max);repmat(a3_summer,744,column_max);repmat(a3_autumn,720,column_max);repmat(a3_autumn,744,column_max);repmat(a3_autumn,720,column_max);repmat(a3_winter,744,column_max)];
a2=0.427; % paramter for heat storage (impervious), averaged value of paved and buildings


Ts_input=LST_initial; % Jan LST
Ts_now_2=Ts_input(751,:); % modeled LST of Jan last No.2 hour
Ts_now_1=Ts_input(752,:); % modeled LST of Jan last No.1 hour
LST_first2hours(1,:)=Ts_input(751,:);  % modeled LST of Jan last No.2 hour
LST_first2hours(2,:)=Ts_input(752,:); % modeled LST of Jan last No.1 hour
LST_result=LST_first2hours; % modeled LST of April last 2 hours

conversion_factor=3600000; % for PET calculation, convert PET unit from m/s to mm/h, 1000mm/(1/3600h)
PET_result=PET_initial; % Jan PET
PET_result(isnan(PET_result))=0;
oumiga=2.1; % model parameter for beta calculation

P_input=P_input; % Precipitation of Jan-Dec, 2011
P_input(P_input<0)=0;

beta_imp_initial=beta_imp_result;
beta_imp_result=beta_imp_initial(751:752,:); % beta for Jan last 2 hours, no precipitation, then equals to 0

if mean(mean(LST_result))>273.15
    esat_result=0.61078.*exp(17.27.*(LST_result-273.15)./((LST_result)));
else
    esat_result=0.61078.*exp(21.875.*(LST_result-273.15)./((LST_result-273.15)+265));
end
qsat_result=0.622.*esat_result./(patm-0.378.*esat_result);


AH_spring=25;
AH_summer=50;
AH_autumn=25;
AH_winter=25;
working_ratio=0.75;
non_working_ratio=1-working_ratio;
working_time_start=1; % UTC,BJ time 9am
working_time_end=10; % UTC,BJ time 18pm

Feb_days=29;
spring_days=31+30+31;
summer_days=30+31+31;
autumn_days=30+31+30;
Dec_days=31;
AH_spring_one_day=[repmat(AH_spring*non_working_ratio,(working_time_start-0),column_max);repmat(AH_spring*working_ratio,(working_time_end-working_time_start),column_max);repmat(AH_spring*non_working_ratio,(23-working_time_end+1),column_max)];
AH_summer_one_day=[repmat(AH_summer*non_working_ratio,(working_time_start-0),column_max);repmat(AH_summer*working_ratio,(working_time_end-working_time_start),column_max);repmat(AH_summer*non_working_ratio,(23-working_time_end+1),column_max)];
AH_autumn_one_day=[repmat(AH_autumn*non_working_ratio,(working_time_start-0),column_max);repmat(AH_autumn*working_ratio,(working_time_end-working_time_start),column_max);repmat(AH_autumn*non_working_ratio,(23-working_time_end+1),column_max)];
AH_winter_one_day=[repmat(AH_winter*non_working_ratio,(working_time_start-0),column_max);repmat(AH_winter*working_ratio,(working_time_end-working_time_start),column_max);repmat(AH_winter*non_working_ratio,(23-working_time_end+1),column_max)];
AH=[AH_winter_one_day(end-1:end,:);repmat(AH_spring_one_day,Feb_days,1);repmat(AH_spring_one_day,spring_days,1);repmat(AH_summer_one_day,summer_days,1);repmat(AH_autumn_one_day,autumn_days,1);repmat(AH_winter_one_day,Dec_days,1)];
%% Part 1 Net shortwave radiation
SW=SW_input(743:end,:); % Jan last 2 hours+Feb-Dec
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
Albedo=0.177; % albedo of impevious
SW_net=SW_aftercorrection.*(1-Albedo);
%% Part 2 Net longwave radiation
LW_downwards=LW_input(743:end,:); % Jan last 2 hours+Feb-Dec

LW_downwards_corrected=LW_downwards;
%% Part 3 Sensible heat
z0m=repmat(z0m_input,8042,1);
Ta=t2m(743:end,:); % in Kelvin; t2m of Jan last 2 hours+Feb-Dec
wind_speed=wind_speed_10m(743:end,:); % wind speed of Jan last 2 hours+Feb-Dec
wind_speed(wind_speed<1)=1; % wind speed correction
wind_speed(wind_speed>10)=10; % wind speed correction

Height=repmat(Height,8042,1);
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
P_sum_result=repmat(0,2,column_max);

for i=3:8042
    %% Part 4 Latent heat
    qa=qa_1000hpa(743:end,:); % specific humidity of Jan last 2 hours+Feb-Dec
    P_sum=sum(P_input(((718+i):(741+i)),:),1); % sum of precipitation during the past 24hours
    P_sum_result=[P_sum_result;P_sum];
    P=P_sum; % keep 5 digits
    [row2,col2]=size(P);
    for index=1:col2 % if accumulated P over past 24 hours >10mm, then beta_impervious=1; otherwise, beta_imp=0
        if P(1,index)>10
            beta_imp(1,index)=1;
        else
            beta_imp(1,index)=0;
        end
    end
    beta_imp_result=[beta_imp_result;beta_imp];
    %% Part 5 Ground heat
    SW_net_now=SW_net(i,:); % the current hour
    SW_net_now_1=SW_net(i-1,:); % 1 hour ago, now-1
    LW_downwards_corrected_now=LW_downwards_corrected(i,:);
    LW_downwards_corrected_now_1=LW_downwards_corrected(i-1,:);
    Qnow_1=SW_net_now_1+LW_downwards_corrected_now_1-(1-surface_emissivity).*LW_downwards_corrected_now_1-surface_emissivity.*Stefan_Boltzmann_constant.*(Ts_now_1.^4); % net radiation of 2 hour ago
    Qnow_part1=SW_net_now+LW_downwards_corrected_now-(1-surface_emissivity).*LW_downwards_corrected_now; % SW_net + LW_downwards-LW_upwards_part1 of current hour
    %% Part 6 Solve for LST
    Ta_now=Ta(i,:); % the current hour
    qa_now=qa(i,:);
    Lamda_now=Lamda(i,:);
    beta_imp_now=beta_imp_result(i,:);
    a1_now=a1(i,:);
    a3_now=a3(i,:);
    AH_now=AH(i,:);
    %%
    %Stability correction
    if i>3
        H_temp_dyn(i,:)=(LST-Ta_now)./ra_now;
        u_star(i,:)=kv.*wind_speed(i,:)./ra_part3(i,:);
        L_MO(i,:)=-u_star(i,:).^3./(kv.*9.81./Ta(i,:).*H_temp_dyn(i,:));
        zoL(i,:)= min(max(H_humidity_meansure./L_MO(i,:),-10),1);
        zoL_x_m(i,:)=(1-15.*zoL(i,:)).^0.25;
        zoL_x_h(i,:)=(1-9.*zoL(i,:)).^0.25;
        Phi_m(i,:)=...
            zoL(i,:).*(zoL(i,:)>0).*4.7+...
            (zoL(i,:)<0).*(-2.*log((1+zoL_x_m(i,:))/2)-log((1+zoL_x_m(i,:).^2)/2)+2.*atan(zoL_x_m(i,:))-pi/2);
        Phi_h(i,:)=...
            zoL(i,:).*(zoL(i,:)>0).*4.7+...
            (zoL(i,:)<0).*(-2.*log((1+zoL_x_h(i,:))/2)-log((1+zoL_x_h(i,:).^2)/2)+2.*atan(zoL_x_h(i,:))-pi/2);
        ra_part3_m(i,:)=ra_part3(i,:)-Phi_m(i,:);
        ra_part4_m(i,:)=ra_part4(i,:)-Phi_h(i,:);
        ra(i,:)=ra_part3_m(i,:).*ra_part4_m(i,:)./ra_part5(i,:);
    end
    ra_now=ra(i,:);
    ra_now(ra_now==inf)=50;
    ra_now(ra_now==-inf)=50;
    ra_now(ra_now==0)=50;


    y=@(Ts,AH_now,SW_net_now,LW_downwards_corrected_now,Ta_now,ra_now,Lamda_now,beta_imp_now,qa_now,a1_now,a3_now,Qnow_part1,Qnow_1)AH_now+SW_net_now...
        +LW_downwards_corrected_now-(1-surface_emissivity).*LW_downwards_corrected_now-surface_emissivity.*Stefan_Boltzmann_constant.*(Ts.^4)...
        -rho*Cp.*(Ts-Ta_now)./ra_now-Lamda_now.*beta_imp_now.*rho*(0.622*0.6108*exp(17.27*(Ts-273.15)./((Ts-273.15)+237.3))./(patm-0.378*0.6108*exp(17.27*(Ts-273.15)./((Ts-273.15)+237.3)))-qa_now)./ra_now...
        -(a1_now.*(Qnow_part1-surface_emissivity.*Stefan_Boltzmann_constant.*(Ts.^4))+a2*((Qnow_part1-surface_emissivity.*Stefan_Boltzmann_constant.*(Ts.^4))-Qnow_1)+a3_now);


    LST=zeros(1,column_max);

    for j=1:column_max
        disp(['Current values: LST_result(' num2str(i-1), ',', num2str(j), ',' ,num2str(row_no), ') = ', num2str(LST_result(i-1,j))]);
        LST(1,j)=fzero(@(Ts)y(Ts,AH_now(1,j),SW_net_now(1,j),LW_downwards_corrected_now(1,j),Ta_now(1,j),ra_now(1,j),Lamda_now(1,j),beta_imp_now(1,j),qa_now(1,j),a1_now(1,j),a3_now(1,j),Qnow_part1(1,j),Qnow_1(1,j)),273); % fzero initial value is assigned the previous hour's LST
    
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