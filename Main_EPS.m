
% Power generation/consumption simulation
% Body -Z축 속도방향 정렬, Body +Y축 지구 반대방향 지향
% Body +Y축에 GPS 안테나 부착함 (차이 없음)

clear all; close all; clc;

%% Constants

DtoR = pi/180 ;
RtoD = 180/pi ;

mu = 398601.2 * (10^3)^3 ;  % m^3/sec^2
R_Earth = 6378.145 * 10^3 ;     % m, radius of the Earth
R_Sun = 1.392 * 10^6 * 10^3 / 2 ; % m, radius of the Sun
from_Sun_to_Earth = 1.496 * 10^8 * 10^3 ;   % m, distance between the Sun and the Earth (1AU)

MJD_VernalEq = julianDate( [ 2015, 3, 21, 12, 0, 0 ] ) - 2400000.5 ;   % modified Julian date
MJD = julianDate( [ 2015, 9, 21, 12, 0, 0 ] ) - 2400000.5 ;   % modified Julian date

%% Orbit elements

a = 200 * 10^3 + R_Earth ;  % 고도 200km, EOL 가정
e = 0.00 ;     % eccentricity
Inclination = 98 * DtoR ;
Right_Ascension = 0 * DtoR ;
Argument_of_Perigee = 0 * DtoR ;

%% (1) S/C position (ECI frame)

Period = 2 * pi * sqrt(a^3/mu) ;

Eccentric_Anomaly = 0 ;
True_Anomaly = 0 ;
       
r = a * (1-e^2) / ( 1 + e*cos(DtoR*True_Anomaly) ) ;

Xdot_0_Orbital = - sqrt(mu*a)/r * sin(Eccentric_Anomaly) ;
Ydot_0_Orbital = sqrt(mu*a*(1-e^2))/r * cos(Eccentric_Anomaly) ;
X_0_Orbital = a * ( cos(Eccentric_Anomaly) - e ) ;
Y_0_Orbital = a * sqrt(1-e^2) * sin(Eccentric_Anomaly) ;

dt = 1 ;
t = 0 : dt : (5*Period) ;   % for 5 periods

X_ECI = zeros(6,length(t)) ;    % ECI coordinate
C = fn_DCM_coord(3,Argument_of_Perigee) * fn_DCM_coord(1,Inclination) * fn_DCM_coord(3,Right_Ascension) ;
X_ECI(1:3,1) = C' * [ X_0_Orbital; Y_0_Orbital; 0 ] ;
X_ECI(4:6,1) = C' * [ Xdot_0_Orbital; Ydot_0_Orbital; 0 ] ;

for i = 2 : length(t)
    X_ECI(:,i) = RK4(@Xdot,X_ECI(:,i-1),dt);
end

%% (2) Earth position w.r.t. the Sun (HAE frame)

mu_Sun = 132712440018 * (10^3)^3 ;
MJD_AOP = julianDate( [ 2014, 12, 22, 12, 0, 0 ] ) - 2400000.5 ;   % modified Julian date

a_E = from_Sun_to_Earth ;
e_E = 0.016 ;

Earth_Orbital = zeros(3,length(t)) ;
Earth_HAE = zeros(3,length(t));

for i = 1 : length(t)   
    Mean_Anomaly = sqrt(mu_Sun/a_E^3) * ( MJD - MJD_AOP ) * 24 * 60 * 60 ;
    Eccentric_Anomaly = Mean_Anomaly + e_E * sin(Mean_Anomaly) ;   % initial guess
    
    del_Eccentric_Anomaly = 1 ; 
    
    % Newton-Raphson Method
    while ( del_Eccentric_Anomaly > 10^(-10) )
        del_Mean_Anomaly = Eccentric_Anomaly - e_E * sin(Eccentric_Anomaly) - Mean_Anomaly ;
        del_Eccentric_Anomaly = - del_Mean_Anomaly / ( 1 - e_E * cos(Eccentric_Anomaly) ) ;
        Eccentric_Anomaly = Eccentric_Anomaly + del_Eccentric_Anomaly ;
    end
    
    Eccentric_Anomaly = mod(Eccentric_Anomaly, 2*pi) ;
    temp = sqrt((1+e_E)/(1-e_E)) * tan(Eccentric_Anomaly/2) ;
    True_Anomaly = 2 * atan(temp) ;
    
    Earth_HAE(1:3,i) = fn_DCM_point(3,True_Anomaly)*[ 0; from_Sun_to_Earth; 0 ] ;
end

%% (3) S/C position (HAE frame)

X_HAE = zeros(3,length(t)) ;
for i = 1 : length(t)
    X_HAE(1:3,i) = Earth_HAE(:,i) + fn_DCM_coord(1,DtoR*23.5)*X_ECI(1:3,i) ;
end

%% (4) Check eclipse condition

from_Earth_to_Apex = R_Earth * from_Sun_to_Earth / ( R_Sun - R_Earth ) ;  % m
rho_c = asin( ( R_Sun - R_Earth ) / from_Sun_to_Earth ) ; % rad

rho_s = zeros(1,length(t)) ;
rho_p = zeros(1,length(t)) ;
theta = zeros(1,length(t)) ;

Vec_from_Sat_to_Sun = zeros(3,length(t)) ;
Vec_from_Sat_to_Earth = zeros(3,length(t)) ;
Check_Eclipse = zeros(1,length(t)) ;

for i = 1 : length(t)    
    Vec_from_Sat_to_Sun(:,i) = - X_HAE(:,i) ;
    Vec_from_Sat_to_Earth(:,i) = Earth_HAE(:,i) - X_HAE(:,i) ;
    rho_s(i) = asin( R_Sun / norm(Vec_from_Sat_to_Sun(:,i)) ) ;  % rad
    rho_p(i) = asin( R_Earth / norm(Vec_from_Sat_to_Earth(:,i)) ) ;    % rad
    theta(i) = acos( dot(Vec_from_Sat_to_Sun(:,i),Vec_from_Sat_to_Earth(:,i)) / norm(Vec_from_Sat_to_Sun(:,i)) / norm(Vec_from_Sat_to_Earth(:,i)) ) ;  % rad
    
    % 2. Total eclipse
    if ( ( from_Sun_to_Earth < norm(Vec_from_Sat_to_Sun(:,i)) ) && ( norm(Vec_from_Sat_to_Sun(:,i)) < from_Sun_to_Earth + from_Earth_to_Apex ) && ( rho_p(i) - rho_s(i) > theta(i) ) )
        Check_Eclipse(i) = 2 ;
    % 4. No eclipse
    else 
        Check_Eclipse(i) = 4 ;
    end    
end

%% (5) Attitude of each panel (in HAE frame)

Att_x_plus = zeros(3,length(t)) ;
Att_x_minus = zeros(3,length(t)) ;
Att_y_plus = zeros(3,length(t)) ;
Att_y_minus = zeros(3,length(t)) ;
Att_z_plus = zeros(3,length(t)) ;
Att_z_minus = zeros(3,length(t)) ;

for i = 1 : length(t)        
    % C2 : from ECI to Body
    % body Y축 위치벡터 정렬
    Velocity_unit = X_ECI(4:6,i)/norm(X_ECI(4:6,i)) ;
    Position_unit = X_ECI(1:3,i)/norm(X_ECI(1:3,i)) ;
    C2 = fn_TRIAD( [ Velocity_unit, Position_unit ], [ 0, 0, -1; 0, 1, 0 ]' ) ;  % body -Z축 속도방향 정렬,  Body Y축 지구 반대방향 지향

    % C : From Body to HAE (Body -> ECI -> HAE)
    C = fn_DCM(1,DtoR*23.5) * C2' ; % 수정함 (2013. 11. 21)
    
    Att_x_plus(:,i) = C * [ 1; 0; 0 ] ;
    Att_x_minus(:,i) = C * [ -1; 0; 0 ] ;
    Att_y_plus(:,i) = C * [ 0; 1; 0 ] ;
    Att_y_minus(:,i) = C * [ 0; -1; 0 ] ;
    Att_z_plus(:,i) = C * [ 0; 0; 1 ] ;
    Att_z_minus(:,i) = C * [ 0; 0; -1 ];    
end

%% (6) Calculate power generation

SolarIntensity = 1300 ; % W/m^2
SolarPanelEfficiency = 0.243 ;  % EOL 가정

% % 태양전지판 면적 (1) 2U 2개 날개처럼 전개
% Area_x_plus = 0.011; 
% Area_x_minus = 0.011;
% Area_y_plus = 0.011 + 0.011*2 - 0.005/2 - 0.001; % 지구방향, 지구반대방향에 카메라 부착, GPS 안테나 부착
% Area_y_minus = 0.011 + 0.011*2 - 0.005/2; % 지구방향, 지구반대방향에 카메라 부착
% Area_z_plus = 0; 
% Area_z_minus = 0;
% % 여기까지

% 태양전지판 면적 (2) 1U 4개 뒤에서 전개
Area_x_plus = 0.011 + 0.005/sqrt(2) ; 
Area_x_minus = 0.011 + 0.005/sqrt(2) ;
Area_y_plus = 0.011 + 0.005/sqrt(2) - 0.005/2 ; % 지구방향, 지구반대방향에 카메라 부착, GPS 안테나 부착
Area_y_minus = 0.011 + 0.005/sqrt(2) - 0.005/2; % 지구방향, 지구반대방향에 카메라 부착
Area_z_plus = 0.005*4/sqrt(2) ; 
Area_z_minus = 0.005*4/sqrt(2) ;
% 여기까지

co_elevation_x_plus = zeros(1,length(t)) ;
co_elevation_x_minus = zeros(1,length(t)) ;
co_elevation_y_plus = zeros(1,length(t)) ;
co_elevation_y_minus = zeros(1,length(t)) ;
co_elevation_z_plus = zeros(1,length(t)) ;
co_elevation_z_minus = zeros(1,length(t)) ;

Power_x_plus = zeros(1,length(t)) ;
Power_x_minus = zeros(1,length(t)) ;
Power_y_plus = zeros(1,length(t)) ;
Power_y_minus = zeros(1,length(t)) ;
Power_z_plus = zeros(1,length(t)) ;
Power_z_minus = zeros(1,length(t)) ;

for i = 1 : length(t)    
    if ( Check_Eclipse(i) == 4 )
        co_elevation_x_plus(i) = acos( dot( Vec_from_Sat_to_Sun(:,i), Att_x_plus(:,i) ) / norm(Vec_from_Sat_to_Sun(:,i)) ) ;
        if ( co_elevation_x_plus(i) < pi/2 )
            Power_x_plus(i) = SolarIntensity*SolarPanelEfficiency*Area_x_plus * cos(co_elevation_x_plus(i)) ;
        end
        
        co_elevation_x_minus(i) = acos( dot( Vec_from_Sat_to_Sun(:,i), Att_x_minus(:,i) ) / norm(Vec_from_Sat_to_Sun(:,i)) ) ;
        if ( co_elevation_x_minus(i) < pi/2 )
            Power_x_minus(i) = SolarIntensity*SolarPanelEfficiency*Area_x_minus * cos(co_elevation_x_minus(i)) ;
        end
        
        co_elevation_y_plus(i) = acos( dot( Vec_from_Sat_to_Sun(:,i), Att_y_plus(:,i) ) / norm(Vec_from_Sat_to_Sun(:,i)) ) ;
        if ( co_elevation_y_plus(i) < pi/2 )
            Power_y_plus(i) = SolarIntensity*SolarPanelEfficiency*Area_y_plus * cos(co_elevation_y_plus(i)) ;
        end
        
        co_elevation_y_minus(i) = acos( dot( Vec_from_Sat_to_Sun(:,i), Att_y_minus(:,i) ) / norm(Vec_from_Sat_to_Sun(:,i)) ) ;
        if ( co_elevation_y_minus(i) < pi/2 )
            Power_y_minus(i) = SolarIntensity*SolarPanelEfficiency*Area_y_minus * cos(co_elevation_y_minus(i)) ;
        end

        co_elevation_z_plus(i) = acos( dot( Vec_from_Sat_to_Sun(:,i), Att_z_plus(:,i) ) / norm(Vec_from_Sat_to_Sun(:,i)) ) ;
        if ( co_elevation_z_plus(i) < pi/2 )
            Power_z_plus(i) = SolarIntensity*SolarPanelEfficiency*Area_z_plus * cos(co_elevation_z_plus(i)) ;
        end
        
        co_elevation_z_minus(i) = acos( dot( Vec_from_Sat_to_Sun(:,i), Att_z_minus(:,i) ) / norm(Vec_from_Sat_to_Sun(:,i)) ) ;      
        if ( co_elevation_z_minus(i) < pi/2 )
            Power_z_minus(i) = SolarIntensity*SolarPanelEfficiency*Area_z_minus * cos(co_elevation_z_minus(i)) ;
        end
    end    
end

%% (7) Power budget

BatterySize = 10.0000 ; % Wh
AttType = 'Normal';
% AttType = 'Safe';

PowerGeneration = zeros(1,length(t)) ;
PowerConsumption = zeros(1,length(t)) ;
BatteryRemained = zeros(1,length(t)) ;
DepthOfDischarge = zeros(1,length(t)) ;

for i = 1 : length(t)
    PowerGeneration(i) = Power_x_plus(i) + Power_x_minus(i) + Power_y_plus(i) + Power_y_minus(i) + Power_z_plus(i) + Power_z_minus(i) ;
    % BCR efficiency 고려 (85%)
    PowerGeneration(i) = 0.85*PowerGeneration(i) ;
  
    if ( strcmp(AttType,'Safe') )
        PowerConsumption(i) = 1.19954 ;    % Safe mode (2013.11.22)
    else
%         PowerConsumption(i) = 1.94954*(85/90) + 3.26732*(5/90) ; % GPS 10%, 5분 통신
        PowerConsumption(i) = 1.94954*(80/90) + 3.26732*(10/90) ; % GPS 10%, 10분 통신
%         PowerConsumption(i) = 2.34954*(80/90) + 3.26732*(10/90) ; % GPS 50%, 10분 통신
%         PowerConsumption(i) = 2.84954*(80/90) + 3.26732*(10/90) ; % GPS 100%, 10분 통신
    end
    
    PowerConsumption(i) = PowerConsumption(i) / 0.95 ;  % Efficiency 고려 (95%)
    
    if ( i == 1 )
        BatteryRemained(i) = BatterySize ;    % Fully charged
    else
        BatteryRemained(i) = BatteryRemained(i-1) + ( PowerGeneration(i-1) - PowerConsumption(i-1) )*(dt/3600) ;
        if ( BatteryRemained(i) > BatterySize )
            BatteryRemained(i) = BatterySize ;
        end
    end
    
    DepthOfDischarge(i) = 100 * (BatterySize - BatteryRemained(i)) / BatterySize ;
end

%% (8) Power generation sum (Wh)
PowerGenerationSum = 0 ;

for i = 1 : length(t)
    PowerGenerationSum = PowerGenerationSum + PowerGeneration(i)*dt/3600 ;
end

PowerGenerationSum % Wh

%% Graphs

figure(1)
plot3( 10^(-3)*Earth_HAE(1,:), 10^(-3)*Earth_HAE(2,:), 10^(-3)*Earth_HAE(3,:), 'rx' )
hold on
plot3( 10^(-3)*X_HAE(1,:), 10^(-3)*X_HAE(2,:), 10^(-3)*X_HAE(3,:), 'b' )
hold off
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
axis equal
title('S/C 3D Orbit in HAE Frame')
grid on
xlim( [-7000+10^(-3)*Earth_HAE(1,1) 7000+10^(-3)*Earth_HAE(1,1)] )
ylim( [-7000+10^(-3)*Earth_HAE(2,1) 7000+10^(-3)*Earth_HAE(2,1)] )
zlim( [-7000+10^(-3)*Earth_HAE(3,1) 7000+10^(-3)*Earth_HAE(3,1)] )

figure(2)
plot3( 10^(-3)*X_ECI(1,:), 10^(-3)*X_ECI(2,:), 10^(-3)*X_ECI(3,:), 'b' )
hold off
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
axis equal
title('S/C 3D Orbit in ECI Frame')
grid on
xlim( [-7000 7000] )
ylim( [-7000 7000] )
zlim( [-7000 7000] )

figure(3)
subplot(3,1,1)
plot( 10^(-3)*X_HAE(2,:), 10^(-3)*X_HAE(3,:), 'b' )
axis equal
xlabel('Y (km)')
ylabel('Z (km)')
title('S/C orbit projection(YZ) in HAE frame')
subplot(3,1,2)
plot( 10^(-3)*X_HAE(1,:), 10^(-3)*X_HAE(2,:), 'b' )
axis equal
xlabel('X (km)')
ylabel('Y (km)')
title('S/C orbit projection(XY) in HAE frame')
subplot(3,1,3)
plot( 10^(-3)*X_HAE(1,:), 10^(-3)*X_HAE(3,:), 'b' )
axis equal
xlabel('X (km)')
ylabel('Z (km)')
title('S/C orbit projection(XZ) in HAE frame')

figure(5)
plot( t, -10^(-3)*Earth_HAE(1,:), 'r' )
hold on
plot( t, -10^(-3)*Earth_HAE(2,:), 'g' )
hold on
plot( t, -10^(-3)*Earth_HAE(3,:), 'b' )
hold off
title('Sun vector XYZ in ECI frame')
legend( 'x', 'y', 'z')

figure(4)
plot( t/3600, Check_Eclipse, 'b*' )
xlabel('Time (hr)')
ylabel('Eclipse condition')
title('1: Partial Eclipse, 2: Total Eclipse, 3: Annular Eclipse, 4: No Eclipse')

figure(5)
plot( t/3600, RtoD*co_elevation_x_plus, 'ko' )
hold on
plot( t/3600, RtoD*co_elevation_x_minus, 'mx' )
hold on
plot( t/3600, RtoD*co_elevation_y_plus, 'b+' )
hold on
plot( t/3600, RtoD*co_elevation_y_minus, 'r*' )
hold on
plot( t/3600, RtoD*co_elevation_z_plus, 'gs' )
hold on
plot( t/3600, RtoD*co_elevation_z_minus, 'cd' )
hold off
xlabel('Time (hr)')
ylabel('Co-elevation angle (deg)')
legend( 'x plus', 'x minus', 'y plus', 'y minus', 'z plus', 'z minus')

figure(6)
plot( t/3600, Power_x_plus, 'ko' )
hold on
plot( t/3600, Power_x_minus, 'mx' )
hold on
plot( t/3600, Power_y_plus, 'b+' )
hold on
plot( t/3600, Power_y_minus, 'r*' )
hold on
plot( t/3600, Power_z_plus, 'gs' )
hold on
plot( t/3600, Power_z_minus, 'cd' )
hold off
xlabel('Time (hr)')
ylabel('Power generation (W)')
legend( 'x plus', 'x minus', 'y plus', 'y minus', 'z plus', 'z minus')

figure(7)
plot( t/3600, PowerGeneration, 'b-' )
hold on
plot( t/3600, PowerConsumption, 'r-' )
hold off
xlabel('Time (hr)')
ylabel('Power (W)')
legend( 'Power generation (W)', 'Power consumption (W)' )

figure(8)
subplot(2,1,1)
plot( t/3600, BatteryRemained )
xlabel('Time (hr)')
ylabel('Battery remained (Wh)')
subplot(2,1,2)
plot( t/3600, DepthOfDischarge, 'b' )
hold on
plot( t/3600, 20, 'r')
hold off
xlabel('Time (hr)')
ylabel('DoD (%)')
legend( 'DoD history', 'recommended' )



