%% Question 5 Harikrishnan R N, 18CHE147 
% Calculation of Terminal Velocity 
%% Constants : 
% Defining the constants used : 
g = 9.80665;  % m/s^2, acceleration due to gravity 
rhop = 1800;   % kg/m^3 , Particle Density 
rhof = 994.6; % kg/m^3, Fluid Density
mu = 8.931e-4; % kg/m-s, Fluid Viscosity 
Dp = 0.208e-3;  % m, Particle Diametere
T = 298.15;      %  Kelvin, Temperature
const = rhof*Dp/mu; % Constant used to calculate Re, Re = const*v 
const1 = (4*g*(rhop-rhof)*Dp/3/rhof)^0.5 ; 
vt = 1 ;   % m/s, Terminal velocity guess value
e = 1e-6   ;       % Error Margin
%% Functions Used :
% Below I have defined the Functions that I have used 
Re = @(v) const*v ; % Reynolds Number 
Cd_1 = @(v) 24/(const*v) ; % Drag Coefficient for Re < 0.1 ; 
Cd_2 = @(v) (24/(const*v))*(1 + 0.14*((const*v)^0.7)); % Drag Coefficient for 0.1 =< Re =< 1000 ;
Cd_3 = @(v) 0.44; % Drag Coefficient for 1000 =< Re =< 350,000 ;
Cd_4 = @(v) 0.19 - (8e4)/(const*v); % Drag Coefficient for Re > 350,000 ;
v_t = @(Cd) const1*(Cd^-0.5) ; % Terminal Velocity in terms of Drag Coefficient
%% Solving the First Part : 
% Here I have shown how to solve the first part : 
%
% I will run a loop where first we calculate the Reynolds Number with a
% guess value of vt. Using this we obtain a new value of vt and check if
% the guess value matches. If  not, we take the obtained vt as the guess
% value in the next loop. 
%
% The program will print the guess and obtained
% guess as the loop runs so that we can observe the convergance.
i = 0; disp_vt = [vt]; disp_vtt = [0];
while true
    Re_g = Re(vt) ; 
    if Re_g < 0.1 
            Cd = Cd_1(vt) ;

    elseif (Re_g >= 0.1)||(Re_g <= 1000)
            Cd = Cd_2(vt) ;

    elseif (Re_g <= 350000)||(Re_g > 1000)
            Cd = Cd_3(vt) ; 

    elseif Re_g > 350000
            Cd = Cd_4(vt) ;
    end  
    vt_t = v_t(Cd);
    disp_vt = [disp_vt;vt]; disp_vtt = [disp_vtt;vt_t]; 
    
    if abs(vt - vt_t) < e
        break
    else
        vt = vt_t;
        i = i + 1;
    end
end
ind = [0:i+1];
T1 = table(ind',disp_vt,disp_vtt,'VariableNames',{'Iteration Number','Guess Value','Obtained Guess'});
disp(T1);
disp(['The Terminal velocity is : ',num2str(vt_t ),' m/s'])
disp(['The Reynolds Number is : ',num2str(vt*const)])
disp(['The Drag Coefficient is : ',num2str(Cd)]);
%% Solving for the Second Part : 
% Here we need to make a small change.
%
% Instead of accelration due to
% gravity, our g changes from g to what is given i.e 30*g. 
%
% Redefining the
% functions to include the changes
g_new = 30*g ; 
const2 = (4*g_new*(rhop-rhof)*Dp/3/rhof)^0.5 ;
v_t_c = @(Cd) const2*(Cd^-0.5) ; % Terminal Velocity in terms of Drag Coefficient with different g
% Again solving as done before:
% First Redefine the initial value of vt !!IMPORTANT!!
vt_c = 1;
j = 0; disp_vtc = [vt_c]; disp_vttc = [0];
while true
    Re_g = Re(vt_c) ; 
    if Re_g < 0.1 
            Cd = Cd_1(vt_c) ;

    elseif (Re_g >= 0.1)||(Re_g <= 1000)
            Cd = Cd_2(vt_c) ;

    elseif (Re_g <= 350000)||(Re_g > 1000)
            Cd = Cd_3(vt_c) ; 

    elseif Re_g > 350000
            Cd = Cd_4(vt_c) ;
    end  
    vt_t_c = v_t_c(Cd);
    disp_vtc = [disp_vtc;vt_c]; disp_vttc = [disp_vttc;vt_t_c];
    
    if abs(vt_c - vt_t_c) < e
        break
    else
        vt_c = vt_t_c;
        j = j + 1;
    end
end
ind = [0:j+1];
T2 = table(ind',disp_vtc,disp_vttc,'VariableNames',{'Iteration Number','Guess Value','Obtained Guess'});
disp(T2);
disp(['The Terminal velocity in Centrifugal Separator is : ',num2str(vt_t_c ),' m/s'])
disp(['The Reynolds Number is : ',num2str(vt_c*const)])
disp(['The Drag Coefficient is : ',num2str(Cd)])


    
