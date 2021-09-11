%% Question 2 Harikrishnan R N, 18CHE147
% Here I have shown how to code for the Stead-State Material Balance on a
% Separation Train. 
%% Defining the Constants : 
% Here I have defined the constant compositions that are used to solve the
% equations of Mass Balance : 
%The row vectors are in the order : D1,B1,D2,B2 are the respective
%compositions according to the question
xf1 = [ 15 ; 25 ; 40 ; 20 ]/100 ; 
xd1 = [ 7 ; 4 ; 54 ; 35 ]/100 ; 
xb1 = [ 18 ; 24 ;  42 ; 16 ]/100 ; 
xd2 = [ 15 ; 10 ; 54 ; 21 ]/100 ; 
xb2 = [ 24 ; 65 ; 10 ; 1 ]/100 ; 
F = 70 ;        % mol/min, Feed Flowrate 
%% Solving for the Feed Column : 
% Here I will find the molar flowrates of each species : 
%Solving the System of Linear Equations defined below :
lin_sys = [ xd1 xb1 xd2 xb2 ] ;
flow_sol = F*inv(lin_sys)*xf1 ; 
d1 = flow_sol(1) ; b1 = flow_sol(2);
d2 = flow_sol(3) ; b2 = flow_sol(4);
mol_flow_rates = [d1;b1;d2;b2];
d = d1+d2 ; b = b1+b2 ; 
%% Solving for Second Column : 
% Here I will find the compositions of each species in #2 column :
%Defining the funcion to solve the equations : 
F1 = @(x) [ x(1)*(d1+b1) - xd1(1)*d1 - xb1(1)*b1;
            x(2)*(d1+b1) - xd1(2)*d1 - xb1(2)*b1;
            x(3)*(d1+b1) - xd1(3)*d1 - xb1(3)*b1;
            x(4)*(d1+b1) - xd1(4)*d1 - xb1(4)*b1; ] ;
x0 = [ 0.1 ; 0.1 ; 0.1 ; 0.1] ; 
options = optimoptions('fsolve','Display','None');
xd = fsolve(F1,x0,options) ; 
%% Solving for Third Column : 
% Here I will find the compositions of each species in #3 column :
%Defining the funcion to solve the equations : 
F1 = @(x) [ x(1)*(d2+b2) - xd2(1)*d2 - xb2(1)*b2;
            x(2)*(d2+b2) - xd2(2)*d2 - xb2(2)*b2;
            x(3)*(d2+b2) - xd2(3)*d2 - xb2(3)*b2;
            x(4)*(d2+b2) - xd2(4)*d2 - xb2(4)*b2; ] ;
x0 = [ 0.1 ; 0.1 ; 0.1 ; 0.1] ; 
xb = fsolve(F1,x0,options) ;
%% Displaying Results : 
streams = ['D1';'B1';'D2';'B2' ] ;
tab1 = table(streams,[mol_flow_rates],'VariableNames',{'Stream','Molar Flowrate, mol/min'});
species = ["Xylene";"Styrene";"Toluene";"Benzene"];
tab2 = table(species,xd,xb,'VariableNames',{'Compound','Stream D','Stream B'});
disp(" Below is the Flowrates of each stream : ");
disp(tab1); 
disp(" Below is the Composition of each stream : ");
disp(tab2);

