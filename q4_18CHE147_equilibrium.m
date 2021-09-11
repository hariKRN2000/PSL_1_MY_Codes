%% Question 4 Harikrishnan R N, 18CHE147
% Here I have shown the code to solve question 4 of assignment: 
%
% We need to solve simultenous non-linear equations here.
%% Defining Constants : 
% Here we will define the constant terms of the code: 
Kc1 = 1.06 ; Kc2 = 2.63 ; Kc3 = 5 ; 
CA0 = 1.5 ; CB0 = 1.5 ; 
%% Defining the Reactions and Extent of Reactions : 

%
% For future references: 
%
% Index corresponding to species are : A-1,B-2,C-3,D-4,X-5,Y-6,Z-7
% Uncomment according to use
n = 7 ;% input('Enter number of Species : ') ; % Enter Number of Species 
m = 3 ;% input('Enter Number of Reactions :') ; % Enter Number of Reactions
% The reaction is given as : 
% 
% A + B --> C + D 
% B + C --> X + Y 
% A + X --> Z
%
ICM = zeros(n,1) ; % Initial Concentration Matrix
ICM(1) = CA0 ; ICM(2) = CB0 ; 
CM = zeros(n,1) ;  % Concentration Matrix
%% Defining the function to be solved : 
conc = @(C) [ (C(3)*C(4))-Kc1*(C(1)*C(2)) ;
              (C(5)*C(6))-Kc2*(C(2)*C(3)) ;
              (C(7))-Kc3*(C(3)*C(2)) ;
               C(1) - CA0 + C(4) + C(7) ; 
               C(2) - CB0 + C(4) + C(6) ; 
               C(3) - C(4) + C(6) ; 
               C(6) - C(5) - C(7) ; ];
%% Solving the function with different guess values 
% Here we will see how can we solve the above defined function for different
% guess values as given in the question : 
%The first guess: 
conc0_1 = 0.001*(ones(size(CM))) ;
conc0_1([4,5,7]) = [0 0 0];
conc_1 = fsolve(conc,conc0_1);
% Second Guess : 
conc0_2 = 0.001*(ones(size(CM))) ; 
conc0_2([4,5,7]) = [1 1 1];
conc_2 = fsolve(conc,conc0_2);
% Third Guess : 
conc0_3 = 0.001*(ones(size(CM))) ; 
conc0_3([4,5,7]) = [10 10 10];
conc_3 = fsolve(conc,conc0_3);
%% Displaying Results 
species = ['A';'B';'C';'D';'X';'Y';'Z'];
ind = 1:n;
T = table(ind',species,ICM,conc_1,conc_2,conc_3,'VariableNames',{'Sr.No','Species','Initial Concentration','First Guess','Second Guess','Third Guess'});
disp(T); disp(['Sum of Conc of species = ',num2str(sum(conc_3))])




