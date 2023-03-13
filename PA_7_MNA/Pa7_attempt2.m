%% Define known parameters

% Resistances
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
Ro = 1000;

% Conductances

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
Go = 1/Ro;

% source conversion factor Alpha
alpha = 100;

% conductance
C1 = 0.25;

% Inductance
L1 = 0.2;

%% Making matrices

Vin = 1; % will sweep from -10 V to 10 V

% 6 unknowns:
% x = [Is, V2, IL, V3, V5, I4, V1], note x(4) = V3, and x(5) = V5

% G matrix: (it takes me 4 times to actually get it right)

% G = [R1,1,0,0,0,0 ;...
%      0,-((R1/R2) - 1),-R1,0,0,0 ;...
%      0,(1 - R1/R2),0,(-R1/R3),0,0 ;...
%      0,(1 - R1/R2),0,0,(-R1/alpha),(R1*R4/alpha) ;...
%      0,(1 - R1/R2),0,0,(-(R1/alpha)*(1 - R4/Ro)),0 ;...
%      R1,0,0,1,0,0 ] ;
 
% G matrix attempt 2;
% G = [(R1),(1),0,0,0,0,0 ;...
%      0,(1 - R1/R2),(-R1),0,0,0,0 ;...
%      0,0,0,0,0,0,1 ;...
%      0,0,(-1),(1/R3),0,0,0 ;...
%      0,0,0,(alpha/(R3*R4)),(-1/R4),(-1),0 ;...
%      0,0,0,(alpha/(R3*R4)),(1/Ro - 1/R4),0,0 ;...
%      0,(1),0,(-1),0,0,0 ] ;
 
 % G matrix attempt 3;
% G = [0,0,0,0,0,0,1 ;...
%      1,G1,0,0,0,0,-G1 ;...
%      0,(G1-G2),-1,0,0,0,-G1 ;...
%      0,0,-1,G3,0,0,0 ;...
%      0,0,0,(G4*alpha/R3),(-G4),(-1),0 ;...
%      0,0,0,(G4*alpha/R3),(Go-G4),0,0 ;...
%      0,1,0,-1,0,0,0 ] ;

% 6 unknowns:
% V = [Is, V2, IL, V3, V5, I4, V1]

 % G matrix attempt 4; (wow this one actually works)
G = [0,0,0,0,0,0,1 ;...
     1,G1,0,0,0,0,-G1 ;...
     0,(-G1-G2),-(-1),0,0,0,G1 ;...
     0,0,-(-1),G3,0,0,0 ;...
     0,0,0,(G4*alpha/R3),(-G4),(-1),0 ;...
     0,0,0,(G4*alpha/R3),(Go-G4),0,0 ;...
     0,1,0,-1,0,0,0 ] ;
 
 % F vector:
 F = [Vin; 0; 0; 0; 0; 0; 0];
 
 % C matrix;
 C = [0,0,0,0,0,0,0 ;...
      0,-C1,0,0,0,0,C1 ;...
      0,-C1,0,0,0,0,C1 ;...
      0,0,0,0,0,0,0 ;...
      0,0,0,0,0,0,0 ;...
      0,0,0,0,0,0,0 ;...
      0,0,-L1,0,0,0,0 ] ;
 
%% Time Domain (TD)
i = 1;
for Vin = -10:1:10
    
    F = [Vin; 0; 0; 0; 0; 0; 0];
    
    x = inv(G)*F;
    
    gain_t_v3(i) = x(4); % get gain at V3
    gain_t_v5(i) = x(5); % get gain at V5

    i = i + 1;
    
end

%% Frequency Domain (FD)
% redefine Vin and F vector:
Vin = 1;
F = [Vin; 0; 0; 0; 0; 0; 0];

j = 1; 
for w = 0:100
    Vin = 1;
        
    F = [Vin; 0; 0; 0; 0; 0; 0];
    A = G + 1i*w*C;

    x = inv(A)*F;
    
    gain_f_V3(j) = abs(x(4)/Vin); % get gain at V3
    gain_f_V5(j) = abs(x(5)/Vin); % get gain at V5

    j = j + 1;
    
end

%% Monte Carlo stuff
random_c = 0.25 + 0.05*randn(1000); % base value of 1000,
                                    % 0.05 standard deviation

Vin = 1; % set input voltage as 1V for the AC operating point
w = 3*pi;
F = [Vin; 0; 0; 0; 0; 0; 0]; % redefine F vector with Vin = 1 V                                    

for i = 1:size(random_c) % go through every entry of random_c
    cRanVal = random_c(i); 

    % make new C matrix, C_ran
    C_ran = [0,0,0,0,0,0,0 ;...
          0,-cRanVal,0,0,0,0,cRanVal ;...
          0,-cRanVal,0,0,0,0,cRanVal ;...
          0,0,0,0,0,0,0 ;...
          0,0,0,0,0,0,0 ;...
          0,0,0,0,0,0,0 ;...
          0,0,-L1,0,0,0,0 ] ;

    A = G + 1i*w*C_ran; % solve matrix for
    x = inv(A)*F;
    
    gain_C(i) = abs(x(5)/Vin); % get gain

end

%% Plot stuff

Vin = -10:1:10;
w = 0:100;

figure('name', 'Pa 7: MNA Building')
subplot(2,2,1) % Plot the gain with respect to input voltage Vin (DC)
plot(Vin,gain_t_v3);hold on
plot(Vin,gain_t_v5);hold off
xlabel('V_{in}')
ylabel('V')
legend('V_{3}','V_{out}', 'location', 'southeast')

subplot(2,2,2) % Plot the gain with respect to frequency (AC)
plot(w,gain_f_V3);hold on
plot(w,gain_f_V5);hold off
xlabel('w')
ylabel('Gain')
legend('V_{3}','V_{out}')

subplot(2,2,3) % Plot the random perturbations on C
histogram(random_c);
xlabel('C')
ylabel('Number')

subplot(2,2,4) % Plot the gain with the random perturbations on C
histogram(gain_C);
xlabel('gain')
ylabel('Number')
