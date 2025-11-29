clear all
clc


load('Clough-Initial-0.95miu.mat');
xinitial = xout;
mowei = length(xout{1}.Values.Time);
del = 1:mowei-1;
for i = 1:1:numElements(xout)
xinitial{i}.Values = delsample(xinitial{i}.Values,"Index",del);
end

clear xout SimulationMetadata logsout


% Global property
g= 9.80665;
NumInterFace = 4;
StiModi = 1;


% Unit converter
length_Convert = 0.3048; % meter/foot;
weight_Convert =  0.453592; % kg/lb;
force_Convert =  4448.2216; % N/kip;
inch_Convert = 1/39.37;% meter/inch;
stiffness_Convert = force_Convert/inch_Convert; %(N/m)/(kip/in)
young_Convert = force_Convert/(inch_Convert^2);% (N/m^2)/(kip/inch^2)

% Deadweight property

chang = 4*length_Convert;%4.6*length_Convert;
kuan = 6*length_Convert;
gao = 1/4*chang;
mass = 8000*weight_Convert;
rho_dw = mass/(chang*kuan*gao);

% Column and girder material properties 
rho_steel = 0.28*weight_Convert/(inch_Convert^3); 
E_steel = 29000*young_Convert;
p_steel = 0.26;

% Material damping ratio 
uni_dr = 0.05;

% Column geometry w5x16
col_depth = 5.01*inch_Convert;
col_ft = 0.36*inch_Convert;
col_fw = 5*inch_Convert;
col_wt = 0.24*inch_Convert;
col_area = (col_depth-2*col_ft)*col_wt + 2*(col_fw*col_ft);


% Girder geometry w6x12
beam_depth = 6.03*inch_Convert;
beam_ft = 0.28*inch_Convert;
beam_fw = 4*inch_Convert;
beam_wt = 0.23*inch_Convert;
beam_area = (beam_depth-2*beam_ft)*beam_wt + 2*(beam_fw*beam_ft);


% Length of different components
L_col = 6*length_Convert + 9*inch_Convert;
L_col2 = 5*length_Convert + 4*inch_Convert;
L_col3 = 5*length_Convert + 4*inch_Convert;

L_beam = 4*length_Convert; 
L_beam2 = 2*length_Convert; 
L_inplane = 6*length_Convert;

frame_rot = atan((L_beam*3/2)/(L_col/2));
L_brace =L_beam*3/2 /sin(frame_rot);




maS = mass*3 + col_area*(L_col+L_col2+L_col3)*rho_steel*4+beam_area*(3*L_beam*2+2*L_inplane)*rho_steel*3+4*rho_steel*beam_area*L_brace;

% Stiffness and damping for one horizontal frictional element.

muf = 1e10;
ksin=0.55;                 %Damping ratio.
fh = 50;                 %Normally from 10 to 100 hz.
wn = 2*pi*fh;
k_n = (wn^2)*(maS)/NumInterFace;
c_n=ksin*2*sqrt(k_n*NumInterFace*(maS))/NumInterFace;       %Damping for one frictional element.

% Stiffness and damping for one vertical support element.

dR = 0.34;                   % Damping ratio
kk2=40*stiffness_Convert;   %Stiffnes for one vertical element.
Ksupport_total = kk2*4;
Csupport_total = dR*2*sqrt(Ksupport_total*(maS));
cc1=Csupport_total/NumInterFace;   %Damping for one vertical element.

%Initial conditions
hc = 0;
initialRotation = 0;
% Krot = maS*g*rcg*cos(initialRotation+(pi/2-alphacg))/(alphacg-initialRotation);
% Demand = maS*g*rcg*cos(initialRotation+(pi/2-alphacg));
zc =1* (hc/2 - (maS)*g/(Ksupport_total));


%Input signal

gamma=0;
Critical = [];
Critical2 = [];
% Solver tolerance
relTol = 1e-4; %was 1e-7
absTol = 1e-4;


% For loop over pulse frequency and amplitude


t_end = 8.93;%Cc+2*pi/omegag;
cr = 1;

% set_param('CloughExperiment','LoadInitialState','on','InitialState','xinitial')
% tic;
% sim('CloughExperiment.slx'); %Run the simulink model.
% toc;
% 

