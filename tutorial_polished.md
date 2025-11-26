# Tutorial Documentation

## Overview

## Pre-requisites

2. It is recommended that the readers read this well-written [matlab example](https://uk.mathworks.com/help/sm/ug/model-excavator-dipper-arm.html) to be familiar with the workflow of using reduced order flexible solid in Simulink.

## Why This Separate Documentation

First, the provided codes additionally add two useful functionalities to the original matlab implementation: 1. Allows to keep any desired fixed-interface modes while the original implementation only keeps fixed-interface modes within one specific frequency range. 2. Allows to reconstruct the full-order deformation, strain, stress and global displacements (add credits to houmain, inform Manolis on this, GPT do not change this line) from the Simulink simulation results while the original implementation only returns global displacements of the interface nodes. The first example will cover this.

Secondly, the Matlab example is not a rocking problem and hence it does not provide guides for using the OxContact elements (an in-house toolbox for the contact modelling in rocking problems) and the reduced order flexible solid together. Both examples will cover this.

Thirdly, the Matlab example is directly applicable only for structural components where the interface is stiff and can be safely assumed to be rigid. For non-rigid type interface, for example rocking components without interface reinforcement, the implementation requires extra care. The first example will exclusively cover this.

Fourthly, the Matlab example involves only few physical components whereas civil engineering structures such as frame structures, normally include many inter-connected components. Building such a structure in Simulink, which is based on floating frame of reference formulation, requires much more efforts than and is not as straightforward as traditional civil engineering software. The second example will exclusively cover this.

Lastly, the floating frame of reference formulation and the craig-bampton reduction method, served as the workhorse in the numerical framework proposed in the paper, remain a relative stranger for civil engineering engineers and researchers, despite its great potentials in addressing rocking problems. The examples covering the implementation of different types of rocking structure should improve their familiarity.

## General Workflow

Craig-bampton is carried out in matlab, then the information (importantly the reduced order stiffness, mass and damping matrices) is put to the reduced order flexible solid in Simulink, which uses floating frame of reference formulation, for dynamic simulation. Finally, reconstruct full order solution.

The first example covers the detailed workflow and the steps that are not covered by the matlab example are underlined.

The first example models a solid marble column with four rigid corner feet standing on a high grip, stiff support medium; the height of the corner feet is negligible. Figure shows the geometry of the column (insert the figure, GPT do not change this line):

Steps for Craig-Bampton reduction in matlab

Define the parameters for the structure and other common parameters

```matlab
clear all

clc

NumInterFace = 4;

StiModi = 1;

g = 9.80665;

Width = 1;

Depth = 1*1.05;

E = 5.42538e09;

rho = 2256.86;

nu = 0.3;

L = 10;

A = Width*Depth;

ms = rho*L*A;

CrossSize = [Width/2 Depth/2; -Width/2 Depth/2; -Width/2 -Depth/2; Width/2 -Depth/2];

% Define the length of tetrahedral model and the beam model

L_te =L;

Subdi = 10; % Subdivide the beam element into 10, although there are only 8 beam elements. The base is

% replaced with the tetraheral element which is 20% height.

Le = L/Subdi;

% Global properties

distanceToEdge = 0.2;

I = ms*(Width^2+L^2)/12;

hcg = L/2;

rcg = sqrt((hcg^2) + (Width/2*(1-distanceToEdge))^2);

p = sqrt(ms*g*rcg/I);

maS = ms;               %Total mass of the system.

alphacg = atan((Width/2*(1-distanceToEdge))/hcg);

% Stiffness and damping for one horizontal frictional element.

muf = 0.8;

fh = 60; %21                %Normally from 10 to 100 hz.

wN=2*pi*fh;

ksin=0.5;                 %Damping ratio.

k_n=(maS)*wN^2/NumInterFace/StiModi;           %Stiffnes for one frictional element.

c_n=ksin*2*sqrt(k_n*NumInterFace*(maS))/NumInterFace;      %Damping for one frictional element.

% Stiffness and damping for one vertical support element.

fv = 60;                 %Normally from 10 to 100 hz.

wn = 2*pi*fv;

dR = 0.5;

Ksupport_total = (wn^2)*(maS)/StiModi;%Sqrt(k/m) gives rad/s

Csupport_total = dR*2*sqrt(Ksupport_total*(maS));

kk2=Ksupport_total/NumInterFace;   %Stiffnes for one vertical element.

cc1=Csupport_total/NumInterFace;   %Damping for one vertical element.

%Initial conditions

hc = 0;

initialRotation = 0;

Krot = maS*g*rcg*cos(initialRotation+(pi/2-alphacg))/(alphacg-initialRotation);

Demand = maS*g*rcg*cos(initialRotation+(pi/2-alphacg));

zc = 1*(hc/2 - (maS)*g/(Ksupport_total));

```matlab
Generate or import geometry of the component for which the craig-bampton reduction will be needed. Details can found in this [example]( https://uk.mathworks.com/help/sm/ug/model-excavator-dipper-arm.html).

```matlab
gm = multicuboid(Width,Depth,L_te);

LF_coords = [-Width/2+distanceToEdge*Width/2;-Depth/2+distanceToEdge*Depth/2;0];

RF_coords = [Width/2-distanceToEdge*Width/2;-Depth/2+distanceToEdge*Depth/2;0];

LB_coords = [-Width/2+distanceToEdge*Width/2;Depth/2-distanceToEdge*Depth/2;0];

RB_coords = [Width/2-distanceToEdge*Width/2;Depth/2-distanceToEdge*Depth/2;0];

RF_dis = RF_coords - LF_coords;

LB_dis = LB_coords - LF_coords;

RB_dis = RB_coords - LF_coords;

origins = [LF_coords';RF_coords';LB_coords';RB_coords'];

numFrames = size(origins,1);

addVertex(gm,"Coordinates",origins);

PoI = [-Width/2 0 0.1*L_te;-Width/2 0 0.2*L_te;-Width/2 0 0.3*L_te;-Width/2 0 0.4*L_te;-Width/2 0 0.5*L_te

-Width/2 0 0.6*L_te;-Width/2 0 0.7*L_te;-Width/2 0 0.8*L_te;-Width/2 0 0.9*L_te;-Width/2 0 L_te

Width/2 0 0.1*L_te;Width/2 0 0.2*L_te;Width/2 0 0.3*L_te;Width/2 0 0.4*L_te;Width/2 0 0.5*L_te;Width/2 0 0.6*L_te

Width/2 0 0.7*L_te;Width/2 0 0.8*L_te;Width/2 0 0.9*L_te;Width/2 0 L_te];

for i = 1:1:size(PoI,1)

addVertex(gm,"Coordinates",PoI(i,:));

end

model = createpde('structural','modal-solid');

model.Geometry = gm;

pdegplot(model,'EdgeLabels','on','FaceAlpha',0.5);

```matlab
Mesh the component. This is where the position of interface node is decide. Mesh first time to know position of the nearby vertex, then add the vertex and mesh again. Add vertex type multipoint constrains. The multipoint constrain is needed because tetrahedral elements node only has three dofs, however an interface should have six dofs in Simulink, so the multipoint constrain essentially forms a small rigid plane that has six dofs.

```matlab
%Specify structural properties

structuralProperties(model,"YoungsModulus",E, ...

"PoissonsRatio",nu, ...

"MassDensity",rho);

% Create finite-element mesh

hmax = 0.2;

msh=generateMesh(model,"Hface",{1,0.06,2,0.2},Hmax=hmax,GeometricOrder="linear",Hgrad=1.2);

```matlab
<ins>Locate the node positions that are within the area of the feet</ins>

```matlab
nodes = msh.Nodes;

elements = msh.Elements;

LF_Node_Index = 9;

RF_Node_Index = 10;

LB_Node_Index = 11;

RB_Node_Index = 12;

ii = 1;

for i = 1:1:size(msh.Elements,2)

if sum(find (msh.Elements(:,i) == LF_Node_Index)) ~= 0

ElemInLF(ii) = i;

ii = ii+1;

else

end

end

ii = 1;

for i = 1:1:size(msh.Elements,2)

if sum(find (msh.Elements(:,i) == RF_Node_Index)) ~= 0

ElemInRF(ii) = i;

ii = ii+1;

else

end

end

ii = 1;

for i = 1:1:size(msh.Elements,2)

if sum(find (msh.Elements(:,i) == LB_Node_Index)) ~= 0

ElemInLB(ii) = i;

ii = ii+1;

else

end

end

ii = 1;

for i = 1:1:size(msh.Elements,2)

if sum(find (msh.Elements(:,i) == RB_Node_Index)) ~= 0

ElemInRB(ii) = i;

ii = ii+1;

else

end

end

Ele_link_LF_index = elements(:,ElemInLF); %Node index of elements that contains the reference node

Ele_link_LF_index2 = Ele_link_LF_index(1:4,:); % Only include the element vertex. The others are midpoint node

adc = unique(Ele_link_LF_index2);%Get rid of the repeated element vertex.

adc(find(adc == LF_Node_Index)) = [];%Get rid of the reference node

LF_surroundingNodes_Index = adc(find(nodes(3,adc)==0));

LF_surroundingNodes_coords = nodes(:,LF_surroundingNodes_Index);%Find the node coordinates that is one element away from the reference node. This will be used to add vertex to the geometry in the next step.

Ele_link_RF_index = elements(:,ElemInRF); %Node index of elements that contains the reference node

Ele_link_RF_index2 = Ele_link_RF_index(1:4,:); % Only include the element vertex. The others are midpoint node

adc = unique(Ele_link_RF_index2);%Get rid of the repeated element vertex.

adc(find(adc == RF_Node_Index)) = [];%Get rid of the reference node

RF_surroundingNodes_Index = adc(find(nodes(3,adc)==0));

RF_surroundingNodes_coords = nodes(:,RF_surroundingNodes_Index);%Find the node coordinates that is one element away from the reference node. This will be used to add vertex to the geometry in the next step.

Ele_link_LB_index = elements(:,ElemInLB); %Node index of elements that contains the reference node

Ele_link_LB_index2 = Ele_link_LB_index(1:4,:); % Only include the element vertex. The others are midpoint node

adc = unique(Ele_link_LB_index2);%Get rid of the repeated element vertex.

adc(find(adc == LB_Node_Index)) = [];%Get rid of the reference node

LB_surroundingNodes_Index = adc(find(nodes(3,adc)==0));

LB_surroundingNodes_coords = nodes(:,LB_surroundingNodes_Index);%Find the node coordinates that is one element away from the reference node. This will be used to add vertex to the geometry in the next step.

Ele_link_RB_index = elements(:,ElemInRB); %Node index of elements that contains the reference node

Ele_link_RB_index2 = Ele_link_RB_index(1:4,:); % Only include the element vertex. The others are midpoint node

adc = unique(Ele_link_RB_index2);%Get rid of the repeated element vertex.

adc(find(adc == RB_Node_Index)) = [];%Get rid of the reference node

RB_surroundingNodes_Index = adc(find(nodes(3,adc)==0));

RB_surroundingNodes_coords = nodes(:,RB_surroundingNodes_Index);%Find the node coordinates that is one element away from the reference node. This will be used to add vertex to the geometry in the next step.

```matlab
<ins>Add the node positions as vertex to the original geometry and re-mesh the new geometry.</ins>

```matlab
model2 = createpde("structural","modal-solid");

gm2 = multicuboid(Width,Depth,L_te);

addVertex(gm2,"Coordinates",origins);

addVertex(gm2,"Coordinates",[LF_surroundingNodes_coords';RF_surroundingNodes_coords';LB_surroundingNodes_coords';RB_surroundingNodes_coords']);

for i = 1:1:size(PoI,1)

addVertex(gm2,"Coordinates",PoI(i,:));

end

model2.Geometry = gm2;

generateMesh(model2,"Hface",{1,0.06,2,0.2},Hmax=hmax,GeometricOrder="linear",Hgrad=1.2); %need to be kept the same as the msh method in model in the last section

structuralProperties(model2,"YoungsModulus",E, ...

"PoissonsRatio",nu, ...

"MassDensity",rho);

```matlab
<ins>Apply vertex type multipoint-constrain to define the interface, the degrees of freedom of which will be retained after Craig-Bampton reduction.</ins>

```matlab
LF_surroundingVertex_Index = [1:1:size(LF_surroundingNodes_Index,1)]+12;

RF_surroundingVertex_Index = [size(LF_surroundingNodes_Index,1)+1:1:size(RF_surroundingNodes_Index,1)+size(LF_surroundingNodes_Index,1)]+12;

LB_surroundingVertex_Index = [size(RF_surroundingNodes_Index,1)+size(LF_surroundingNodes_Index,1)+1:1:size(RF_surroundingNodes_Index,1)+size(LF_surroundingNodes_Index,1)+size(LB_surroundingNodes_Index,1)]+12;

RB_surroundingVertex_Index = [size(RF_surroundingNodes_Index,1)+size(LF_surroundingNodes_Index,1)+size(LB_surroundingNodes_Index,1)+1:1:size(RF_surroundingNodes_Index,1)+size(LF_surroundingNodes_Index,1)+size(LB_surroundingNodes_Index,1)+size(RB_surroundingNodes_Index,1)]+12;

LF_surroundingVertex_Index = [9 LF_surroundingVertex_Index];

RF_surroundingVertex_Index = [10 RF_surroundingVertex_Index];

LB_surroundingVertex_Index = [11 LB_surroundingVertex_Index];

RB_surroundingVertex_Index = [12 RB_surroundingVertex_Index];

structuralBC(model2, ...

'Vertex',LF_surroundingVertex_Index, ...

'Constraint','multipoint', ...

'Reference',origins(1,:));

structuralBC(model2, ...

'Vertex',RF_surroundingVertex_Index, ...

'Constraint','multipoint', ...

'Reference',origins(2,:));

structuralBC(model2, ...

'Vertex',LB_surroundingVertex_Index, ...

'Constraint','multipoint', ...

'Reference',origins(3,:));

structuralBC(model2, ...

'Vertex',RB_surroundingVertex_Index, ...

'Constraint','multipoint', ...

'Reference',origins(4,:));

figure(1)

pdegplot(model2,'Vertexlabels','on','FaceAlpha',0.5);

Q = meshQuality(model2.Mesh,"aspect-ratio");

elemIDs = find(Q < 0.5);

figure(2)

pdemesh(model2.Mesh,"FaceAlpha",0.5)

hold on

pdemesh(model2.Mesh.Nodes,model2.Mesh.Elements(:,elemIDs), ...

"FaceColor","blue","EdgeColor","blue")

mv03_percent = volume(model2.Mesh,elemIDs)/volume(model2.Mesh)*100

figure(3)

hist(Q);

xlabel('Element Shape Quality', 'fontweight','b');

ylabel('Number of Elements', 'fontweight','b');

model2.SolverOptions.MaxShift = 500;

```matlab
Check the fixed-interface modes of the component and take notes of the frequency of the modes you want to keep. <ins>Optionally, you can run the modal analysis as a live script and use this nice live task  [‘visualize PDE results’]( ) to visualize the mode shapes</ins>.

```matlab
%% Modal analysis

modalresults = solve(model2,FrequencyRange=[0,3e3]); %in rad/s.

Freq = modalresults.NaturalFrequencies;

```matlab
(insert the modeshape figure here, GPT do not change this line.)

Apply the Craig-Bampton reduction to obtain the reduced stiffness and mass matrices. <ins>Put the list of frequency range corresponding to the modes you want to keep.</ins>

```matlab
%% Reduced order

% Pseudo-full order of nonlinear stress model

FreqRange = [0 Freq(20)+2;

Freq(24)+2 Freq(26)-2;

Freq(29)+2 Freq(31)-2;

Freq(35)+2 Freq(37)-2;

Freq(40)+2 Freq(42)-2];

R = reduce(model2,"FrequencyRange",FreqRange); %The frequency rang is in rad/s

Reduced.K = (R.K+R.K')/2;                    % Reduced stiffness matrix

Reduced.M = (R.M+R.M')/2;                    % Reduced mass matrix

Reduced.P = R.ReferenceLocations';

frmPerm = zeros(numFrames,1);  %was zeros(numFrames,1)  % Frame permutation vector

dofPerm = 1:size(Reduced.K,1);       % DOF permutation vector

%assert(size(Reduced.P,1) == numFrames);

for i = 1:numFrames %Was 1:numFrames

for j = 1:numFrames %Was 1:numFrames

if isequal(Reduced.P(j,:),origins(i,:))

frmPerm(i) = j;

dofPerm(6*(i-1)+(1:6)) = 6*(j-1)+(1:6);

continue;

end

end

end

Reduced.P = Reduced.P(frmPerm,:);

Reduced.K = Reduced.K(dofPerm,:);

Reduced.K = Reduced.K(:,dofPerm);

Reduced.M = Reduced.M(dofPerm,:);

Reduced.M = Reduced.M(:,dofPerm);

```matlab
Manually calculate the reduced order damping matrix. The following code allows to assign different modal damping ratios to different modes. Optionally, you do not need to calculate the damping matrix, the Simulink reduced order flexible solid can calculate the damping matrix for uniform modal damping model and Rayleigh damping model.

```matlab
%% Manual calculation of damping matrix

MoI_Length = ones(1,8);

MoI_Dr = 0.8;

OtherModeDr = 0.8;

MaterialDRMatrix = [0 0 0 0 0 0 MoI_Dr*MoI_Length OtherModeDr*ones(1,length(Reduced.K)-6-length(MoI_Length))];

MaterialDRMatrix = diag(MaterialDRMatrix);

[V,D] = eig(Reduced.K,Reduced.M);

[d,sortIdxs] = sort(diag(D));

V = V(:,sortIdxs);

% Due to small numerical errors, the six eigenvalues associated with the

% rigid-body modes may not be exactly zero. To avoid numerical issues,

% check that the first six eigenvalues are close enough to zero. Then

% replace them with exact 0 values.

assert(all(abs(d(1:6))/abs(d(7)) < 1e-9),'Error due to "zero" eigenvalues.');

d(1:6) = 0;

% Vectors of generalized masses and natural frequencies

MV = Reduced.M*V;

generalizedMasses = diag(V'*MV);

naturalFrequencies = sqrt(d);

Reduced.C = MV * diag(2*MaterialDRMatrix*naturalFrequencies./generalizedMasses) * MV';

Reduced.C = (Reduced.C+Reduced.C')/2;

Reduced.C = Reduced.C(dofPerm,:);

Reduced.C = Reduced.C(:,dofPerm);

```matlab