filename2 = ['Column-Results-26Modes-' num2str(freqRatio) '-' num2str(AmpRatio) '.mat']; 
load(filename2);                                    % Load simulation results file.

filename3 = 'Column-Reconstruction-26Modes.mat';
load(filename3);                                    % Load the reconstruction file.

Num_Interface = 4;                                  % Number of interfaces.  
t = xout{1}.Values.Time;                            % Total time steps 
PoI = [0.5 -0.366 1.008; 0.5 -0.525 0.3];           % Point of interest, as defined before. If different from the PoI in the reconstruction file, the actual positions inquired will be the node positions closest to the PoI. 
                                                                            
                                                  
% Find the node positions closest to the PoI. 

NodesCoord = model2.Mesh.Nodes;
for i = 1:1:size(PoI,1)
    NodePosition = PoI(i,:);
    XX = NodesCoord(1,:)-NodePosition(1);
    YY = NodesCoord(2,:)- NodePosition(2);
    ZZ = NodesCoord(3,:)- NodePosition(3);
    NodeIndex(i) = intersect(intersect(find(abs(XX)<0.02),find(abs(YY)<0.02)),find(abs(ZZ)<0.02));
end
NodePosition = NodesCoord(:,NodeIndex)';



% Pre-define the size of the variables. 

time = zeros(1,length(t));
StressZZ = zeros(size(PoI,1),length(time)); 
StressXX = zeros(size(PoI,1),length(time)); 
StressYY = zeros(size(PoI,1),length(time)); 
StressYZ = zeros(size(PoI,1),length(time));
StressXZ = zeros(size(PoI,1),length(time));
StressXY = zeros(size(PoI,1),length(time));
DisX = zeros(1,length(time));
DisY = zeros(1,length(time));
DisZ = zeros(1,length(time));


% Break the total time steps into number of segNum pieces, and deal with them one at a time so the
% post-processing does not consume all the ram.

segNum = 500;                                       
segment = cell(segNum,1);
startId = 1;
for i = 2:1:segNum-1
    segment{i,1} = round((length(t)-startId)*(i-1)/segNum)+startId+1:round((length(t)-startId)*i/segNum)+startId;
end
segment{1,1} = startId:round((length(t)-startId)/segNum)+startId;
segment{end,1} = round((length(t)-startId)*(segNum-1)/segNum)+startId+1:length(t);





% Reconstruct the full order solutions

for abcd = 1:1: size(segment,1)                                 % Post-process one piece of time steps during each loop.
    Interval = segment{abcd,1};
    BodyFrameSolu = zeros(size(Interval,2),6);                  % Relative displacements, velocities and accelerations of the FFR origin to itself, must be zeros for all the six DoFs.  

    interfaceDofSolu = [BodyFrameSolu xout{19}.Values.Data(Interval,1:(Num_Interface-1)*6)];    % Inferface DoFs solutions. xout{19} is the state solution of the reduced order flexible solid.
    internDofSolu =  xout{19}.Values.Data(Interval,(Num_Interface-1)*6+1:end);                  % Internal DoFs solutions.
    simSolution = [interfaceDofSolu internDofSolu];                                             

    interfaceDofSoluDot =  [BodyFrameSolu xout{20}.Values.Data(Interval,1:(Num_Interface-1)*6)];% The first derivative of the inferface DoFs solutions. xout{20} is the first derivative of the state solution of the reduced order flexible solid.
    internDofSoluDot =  xout{20}.Values.Data(Interval,(Num_Interface-1)*6+1:end);               % The first derivative of the internal DoFs solutions.                  
    simSolutionDot = [interfaceDofSoluDot internDofSoluDot];

    u = simSolution';
    ut = simSolutionDot';
    utt = zeros(size(u,1),size(u,2));                           % The second derivative of the DoFs. Not needed for reconstructing displacement and velocity solutions, but syntax requires it so all zeros with the compatible size.
    tlist = t(Interval);
    RTrom = reconstructSolution(R,u,ut,utt,tlist);              % Run the reconstruction function.


    % strain = evaluateStrain(RTrom);
    stress = evaluateStress(RTrom);                             % Evaluate the stress.
    time(1,Interval-startId+1) = tlist;
    for aq = 1:1:length(NodeIndex)
        StressZZ(aq,Interval-startId+1) = stress.zz(NodeIndex(aq),:);
        StressXX(aq,Interval-startId+1) = stress.xx(NodeIndex(aq),:);
        StressYY(aq,Interval-startId+1) = stress.yy(NodeIndex(aq),:);
        StressXY(aq,Interval-startId+1) = stress.xy(NodeIndex(aq),:);
        StressXZ(aq,Interval-startId+1) = stress.xz(NodeIndex(aq),:);
        StressYZ(aq,Interval-startId+1) = stress.yz(NodeIndex(aq),:);
    end
end

% Save the reconstruction results

filename4 = ['Column-Stress-26Modes-' num2str(freqRatio) '-' num2str(AmpRatio) '.mat']; 
save(filename4,'time','StressZZ','NodePosition','RTrom');

