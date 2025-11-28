function [Kr, Mr, W, Tr, Np] = craigBamptonImpl(self,Ko,Mo,AttachmentDoF,Tp,numModesORFreqRange)
% craigBamptonImpl - Internal function to perform Craig-Bampton reduction.

% Copyright 2022-2023 The MathWorks, Inc.

% ZY modified version so the numModesOrFreqRange accepts a list of
% frequency range, whereas the original implementation only accpets one
% specific range. 

% Find the changes at line 45.  

dim = size(self.Mesh.Nodes,1);
numNodes = size(self.Mesh.Nodes,2);
CondensedDoF = true(dim*numNodes, 1);
CondensedDoF(AttachmentDoF) = false;

kcc = Ko(CondensedDoF,CondensedDoF); %Internal physical dofs K matrix Ko has no MPC yet 
mcc = Mo(CondensedDoF,CondensedDoF); %Internal physical dofs M matrix Mo has no MPC yet 
dK = [];
p = [];
if isscalar(numModesORFreqRange)
    % Number of modes requested
    [dK, p] = setupDecomposition(kcc);
    if ~isempty(p)
        [Phi,l_fim] = eigs(dK,size(kcc,1),mcc(p,p),...
          numModesORFreqRange,'smallestabs',...
          'IsFunctionSymmetric', true, 'IsSymmetricDefinite', true);
        Phi(p,:) = Phi;
    else
        [Phi,l_fim] = eigs(dK,size(kcc,1),mcc,...
          numModesORFreqRange,'smallestabs',...
          'IsFunctionSymmetric', true, 'IsSymmetricDefinite', true);
    end
    l_fim = diag(l_fim);
else
    opts.FrequencyRange = numModesORFreqRange;
    if ~isempty(self.SolverOptions.MaxShift)
        opts.MaxShift = self.SolverOptions.MaxShift;
    end

    if ~isempty(self.SolverOptions.BlockSize)
        opts.BlockSize = self.SolverOptions.BlockSize;
    end
    Opts = opts;
    % ZY Changed the following so phi includes multiple ranges of mode
    % shapes
    for i = 1:1:size(opts.FrequencyRange,1)
        Opts.FrequencyRange = opts.FrequencyRange(i,:);
        [Phi,l_fim] = matlab.internal.math.lanczos(kcc,mcc,Opts);
        phi{i} = Phi;
        L_fim{i} = l_fim;
    end
        Phi = [];
        l_fim = [];
    for i = 1:1:size(phi,2)
        Phi = [Phi phi{i}];
        l_fim = [l_fim;L_fim{i}];
    end
end


if isempty(Phi)
    Phi = zeros(size(kcc,1),0);
end

if ~isempty(Tp)
    %  Form blocks for Kr
    kcaT = Ko(CondensedDoF,AttachmentDoF)*Tp;
    if isempty(dK)
        Psi = -(kcc\full(kcaT));
    else
        if isempty(p)
            Psi = -dK(full(kcaT));
        else
            Psi(p,:) = -dK(full(kcaT(p,:)));
        end
    end



    Kr11 = Tp'*(Ko(AttachmentDoF,AttachmentDoF)*Tp)+kcaT'*Psi;
    %  Form blocks for Mr
    mcaT = Mo(CondensedDoF,AttachmentDoF)*Tp;
    Mac = mcaT'*Psi;
    Mr11 = Tp'*(Mo(AttachmentDoF,AttachmentDoF)*Tp)+Mac+Mac'+Psi'*(mcc*Psi);
    Mr12 = mcaT'*Phi+Psi'*(mcc*Phi);
    %  Form Tr
    Tr = matlab.internal.math.blkdiag(Tp, speye(size(kcc)));
else
    %  Form blocks for Kr
    kca = Ko(CondensedDoF,AttachmentDoF);
    if isempty(dK)
        Psi = -(kcc\(full(kca)));
    else
        if isempty(p)
            Psi = -dK(full(kca));
        else
            Psi(p,:) = -dK(full(kca(p,:)));
        end
    end
    Kr11 = Ko(AttachmentDoF,AttachmentDoF)+kca'*Psi;
    %  Form blocks for Mr
    mca = Mo(CondensedDoF,AttachmentDoF);
    Mac = mca'*Psi;
    Mr11 = Mo(AttachmentDoF,AttachmentDoF)+Mac+Mac'+Psi'*(mcc*Psi);
    Mr12 = mca'*Phi+Psi'*(mcc*Phi);
    %  Form Tr
    Tr = [];
end

%  Form projections.  Explicitly symmetrize diagonal blocks.
Kr22 = Phi'*(kcc*Phi);
Kr = blkdiag((Kr11+Kr11')/2, (Kr22+Kr22')/2);

Mr22 = Phi'*(mcc*Phi);
Mr = [(Mr11+Mr11')/2  Mr12; Mr12'  (Mr22+Mr22')/2];

%  Form W
Np = numel(l_fim);
W  = [eye(size(Psi,2), size(Psi,2)+Np); Psi Phi];

end

function [fh, p] = setupDecomposition(kcc)
% Use nested dissection instead of the default and try Cholesky explicitly
p = dissect(kcc);
dK = matlab.internal.decomposition.SparseCholesky(kcc(p,p), false);
if success(dK)
    fh = @(x)solve(dK, x);
else
    % Cholesky failed, but we expect that kcc is symmetric still.
    % Revert to LDL
    p = [];
    dK = decomposition(kcc, 'ldl', 'CheckCondition', false, ...
      'AllowIterativeRefinement', false);
    fh = @(x)dK\x;
end
end
