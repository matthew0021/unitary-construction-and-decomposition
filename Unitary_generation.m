function out = Unitary_generation(varargin)
%
% -------------------------------------------------------------------------
% Brief description: 
% 
% To compute the Unitary operator based on Phys. Rev. Lett. 120, 240502
% (2018) using inferred quantum memory states based on Physical Review A 
% 101 (3), 032327 (2020).
%
% This code was written as part of the Quantum Inference Project to build 
% a suitable unitary operator for the inferred quantum memory states.
% The code was written by Matthew Ho, a PhD candidate at the School of
% Physical and Mathematical Sciences at the Nanyang Technological
% University (NTU), Singapore & Complexity Institute, NTU, Singapore.
% 
% The code combines building quantum memory states and using the quantum
% memory states to build the unitary operator. The quantum memory states
% are merged if they are deemed to be 'equivalent', i.e. having same
% conditional probabilities of leading to the same futures. This will
% reduce the dimensions of the unitary operator. The default tolerance for
% merging quantum memory states is delta=1/(2*sqrt(N)) but one has freedom 
% to tweak it to any value.
% 
% The following article shall be cited when this code is used:
% Ho, Matthew, Ryuji Takagi, and Mile Gu. "Enhancing quantum models of
% stochastic processes with error mitigation." arXiv preprint
% arXiv:2105.06448 (2021).
% 
% Email:
% ho.matthew.0015@gmail.com
% 
% -------------------------------------------------------------------------
% Code description:
% 
% outputs = Unitary_generation(L, bitstream)
%
% Inputs:
% 1. L = choice of history length. If L is not given, the algorithm will
%        choose an L based on the length of bitstream.
% 2. bits = bitstream {0,1}
% 3. delta (If not given, default to 1/(2*sqrt(N) -- proof in 
%           arXiv:2105.06448).
%
% Outputs:
% out{1,2} = sprintf('Inferred Cq(L=%d) - non-merged states',L);
% out{2,2} = sprintf('Inferred Cq(L=%d) - merged states',L);
% out{3,2} = sprintf('Quantum memory states (merged)');
% out{4,2} = sprintf('Topological complexity Dq(L=%d,delta=%.6f)',L,delta);
% out{5,2} = sprintf('Unitary operator');
% out{6,2} = sprintf('Transition Matrix');
% out{7,2} = sprintf('Delta');
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Technical Logs:
% 2020-04-23 - 2247 hrs - Code is completed!
% 2020-06-10 - 1633 hrs - Inputs changed to varargin with error checks.
% 2020-06-11 - 2021 hrs - Added weights to merged conditional probabilities
%                         of futures
% 2020-07-24 - 1520 hrs - Checked for the rounding errors. No rounding
%                         errors. Rounding tolerance was already set to
%                         rounding_DP = 1E-15
% 2020-08-15 - xxxx hrs - Assigning rand numbers to only 2 rows per
%                         zero-col instead of all rows per zero-col.
%                         Relabelled new section as [METHOD 1] and old 
%                         section as [METHOD 2].
%                       - We use [METHOD 2].
% -------------------------------------------------------------------------

disp('======================== UNITARY GENERATION =========================')
disp('         CHECK INPUTS:')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT PARAMS L, BITS, AND DELTA %%%
%%% AND ERROR CHECKS                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    
    case 2
        L = varargin{1};
        bits = varargin{2};
        delta = sqrt(1/(2*length(bits)));
        fprintf('                [Inputs L = %d, %d bits, delta (default) = %.6f] \n',L,length(bits),delta)
        
    case 3
        L = varargin{1};
        bits = varargin{2};
        delta = varargin{3};
        fprintf('                [Inputs L = %d, %d bits, delta = %.6f] \n',L,length(bits),delta)
        
    otherwise
        out = NaN;
        fprintf('Error 0, check the inputs \n')
        disp('================== ERROR WITH UNITARY GENERATION ====================')
        return
        
end

% Check bit length is ok
Lmax = floor(log2(length(bits)/1000));
if Lmax <= 0
    out = NaN;
    fprintf('Error 1, check the input for data stream \n')
    disp('================== ERROR WITH UNITARY GENERATION ====================')
    return
    
elseif L < 1
    out = NaN;
    fprintf('Error 2, check the input for L \n')
    disp('================== ERROR WITH UNITARY GENERATION ====================')
    return
    
elseif delta < 0 || delta >= 1
    out = NaN;
    fprintf('Error 3, check the input for delta \n')
    disp('================== ERROR WITH UNITARY GENERATION ====================')
    return
    
end

% L possible values
if L > Lmax
    fprintf('                Input L = %d exceeds computational capabilities. Set L = Lmax = %d\n',L,Lmax)
    L = Lmax; 
    if L > 8
        L = 8;
    end
elseif L <= Lmax
    L = L;
end

fprintf('                The code will now use\n')
fprintf('                    L = %d, %d bits \n',L,length(bits))
fprintf('                    delta = %.6f \n',delta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                    %%%
%%%    THE UNITARY CODE STARTS HERE    %%%
%%%                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('         PART 1: BUILDING AND MERGING CAUSAL STATES')

tic

%%% LOOP UNTIL FAIL = 0 
%%% THIS ADJUSTS DELTA TOLERANCE
fail = 1;
while fail == 1

    %%%%%%%%%%%%%
    %%% BITSTREAM 
    %%%%%%%%%%%%%

    % Flip it to a row vector if necessary.
    [r, c] = size(bits);
    if r > 1 && c == 1
        bits = bits';
    end
    bitsSize = size(bits);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% USER INPUT LENGTH OF BITSTREAM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    len = L; %length of bitstream

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TO CONSTRUCT MATRIX WITH INCREASING TIMESTEPS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    bitsMatrix = zeros(bitsSize(2)-(2*len)+len,(len+1)); % Edited on 24-Mar-2018

    for i=1:1:bitsSize(2)-(2*len)+len  % Edited on 24-Mar-2018
        bitsMatrix(i,:) = bits(i:(len+1)+i-1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CONDITIONAL PROBABILITIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    condCounts = zeros(2^len,2);
    for i=1:1:size(bitsMatrix,1)
        r = bin2dec(num2str(bitsMatrix(i,1:len)))+1;
        c = bitsMatrix(i,len+1)+1;
        condCounts(r,c) = condCounts(r,c)+1;
    end

    condCounts;
    condProb = zeros(2^len,2);
    for i=1:1:2^len
        for j=1:1:2
            condProb(i,j) = condCounts(i,j)/sum(condCounts(i,:));
        end
    end

    for i=1:1:2^len
        for j=1:1:2
            if isnan(condProb(i,j)) == 1
                condProb(i,j) = 0;
            end
        end
    end
    condProb;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TO GENERATE MATRIX OF SIZE (2^L, L)
    %%% EACH ROW CORRESPONDS TO THE BINARY 
    %%% REPRESENTATION OF EACH ROW NUMBER+1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generating the past and future combinations
    base = 2;
    maxdecimal = base^len;
    nums = linspace(1,maxdecimal,maxdecimal);
    past = [];
    for i=1:1:length(nums)
        str = dec2base(nums(i)-1,base);
        vector = str -'0';
        while length(vector)<len
            vector = [0 vector];
        end
        past = vertcat(past,vector);
    end
    past;
    future = past;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% EVALUATING CONDITIONAL PROBABILITIES - CONCATENATED %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    condprobs_concatenated = zeros(2^len,2^len);
    %tic
    for i=1:1:2^len
        for j=1:1:2^len

            % step (1)
            temp = [past(i,:) future(j,:)];
    
            % step (2)
            mat = zeros(len,len+1);
            for k=1:1:len
                mat(k,:) = temp(k:k+len);
            end
            mat;
    
            % step (3)
            for k=1:1:len
                row_to_access = bin2dec(num2str(mat(k,1:len)))+1;
                col_to_access = mat(k,len+1)+1;
                prob_from_mat(k) = condProb(row_to_access,col_to_access);
            end
    
            condprobs_concatenated(i,j) = prod(prob_from_mat);
    
        end
    end
    condprobs_concatenated = condprobs_concatenated'; % FLIP TO COLUMN-WISE STATES
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% QUANTUM MEMORY STATES %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    probamps = sqrt(condprobs_concatenated); % COLUMN-WISE STATES
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% BUILDING / MERGING STATES %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sigma_states = zeros(1,size(condprobs_concatenated,1));
    for i=1:1:length(sigma_states)
        if sum(condprobs_concatenated(i,:)) > 0
            sigma_states(i) = max(sigma_states)+10;
        end
    end
    sigma_states;

    indices = [];
    for i=1:1:length(sigma_states)
        if sum(condprobs_concatenated(i,:)) > 0
            indices = [indices, i];
        end
    end
    indices;
    
    %%% MERGE STATES USING INNER PRODUCTS <sigma_i|sigma_j> 
    % If <o|o> >= 1-delta, same state, else different state
    for i=1:1:length(indices)-1
        for j=i+1:1:length(indices)
            if probamps(:,indices(i))'*probamps(:,indices(j)) >= 1-delta
    %             fprintf('i=%d,j=%d, <o|o> = %.4f \n',indices(i),indices(j),initialStates(:,indices(i))'*initialStates(:,indices(j)))
                sigma_states(indices(j)) = sigma_states(indices(i));
            end
        end
    end
    sigma_states;
    
    % Reducing to consecutive state numbers
    stateno = zeros(length(sigma_states));
    for i=1:1:max(max(sigma_states))/10
        x = find(sigma_states==i*10);
        if isempty(x) == 0
            if length(x) > 1
                temp = max(max(stateno))+1;
                for j=1:1:length(x)
                    stateno(x(j),:) = temp;
                end
            elseif length(x) == 1
                stateno(x,:) = max(max(stateno))+1;
            end
        end
    end
    stateno';
    sigma_states = stateno(:,1)';
    sigma_states; % After merging states
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% STATE TO STATE TRANSITIONS -- GET TRANSITION MATRIX %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is needed because now we have probability of outputting 0 and 1 but
    % >1 state to transition to. We only need 1 state to transition.
    % Extract from Cmu_generate code since it works already.

    %%% To parse through bitstream again, now to see state-to-state transitions
    %fprintf('Getting transition matrix... ')
    %tic
    bitsMatrix2 = zeros(bitsSize(2)-2*L+1,2*L); %no. of cols = 2*L
    for i=1:1:bitsSize(2)-2*L+1
        bitsMatrix2(i,:) = bits(i:2*L+i-1);
    end
    state_trans = zeros(length(bitsMatrix2),1);

    stateno';
    for i=1:1:length(bitsMatrix2)
        % Build state_trans
        idx_r = bin2dec(num2str(bitsMatrix2(i,1:L)))+1;     % past
        idx_c = bin2dec(num2str(bitsMatrix2(i,L+1:2*L)))+1; % future
        state_trans(i) = stateno(idx_r,idx_c);              % state-to-state transitions in a time series format
    end
    state_trans_without_zeros = state_trans(state_trans ~= 0); % This does what the previous 10 lines of code does.
    state_trans = state_trans_without_zeros; % I just made it the same variable as previously done.

    % Transition matrix!!
    Tcount = zeros(max(max(stateno)));
    for i=1:1:length(state_trans)-1
        row = state_trans(i);    % past state
        col = state_trans(i+1);  % future state
        Tcount(row,col) = Tcount(row,col)+1; %Count the transitions from State(i) to State(i+1) given by future(i) to future(i+1)
    end
    Tcount';

    % Past represented by row numbers, futures represented by column numbers.
    % Transition probabilities. Past = rows, Futures = columns
    [r, c] = size(Tcount);
    T = zeros(r,c);
    for i=1:1:r
        for j=1:1:c
            T(i,j) = Tcount(i,j)/sum(Tcount(i,:));
        end
    end       
    T; 
    for i=1:1:length(T)
        for j=1:1:length(T)
            if isnan(T(i,j)) == 1
                T(i,j) = 0; %Get rid of NaNs before diagonalising
            end
        end
    end
    transition_matrix = T'; % Each col is idx of a state, each row is idx of next state, elements are probabilities of transitions

    % TECHNICAL NOTE: IF STATE TRANSITS TO >2 OTHER STATES, IT IMPLIES
    % NON-UNIFILARITY. THEREFORE, ADJUST DELTA UNTIL EACH STATE TRANSITS TO <= 2
    % OTHER STATES.

    %%% Check columns of T for >2 entries per column. If true, repeat but
    %%% increase delta.

    entry_count = zeros(1,size(T,2));
    for j=1:1:size(T,2)
        for i=1:1:size(T,1)
            if T(i,j) ~= 0
                entry_count(i) = entry_count(i)+1;
            end
        end
    end
    entry_count;

    if max(entry_count) > 2
        fail = 1;
        olddelta = delta;
        delta = delta + sqrt(1/length(bits));
        fprintf('                 Initial delta: %.6f failed. Increasing delta to %.6f \n',olddelta,delta)
    else
        fail = 0;
        fprintf('                 Final delta used: %.6f \n',delta)
    end


end
fprintf('                 Number of causal states: %d \n',max(sigma_states))

% AFTER DELTA HAS BEEN SORTED, PROCEED.
% fprintf('TRANSITION MATRIX SETTLED \n')

sigma_states;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FINDING PROBABILITY OF EACH UNMERGED STATE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bitsMatrix3 = zeros(bitsSize(2)-L+1,L);
for i=1:1:bitsSize(2)-len+1
    bitsMatrix3(i,:) = bits(i:L+i-1);
end
[r6, ~] = size(bitsMatrix3);

stateProbCount = zeros(2^L,1);
for h=1:1:r6
    r = bin2dec(num2str(bitsMatrix3(h,:))) + 1;
    stateProbCount(r) = stateProbCount(r)+1;
end
% Calculating probability vector
stateProbCount;
stateProbVec = stateProbCount/sum(stateProbCount);
stateProbVec = stateProbVec';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SQRTED CONDITIONAL FUTURES FOR THE STATES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIND WEIGHTED SQRT COND FUTS - 2020-06-11, 2021 hrs
independent_sigma_states = []; 
condprobs_concatenated_independent_sigma = zeros(2^L,max(sigma_states));
for i=1:1:max(sigma_states)
    independent_sigma_states(i) = i;
    vec = find(sigma_states == i);
    ttl_prob_for_that_state = sum(stateProbVec(vec));
    temp = [];
    % Weight according to stateProbVec
    for j=1:1:length(vec)
        temp = (stateProbVec(vec(j))/ttl_prob_for_that_state) * condprobs_concatenated(:,vec(j));
        condprobs_concatenated_independent_sigma(:,i) = condprobs_concatenated_independent_sigma(:,i) + temp;
    end
end
independent_sigma_states;
condprobs_concatenated_independent_sigma;
probamps_independent_sigma = sqrt(condprobs_concatenated_independent_sigma);

% TECHNICAL NOTE:
% If collapse to future = 001, then the state will be = sigma(dec2bin(num2str('001'))+1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FINDING PROBABILITY OF EACH STATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bitsMatrix4 = zeros(bitsSize(2)-L+1,L);
for i=1:1:bitsSize(2)-L+1
    bitsMatrix4(i,:) = bits(i:L+i-1);
end
[r8, ~] = size(bitsMatrix4);

stateProbCount = zeros(1,2^L);
for h=1:1:r8
    r9 = bin2dec(num2str(bitsMatrix4(h,:))) + 1;
    stateProbCount(r9) = stateProbCount(r9)+1;
end

stateProbCount;
stateProbVec = stateProbCount/sum(stateProbCount);
prob_of_independent_sigma = zeros(1,max(sigma_states));
for i=1:1:max(sigma_states)
    vec2 = find(sigma_states == i);
    vec3 = [];
    for j=1:1:length(vec2)
        coltoextract2 = vec2(j);
        vec3(j) = stateProbVec(coltoextract2);
    end
    vec3;
    prob_of_independent_sigma(i) = sum(vec3);%stateProbVec(coltoextract2);
end
prob_of_independent_sigma; % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DENSITY MATRIX OF MERGED STATES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho_merged_states = zeros(2^L,2^L);
for i=1:1:max(sigma_states)
    rho_merged_states = rho_merged_states + prob_of_independent_sigma(i) * probamps_independent_sigma(:,i) * probamps_independent_sigma(:,i)';
end
rho_merged_states;

% Calculating Cq using initial \rho
[~, eigval_merged] = eig(rho_merged_states);
log2eigval_merged = real(log2(eigval_merged));  
eigval_merged;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATING CQ - MERGED STATES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r4, ~] = size(log2eigval_merged);
for i=1:1:r4
    for j=1:1:r4
        if log2eigval_merged(i,j) == -Inf
            log2eigval_merged(i,j) = 0;
        end
    end
end
Cq_merged = abs(-trace(eigval_merged*log2eigval_merged));

%%%%%%%%%%%%%%%%%%%%
%%% GRAM SCHMIDT %%%
%%%%%%%%%%%%%%%%%%%%
% Since condfut_independent doesn't have empty columns, do a normal Gram
% Schmidt to obtain the basis wrt first column.

VVV = probamps_independent_sigma;
[n, k] = size(VVV);
UUU = zeros(n,k);
UUU(:,1) = VVV(:,1)/sqrt(VVV(:,1)'*VVV(:,1));
for i = 2:k
    UUU(:,i) = VVV(:,i);
    for j = 1:i-1
        UUU(:,i) = UUU(:,i) - ( UUU(:,i)'*UUU(:,j) )/( UUU(:,j)'*UUU(:,j) )*UUU(:,j);
    end
    UUU(:,i) = UUU(:,i)/sqrt(UUU(:,i)'*UUU(:,i));
end
ortho_basis_sigma = UUU;

%%%%%%%%%%%%%%%%%%
%%% Cij MATRIX %%% -- To use L-future timesteps or 1-future timestep?
%%%%%%%%%%%%%%%%%%

% L future timesteps (Each quantum state is: probamps_independent_sigma)
Cij = zeros(size(probamps_independent_sigma,2),size(probamps_independent_sigma,2));
for i=1:1:size(probamps_independent_sigma,2)
    for j=1:1:size(probamps_independent_sigma,2)
        Cij(i,j) = probamps_independent_sigma(:,i)' * probamps_independent_sigma(:,j);
    end
end
Cij;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FINDING COEFFICIENTS USING \ or mldivide %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coefficients = zeros(size(probamps_independent_sigma,2),size(probamps_independent_sigma,2));
forwardslash_vec = [];

for idx = 1:1:size(probamps_independent_sigma,2)

    LHS = [probamps_independent_sigma(:,1:idx)];
    RHS = ortho_basis_sigma(:,idx);


    AugmentedMatrix = [LHS, RHS];            % Augmented matrix
    [R,~] = rref(AugmentedMatrix);           % To obtain row-reduced echelon form
    coeff_vec = R(:,end)';                   % There will always be 2^L futures but may be >= 2^L states.
                                             % So, remove #ofstates+1:end of coeff_vec
    coeff_vec(max(sigma_states)+1:end) = []; % Then assign to coefficients matrix
    coefficients(idx,:) = coeff_vec;

    % MAYBE MORE EFFICIENT TO USE BACKSLASH \ OPERATOR INSTEAD OF RREF TO
    % SOLVE SYSTEM OF LINEAR EQUATIONS
    % LHS\RHS == mldivide(LHS,RHS)
    forwardslash_vec = [LHS\RHS]';

    coefficients(idx,1:length(forwardslash_vec)) = forwardslash_vec;
end
coefficients;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SINGLE FUTURE CONDITIONAL PROBABILITIES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add up P(000|sigma_j) + P(001|sigma_j) + P(010|sigma_j) + P(011|sigma_j) to be probability of outputting 0
% Add up P(100|sigma_j) + P(101|sigma_j) + P(110|sigma_j) + P(111|sigma_j) to be probability of outputting 1
% Note: Single future = 2 outcomes only --> 2 rows.

condprobs_singletimestep = zeros(2,max(sigma_states));
for j=1:1:max(sigma_states)
    condprobs_singletimestep(1,j) = sum(condprobs_concatenated_independent_sigma(1:(2^L)/2,j));
    condprobs_singletimestep(2,j) = sum(condprobs_concatenated_independent_sigma((2^L)/2+1:2^L,j));
end
condprobs_singletimestep; % 1st row = output 0;
                        % 2nd row = output 1.
                        % 1st col = sigma_1;
                        % 2nd col = sigma_2, ... etc.
                        % matrix elements = probabilities

disp('         PART 2: COMPUTING THE INFERRED QUANTUM STATISTICAL MEMORY')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTING CQ -- 2^L STATES, NOT MERGED  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% -- COPY-PASTED FROM CQ_GENERATE_CONCATENATED

sqrt_condFuture_probs = probamps';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FINDING PROBABILITY OF EACH STATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bitsMatrix_Cq = zeros(bitsSize(2)-len+1,len);
for i=1:1:bitsSize(2)-len+1
    bitsMatrix_Cq(i,:) = bits(i:len+i-1);
end
[r6, ~] = size(bitsMatrix_Cq);

stateProbCount = zeros(2^len,1);
for h=1:1:r6
    r = bin2dec(num2str(bitsMatrix_Cq(h,:))) + 1;
    stateProbCount(r) = stateProbCount(r)+1;
end

% Calculating probability vector
stateProbCount;
stateProbVec = stateProbCount/sum(stateProbCount);

%%%%%%%%%%%%%%%%
%%% FINDING \rho
%%%%%%%%%%%%%%%%

rho = zeros(2^len,2^len);
for i=1:1:2^len
    rho = rho + stateProbVec(i) * sqrt_condFuture_probs(i,:)' * sqrt_condFuture_probs(i,:);
end
rho;

% Calculating Cq using initial \rho
[~, eigenvalues] = eig(rho);
log2eigenvalues = real(log2(eigenvalues));
eigenvalues;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATING CQ - NOT MERGED %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r4, ~] = size(log2eigenvalues);
for i=1:1:r4
    for j=1:1:r4
        if log2eigenvalues(i,j) == -Inf
            log2eigenvalues(i,j) = 0;
        end
    end
end
Cq_concatenated = abs(-trace(eigenvalues*log2eigenvalues));

% For dimensions of the system, use topological complexity
%topological_complexity = log2(count);
topological_complexity = log2(max(sigma_states));


%%% END OF COMPUTING CQ - NOT MERGED %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('         PART 3: GETTING NEXT STATES @ SINGLE TIME STEP')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GETTING NEXTSTATE_SINGLETIMESTEP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

transition_matrix = T';
condprobs_singletimestep;
probamp_singletimestep = sqrt(condprobs_singletimestep);

% ASSIGN TEMP VARIABLES
temp_transition_matrix = transition_matrix;
for idx_i = 1:1:size(temp_transition_matrix,1)
    for idx_j = 1:1:size(temp_transition_matrix,2)
        if temp_transition_matrix(idx_i,idx_j) == 0
            temp_transition_matrix(idx_i,idx_j) = NaN; % To prevent round(value,0) = 0 from coinciding with other exact 0 values. This affects finding nextstate_singletimestep.
        end
    end
end
temp_transition_matrix;
temp_condfut_singletimestep = condprobs_singletimestep;

% disp('==================================')

nextstate_singletimestep = zeros(size(temp_condfut_singletimestep,1),size(temp_condfut_singletimestep,2));
% GO THROUGH EACH COLUMN. IF EACH COLUMN FAILS, ADJUST ROUNDING_DP
for j = 1:1:size(temp_condfut_singletimestep,2)
    rounding_DP = 15;
    does_this_fail = 1;
    while does_this_fail == 1
        for i=1:1:size(temp_condfut_singletimestep,1)
            if temp_condfut_singletimestep(i,j) ~= 0 % Nonzero means it has state transitions
                
                temp_transition_matrix(:,j);
                temp_condfut_singletimestep(i,j);
                at_row_of_transition_matrix = find(round(temp_transition_matrix(:,j),rounding_DP) == round(temp_condfut_singletimestep(i,j),rounding_DP));
                if isempty(at_row_of_transition_matrix) == 1
                    nextstate_singletimestep(i,j) = NaN;

                % ADDED THIS 'ELSEIF' FOR THE CASE WHERE TRANSITION MATRIX
                % HAS COLUMNS [0.5, 0.5]', MEANING THAT SAME PROBABILITY OF
                % TRANSITIONING TO STATE 1 AND 2.
                elseif length(at_row_of_transition_matrix) == 2
                    if temp_transition_matrix(at_row_of_transition_matrix(1),j) == temp_transition_matrix(at_row_of_transition_matrix(2),j)
                        nextstate_singletimestep(i,j) = at_row_of_transition_matrix(i);
                    end
                    
                else
                    nextstate_singletimestep(i,j) = at_row_of_transition_matrix;
                end
                
            end
        end
        if sum(isnan(nextstate_singletimestep(:,j)) == 1) > 0
            does_this_fail = 1;
%             fprintf('Column %2d: rounding_DP = %2d, failed\n',j,rounding_DP) 
            rounding_DP = rounding_DP-1;
        else
            does_this_fail = 0;
%             fprintf('Column %2d: rounding_DP = %2d, passed \n',j,rounding_DP)
        end
    end
end

nextstate_singletimestep; % 1st row: output 0, 2nd row: output 1

disp('         PART 4: BUILDING THE UNITARY ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BUILDING UNITARY - MODIFIED FROM OLDER CODE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Unitary_nonzeros = zeros(max(sigma_states)*2,max(sigma_states));

for jj=1:1:max(sigma_states)
    Unitary_nonzeros_row = 0;
    for ii = 1:1:max(sigma_states)
        for output = 0:1:1

            output_idx = output+1;
            
            Unitary_nonzeros_row = Unitary_nonzeros_row+1;
            Unitary_nonzeros_col = jj;
            
            inner_matrix = zeros(ii,jj);

            for i=1:1:ii
                for j=1:1:jj

                    % coefficient for idx_i
                    row_coeff_idx_i = ii;
                    T1 = coefficients(row_coeff_idx_i,i);
                    
                    % coefficient for idx_i'
                    row_coeff_idx_iprime = jj;
                    T2 = coefficients(row_coeff_idx_iprime,j);
                    
                    % prob_amp stuff
                    row_prob_amp = output_idx;
                    col_prob_amp = j;
                    T3 = probamp_singletimestep(row_prob_amp,col_prob_amp);
                    
                    % Cij stuff
                    row_Cij = i;
                    col_Cij = nextstate_singletimestep(output_idx,j);
                    if col_Cij == 0
                        T4 = 0;
                    else
                        T4 = Cij(row_Cij,col_Cij);
                    end

                    inner_matrix(i,j) = T1 * T2 * T3 * T4;
                end
            end
            
            inner_matrix;
            Unitary_element = sum(sum(inner_matrix));
            Unitary_nonzeros(Unitary_nonzeros_row,Unitary_nonzeros_col) = Unitary_element;
            
        end
    end
end
Unitary_nonzeros;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NORMALISING Unitary_nonzeros %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:1:size(Unitary_nonzeros,2)
    if norm(Unitary_nonzeros(:,j)) ~= 0
        Unitary_nonzeros(:,j) = Unitary_nonzeros(:,j)/norm(Unitary_nonzeros(:,j));
    elseif norm(Unitary_nonzeros(:,j)) == 0
        Unitary_nonzeros(:,j) = zeros(size(Unitary_nonzeros,1),1);
    end
end

Unitary_nonzeros;



% % [METHOD 1]
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % [ASSIGN RAND NUMBERS TO 2 ROWS PER ZERO-COL]
% % %%% ASSIGN RAND NUMBERS TO ZERO-COLS IN PAIRS OF ROWS
% % %%% AND DO A GRAM SCHMIDT TO ORTHONORMALISE THEM
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('                 [New method of filling in zero-cols]')
% Unitary_notyetreorder = zeros(size(Unitary_nonzeros,1),size(Unitary_nonzeros,1));
% Unitary_notyetreorder(:,1:size(Unitary_nonzeros,2)) = Unitary_nonzeros;
% 
% 
% startingcol = size(Unitary_nonzeros,2)+1;
% for j=startingcol:1:size(Unitary_notyetreorder,2)
%     i = 2*(j-startingcol+1)-1; %1:1:size(Unitary_notyetreorder,1)
%     ii = 2*(j-startingcol+1)-1+1;
% %     fprintf('i  = %d, j=%d \n',i,j);
% %     fprintf('ii = %d, j=%d \n',ii,j);
%         
%     Unitary_notyetreorder(i,j) = +rand;
%     Unitary_notyetreorder(ii,j) = -rand;
%     
% end
% Unitary_notyetreorder;


%%% [METHOD 2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ASSIGN RAND NUMBERS TO ZERO-COLS 
%%% AND DO A GRAM SCHMIDT TO ORTHONORMALISE THEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Unitary_notyetreorder = zeros(size(Unitary_nonzeros,1),size(Unitary_nonzeros,1));
Unitary_notyetreorder(:,1:size(Unitary_nonzeros,2)) = Unitary_nonzeros;

startingcol = size(Unitary_nonzeros,2)+1;
for j=startingcol:1:size(Unitary_notyetreorder,2)
    for i=1:1:size(Unitary_notyetreorder,1)
        Unitary_notyetreorder(i,j) = rand;
    end
end
Unitary_notyetreorder;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GRAM SCHMIDT FOR THE REMAINING ZERO-COLS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nn, kk] = size(Unitary_notyetreorder);
U_temp = zeros(nn,kk);
U_temp(:,1:startingcol-1) = Unitary_nonzeros;

for i = startingcol:1:kk
    U_temp(:,i) = Unitary_notyetreorder(:,i);
    for j = 1:1:i-1
        U_temp(:,i) = U_temp(:,i) - ( U_temp(:,i)'*U_temp(:,j) )/( U_temp(:,j)'*U_temp(:,j) )*U_temp(:,j);
    end
    U_temp(:,i) = U_temp(:,i)/sqrt(U_temp(:,i)'*U_temp(:,i));
    Unitary_notyetreorder(:,i) = U_temp(:,i);
end
Unitary_notyetreorder;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% REORDER THE COLS OF UNITARY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Unitary_notyetreorder;
Unitary_zerocols = Unitary_notyetreorder(:,startingcol:end);

Unitary = zeros(size(Unitary_nonzeros,1),size(Unitary_nonzeros,1));
for i=1:1:size(Unitary_nonzeros,2)
    Unitary(:,i*2-1) = Unitary_nonzeros(:,i);
    Unitary(:,i*2) = Unitary_zerocols(:,i);
end
Unitary;

%%%%% Condition for Unitarity
%%%%% Udagger = Uinverse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GENERATING THE MERGED QUANTUM MEMORY STATES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

QMS = zeros(max(sigma_states),max(sigma_states));
QMS(1,1) = 1;
for j=2:1:max(sigma_states)
    QMS(1,j) = Cij(1,j);
end
if size(QMS,1) > 1 % To avoid error from rand process with 1 QMS only.    
    QMS(2,2) = sqrt(1-QMS(1,2)^2);

    for j=3:1:max(sigma_states) % Always start with 3
       for i=2:1:j-1
           
           vec = zeros(1,i);
           for ii=1:1:1
               vec(ii) = Cij(i,j);
           end
           for ii=2:1:i
               vec(ii) = - QMS(ii-1,i) * QMS(ii-1,j);
           end
           Denominator = QMS(i,i);
           QMS(i,j) = sum(vec)/Denominator;

       end

       vecc = QMS(1:j-1,j);
       vecc_sq = zeros(length(vecc),1);
       for k=1:1:length(vecc)
           vecc_sq(k) = vecc(k).^2;
       end
       vecc_sq;
       vecc_sq = -vecc_sq;
       QMS(j,j) = sqrt(1 + sum(vecc_sq));

    end

end
% ROUND OFF TO 4 DP?
QMS;% = round(QMS,4);

%%%%%%%%%%%%%%
%%% OUTPUT %%%
%%%%%%%%%%%%%%

% OUTPUT UNITARY
% OUTPUT QUANTUM MEMORY STATES
out{1,1} = sprintf('1. Inferred Cq(L=%d) - non-merged states',L);
out{2,1} = sprintf('2. Inferred Cq(L=%d) - merged states',L);
out{3,1} = sprintf('3. Quantum memory states (merged)');
out{4,1} = sprintf('4. Topological complexity Dq(L=%d,delta=%.6f)',L,delta);
out{5,1} = sprintf('5. Unitary operator');
out{6,1} = sprintf('6. Transition Matrix');
out{7,1} = sprintf('7. Delta');
out{8,1} = sprintf('8. Nextstate and outputs (single timesteps). 1st row: output 0, 2nd row: output 1');
out{9,1} = sprintf('9. Quantum memory states (Not merged)');


out{1,2} = Cq_concatenated;
out{2,2} = Cq_merged;
out{3,2} = QMS;
out{4,2} = topological_complexity;
out{5,2} = Unitary;
out{6,2} = transition_matrix;
out{7,2} = delta;
out{8,2} = nextstate_singletimestep;
out{9,2} = probamps;

fprintf('         Elapsed time %.6f seconds. \n',toc)
disp('==================== END OF UNITARY GENERATION ======================')
                        
end









