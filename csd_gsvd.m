function [U1, U2, V1, V2, C, S] = csd_gsvd(Unitary)
%
% csd_gsvd computes the gsvd for a unitary matrix using four ways. Then,
% by comparing how close the reconstructed unitary is with respect to the
% input unitary, the code chooses the best way.
%
% [U1, U2, V1, V2, C, S] = csd_gsvd(Unitary)
%
% Unitary = [g11 g12] = Left * Mid * Right
%           [g21 g22] 
%
% where
% 
% Left = [U1   ],  Mid = [C -S],  Right = [V1   ]
%        [   U2]         [S  C]           [   V2]
%
% Input: Unitary matrix
% Output: [U1, U2, V1, V2, C, S]
%--------------------------------------------------------------------------
% Written by Matthew Ho on 2020-08-29, 2049 hrs
%--------------------------------------------------------------------------


% %%%%%%%%%%%% GENERATE RANDOM UNITARY TO TEST
% %%
% n=4;
% X = (randn(n))/sqrt(2);
% [Q,R] = qr(X);
% R = diag(diag(R)./abs(diag(R)));
% Unitary = Q*R

g11 = Unitary(1:end/2,1:end/2);
g12 = Unitary(1:end/2,end/2+1:end);
g21 = Unitary(end/2+1:end,1:end/2);
g22 = Unitary(end/2+1:end,end/2+1:end);


% disp('METHOD 1')

%%% METHOD 1
[V1_1, Csquared_1] = eig(g11'*g11);
V1_1 = V1_1';
C_1 = sqrt(Csquared_1);
S_1 = sqrt(eye(size(Csquared_1,2))-Csquared_1);

U1_1 = g11 * inv(V1_1) * inv(C_1);
U2_1 = g21 * inv(V1_1) * inv(S_1);
V2_1 = inv(C_1) * inv(U2_1) * g22;
V2_1_2 = - inv(S_1) * inv(U1_1) * g12;

% disp('METHOD 2')

%%% METHOD 2
[V2_2, Csquared_2] = eig(g22'*g22);
V2_2 = V2_2';
C_2 = sqrt(Csquared_2);
S_2 = sqrt(eye(size(Csquared_2,2))-Csquared_2);

U2_2 = g22 * inv(V2_2) * inv(C_2);
U1_2 = -g12 * inv(V2_2) * inv(S_2);
V1_2 = inv(C_2) * inv(U1_2) * g11;
V1_2_2 = inv(S_2) * inv(U2_2) * g21;

% disp('METHOD 3')

%%% METHOD 3
[V1_3, Ssquared_3] = eig(g21'*g21);
V1_3 = V1_3';
S_3 = sqrt(Ssquared_3);
C_3 = sqrt(eye(size(Ssquared_3,2))-Ssquared_3);

U2_3 = g21 * inv(V1_3) * inv(S_3);
U1_3 = g11 * inv(V1_3) * inv(C_3);
V2_3 = inv(C_3) * inv(U2_3) * g22;
V2_3_2 = -inv(S_3) * inv(U1_3) * g12;

% disp('METHOD 4')

%%% METHOD 4
[V2_4, Ssquared_4] = eig(g12'*g12);
V2_4 = V2_4';
S_4 = sqrt(Ssquared_4);
C_4 = sqrt(eye(size(Ssquared_4,2))-Ssquared_4);

U1_4 = -g12 * inv(V2_4) * inv(S_4);
U2_4 = g22 * inv(V2_4) * inv(C_4);
V1_4 = inv(C_4) * inv(U1_4) * g11;
V1_4_2 = inv(S_4) * inv(U2_4) * g21;



%%% DISPLAY RESULTS
[U1_1, NaN(size(U1_1,1),1), U2_1, NaN(size(U1_1,1),1), V1_1, NaN(size(U1_1,1),1), V2_1, NaN(size(U1_1,1),1), C_1, NaN(size(U1_1,1),1), S_1];
[U1_2, NaN(size(U1_1,1),1), U2_2, NaN(size(U1_1,1),1), V1_2, NaN(size(U1_1,1),1), V2_2, NaN(size(U1_1,1),1), C_2, NaN(size(U1_1,1),1), S_2];
[U1_3, NaN(size(U1_1,1),1), U2_3, NaN(size(U1_1,1),1), V1_3, NaN(size(U1_1,1),1), V2_3, NaN(size(U1_1,1),1), C_3, NaN(size(U1_1,1),1), S_3];
[U1_4, NaN(size(U1_1,1),1), U2_4, NaN(size(U1_1,1),1), V1_4, NaN(size(U1_1,1),1), V2_4, NaN(size(U1_1,1),1), C_4, NaN(size(U1_1,1),1), S_4];


LEFT_1 = [U1_1, zeros(size(U1_1,1)); zeros(size(U1_1,1)), U2_1];
MID_1 = [C_1 -S_1; S_1 C_1];
RIGHT_1 = [V1_1, zeros(size(U1_1,1)); zeros(size(U1_1,1)) V2_1];
Result1 = LEFT_1 * MID_1 * RIGHT_1;

LEFT_2 = [U1_2, zeros(size(U1_2,1)); zeros(size(U1_2,1)), U2_2];
MID_2 = [C_2 -S_2; S_2 C_2];
RIGHT_2 = [V1_2, zeros(size(U1_2,1)); zeros(size(U1_2,1)) V2_2];
Result2 = LEFT_2 * MID_2 * RIGHT_2;

LEFT_3 = [U1_3, zeros(size(U1_3,1)); zeros(size(U1_3,1)), U2_3];
MID_3 = [C_3 -S_3; S_3 C_3];
RIGHT_3 = [V1_3, zeros(size(U1_3,1)); zeros(size(U1_3,1)) V2_3];
Result3 = LEFT_3 * MID_3 * RIGHT_3;

LEFT_4 = [U1_4, zeros(size(U1_4,1)); zeros(size(U1_4,1)), U2_4];
MID_4 = [C_4 -S_4; S_4 C_4];
RIGHT_4 = [V1_4, zeros(size(U1_4,1)); zeros(size(U1_4,1)) V2_4];
Result4 = LEFT_4 * MID_4 * RIGHT_4;


results_mat(:,1) = 1:1:4;
results_mat(1,2) = abs(sum(sum(Result1-Unitary)));
results_mat(2,2) = abs(sum(sum(Result2-Unitary)));
results_mat(3,2) = abs(sum(sum(Result3-Unitary)));
results_mat(4,2) = abs(sum(sum(Result4-Unitary)));

fprintf('Mtd 1: %.32f \nMtd 2: %.32f \nMtd 3: %.32f \nMtd 4: %.32f \n',results_mat(1,2),results_mat(2,2),results_mat(3,2),results_mat(4,2))

results_mat;

idx = find(results_mat(:,2) == min(results_mat(:,2)));
fprintf('Pick method %d \n\n',idx)

if idx == 1
    U1 = U1_1;
    U2 = U2_1;
    V1 = V1_1;
    V2 = V2_1;
    C = C_1;
    S = S_1;
elseif idx == 2
    U1 = U1_2;
    U2 = U2_2;
    V1 = V1_2;
    V2 = V2_2;
    C = C_2;
    S = S_2;
elseif idx == 3
    U1 = U1_3;
    U2 = U2_3;
    V1 = V1_3;
    V2 = V2_3;
    C = C_3;
    S = S_3;
elseif idx == 4
    U1 = U1_4;
    U2 = U2_4;
    V1 = V1_4;
    V2 = V2_4;
    C = C_4;
    S = S_4;
end


[U1, U2, V1, V2, C, S];


end












