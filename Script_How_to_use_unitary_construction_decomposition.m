

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% EXAMPLE OF HOW TO USE Unitary_generation.m AND csd_gsvd.m FUNCTIONS 
%%% 1. GENERATE STOCHASTIC PROCESS -- EG. PERTURBED COIN PROCESS
%%% 2. CALL Unitary_generate
%%% 3. CALL csd_gsvd
%%% 4. CALL ZYZ_decomposition
%%% 5. PRINT FOR USE (ON QISKIT)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
clc

% disp('PERTURBED COIN')
% PERTURBED COIN

p = 0.2;
q = p;
tic
bitlength = 10000;
cut = 10000;
bits = zeros(1,bitlength+cut);
state = round(rand)+1;
for i=1:1:length(bits)
    if state == 1
        if rand < p
            bits(i) = 1;
            state = 2;
        else
            bits(i) = 0;
            state = 1;
        end
    elseif state == 2
        if rand < q
            bits(i) = 0;
            state = 1;
        else
            bits(i) = 1;
            state = 2;
        end
    end
end
toc
bits(1:cut) = [];
processtitletext = sprintf('Perturbed Coin Process');
fprintf('Perturbed Coin %d bits generated \n',length(bits))

L=2;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CONSTRUCT UNITARY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

delta = 1/(2*sqrt(bitlength)); % Recommend this, proof in arXiv:2105.06448

output_unitary = Unitary_generation(L,bits,delta) 
Unitary = output_unitary{5,2}

%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DECOMPOSE UNITARY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

[U1,U2,V1,V2,COSINES,SINES] = csd_gsvd(Unitary);

MID = [COSINES, -SINES; SINES, COSINES];

NOT = [0 1 ; 1 0];
CNOT = [eye(2), zeros(2); zeros(2), NOT];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BREAK INTO ZYZ ROTATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_U1 = ZYZ_decomposition(U1);
U1_phase  = output_U1{2,1}(1);
U1_phi    = output_U1{2,1}(2);
U1_theta  = output_U1{2,1}(3);
U1_lamb   = output_U1{2,1}(4);

CU1_phase_gate  = [eye(2), zeros(2); zeros(2), output_U1{2,2}];
CU1_phi_gate    = [eye(2), zeros(2); zeros(2), output_U1{2,3}];
CU1_theta_gate  = [eye(2), zeros(2); zeros(2), output_U1{2,4}];
CU1_lamb_gate   = [eye(2), zeros(2); zeros(2), output_U1{2,5}];

A = CNOT * CU1_phase_gate * CNOT * CU1_phase_gate * CU1_phi_gate * CU1_theta_gate * CU1_lamb_gate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_U2 = ZYZ_decomposition(U2);
U2_phase  = output_U2{2,1}(1);
U2_phi    = output_U2{2,1}(2);
U2_theta  = output_U2{2,1}(3);
U2_lamb   = output_U2{2,1}(4);

CU2_phase_gate  = [eye(2), zeros(2); zeros(2), output_U2{2,2}];
CU2_phi_gate    = [eye(2), zeros(2); zeros(2), output_U2{2,3}];
CU2_theta_gate  = [eye(2), zeros(2); zeros(2), output_U2{2,4}];
CU2_lamb_gate   = [eye(2), zeros(2); zeros(2), output_U2{2,5}];

B = CNOT * CU2_phase_gate * CNOT * CU2_phase_gate * CU2_phi_gate * CU2_theta_gate * CU2_lamb_gate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_V1 = ZYZ_decomposition(V1);
V1_phase  = output_V1{2,1}(1);
V1_phi    = output_V1{2,1}(2);
V1_theta  = output_V1{2,1}(3);
V1_lamb   = output_V1{2,1}(4);

CV1_phase_gate  = [eye(2), zeros(2); zeros(2), output_V1{2,2}];
CV1_phi_gate    = [eye(2), zeros(2); zeros(2), output_V1{2,3}];
CV1_theta_gate  = [eye(2), zeros(2); zeros(2), output_V1{2,4}];
CV1_lamb_gate   = [eye(2), zeros(2); zeros(2), output_V1{2,5}];

C = CNOT * CV1_phase_gate * CNOT * CV1_phase_gate * CV1_phi_gate * CV1_theta_gate * CV1_lamb_gate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_V2 = ZYZ_decomposition(V2);
V2_phase  = output_V2{2,1}(1);
V2_phi    = output_V2{2,1}(2);
V2_theta  = output_V2{2,1}(3);
V2_lamb   = output_V2{2,1}(4);

CV2_phase_gate  = [eye(2), zeros(2); zeros(2), output_V2{2,2}];
CV2_phi_gate    = [eye(2), zeros(2); zeros(2), output_V2{2,3}];
CV2_theta_gate  = [eye(2), zeros(2); zeros(2), output_V2{2,4}];
CV2_lamb_gate   = [eye(2), zeros(2); zeros(2), output_V2{2,5}];

D = CNOT * CV2_phase_gate * CNOT * CV2_phase_gate * CV2_phi_gate * CV2_theta_gate * CV2_lamb_gate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

recon_unitary = kron(NOT,eye(2)) * A * kron(NOT,eye(2)) * B * MID * kron(NOT,eye(2)) * C * kron(NOT,eye(2)) * D;

Unitary_gateset = zeros(4,4,33);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% In order of application on a quantum computer %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Unitary_gateset(:,:,1) = CV2_lamb_gate;
Unitary_gateset(:,:,2) = CV2_theta_gate;
Unitary_gateset(:,:,3) = CV2_phi_gate;
Unitary_gateset(:,:,4) = CV2_phase_gate;
Unitary_gateset(:,:,5) = CNOT;
Unitary_gateset(:,:,6) = CV2_phase_gate;
Unitary_gateset(:,:,7) = CNOT;
Unitary_gateset(:,:,8) = kron(NOT,eye(2));

Unitary_gateset(:,:,9)  = CV1_lamb_gate;
Unitary_gateset(:,:,10) = CV1_theta_gate;
Unitary_gateset(:,:,11) = CV1_phi_gate;
Unitary_gateset(:,:,12) = CV1_phase_gate;
Unitary_gateset(:,:,13) = CNOT;
Unitary_gateset(:,:,14) = CV1_phase_gate;
Unitary_gateset(:,:,15) = CNOT;
Unitary_gateset(:,:,16) = kron(NOT,eye(2));

Unitary_gateset(:,:,17) = MID;

Unitary_gateset(:,:,18) = CU2_lamb_gate;
Unitary_gateset(:,:,19) = CU2_theta_gate;
Unitary_gateset(:,:,20) = CU2_phi_gate;
Unitary_gateset(:,:,21) = CU2_phase_gate;
Unitary_gateset(:,:,22) = CNOT;
Unitary_gateset(:,:,23) = CU2_phase_gate;
Unitary_gateset(:,:,24) = CNOT;
Unitary_gateset(:,:,25) = kron(NOT,eye(2));

Unitary_gateset(:,:,26) = CU1_lamb_gate;
Unitary_gateset(:,:,27) = CU1_theta_gate;
Unitary_gateset(:,:,28) = CU1_phi_gate;
Unitary_gateset(:,:,29) = CU1_phase_gate;
Unitary_gateset(:,:,30) = CNOT;
Unitary_gateset(:,:,31) = CU1_phase_gate;
Unitary_gateset(:,:,32) = CNOT;
Unitary_gateset(:,:,33) = kron(NOT,eye(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECK IF THE DECOMPOSED ELEMENTARY QUANTUM GATES GIVE THE ORIGINAL UNITARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TEST = eye(4);
for i=1:1:33
    TEST = TEST * Unitary_gateset(:,:,33-i+1);
end

fprintf('CHECK: Reconstructing the unitary from the elementary quantum gates \n\n')
disp(real(TEST))

fprintf('Original unitary (before decomposition into elementary quantum gates) \n\n')
disp(Unitary)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXAMPLE OF 2X2 UNITARY OPERATOR'S GATES ON QISKIT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('## FOR APPLICATION ON QISKIT (2X2 UNITARY OPERATOR): \n')
fprintf('## --------------- COPY FROM THIS LINE ---------------------------\n')
fprintf('# Input angles: \n')
fprintf('U1_phase_angle  = %20.16f \n',U1_phase)
fprintf('U1_RZGate_phi   = %20.16f \n',U1_phi)
fprintf('U1_RYGate_theta = %20.16f \n',U1_theta)
fprintf('U1_RZGate_lamb  = %20.16f \n',U1_lamb)
fprintf('U2_phase_angle  = %20.16f \n',U2_phase)
fprintf('U2_RZGate_phi   = %20.16f \n',U2_phi)
fprintf('U2_RYGate_theta = %20.16f \n',U2_theta)
fprintf('U2_RZGate_lamb  = %20.16f \n',U2_lamb)
%
angle1 = atan((MID(3,1))/(MID(1,1)));
angle2 = atan((MID(4,2))/(MID(2,2)));
fprintf('MID_angle1      = %20.16f \n',angle1)
fprintf('MID_angle2      = %20.16f \n',angle2)
%
fprintf('V1_phase_angle  = %20.16f \n',V1_phase)
fprintf('V1_RZGate_phi   = %20.16f \n',V1_phi)
fprintf('V1_RYGate_theta = %20.16f \n',V1_theta)
fprintf('V1_RZGate_lamb  = %20.16f \n',V1_lamb)
fprintf('V2_phase_angle  = %20.16f \n',V2_phase)
fprintf('V2_RZGate_phi   = %20.16f \n',V2_phi)
fprintf('V2_RYGate_theta = %20.16f \n',V2_theta)
fprintf('V2_RZGate_lamb  = %20.16f \n',V2_lamb)

fprintf('multiplex = UnitaryGate(np.array([[np.cos(MID_angle1), 0,                 -np.sin(MID_angle1), 0                  ], \n')
fprintf('                                  [0,                  np.cos(MID_angle2), 0,                  -np.sin(MID_angle2)], \n')
fprintf('                                  [np.sin(MID_angle1), 0,                  np.cos(MID_angle1), 0                  ], \n')
fprintf('                                  [0,                  np.sin(MID_angle2), 0,                   np.cos(MID_angle2)]])) \n\n')

fprintf('# Apply to quantum circuit ''circ'': \n')
fprintf('#circ.ry(1.287002217587,0)# Rotate to get |sigma_2> \n')
fprintf('circ.crz(V2_RZGate_lamb,1,0) \n')
fprintf('circ.cry(V2_RYGate_theta,1,0) \n')
fprintf('circ.crz(V2_RZGate_phi,1,0) \n')
fprintf('circ.cu1(V2_phase_angle,1,0) \n')
fprintf('circ.cx(1,0) \n')
fprintf('circ.cu1(V2_phase_angle,1,0) \n')
fprintf('circ.cx(1,0) \n')
fprintf('circ.x(1) \n')
fprintf('circ.crz(V1_RZGate_lamb,1,0) \n')
fprintf('circ.cry(V1_RYGate_theta,1,0) \n')
fprintf('circ.crz(V1_RZGate_phi,1,0) \n')
fprintf('circ.cu1(V1_phase_angle,1,0) \n')
fprintf('circ.cx(1,0) \n')
fprintf('circ.cu1(V1_phase_angle,1,0) \n')
fprintf('circ.cx(1,0) \n')
fprintf('circ.x(1) \n')
fprintf('circ.unitary(multiplex, [0,1], label=''GSVD MTPLX'') \n')
fprintf('circ.crz(U2_RZGate_lamb,1,0) \n')
fprintf('circ.cry(U2_RYGate_theta,1,0) \n')
fprintf('circ.crz(U2_RZGate_phi,1,0) \n')
fprintf('circ.cu1(U2_phase_angle,1,0) \n')
fprintf('circ.cx(1,0) \n')
fprintf('circ.cu1(U2_phase_angle,1,0) \n')
fprintf('circ.cx(1,0) \n')
fprintf('circ.x(1) \n')
fprintf('circ.crz(U1_RZGate_lamb,1,0) \n')
fprintf('circ.cry(U1_RYGate_theta,1,0) \n')
fprintf('circ.crz(U1_RZGate_phi,1,0) \n')
fprintf('circ.cu1(U1_phase_angle,1,0) \n')
fprintf('circ.cx(1,0) \n')
fprintf('circ.cu1(U1_phase_angle,1,0) \n')
fprintf('circ.cx(1,0) \n')
fprintf('circ.x(1) \n')

fprintf('## --------------- TO THIS LINE ----------------------------------\n')







%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECK IF THE ANGLES FOR THE DECOMPOSED ELEMENTARY QUANTUM GATES GIVE THE ORIGINAL UNITARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U1_phase_gate = phasegate(U1_phase);
U1_global_phase = [0, 1; 1, 0] * U1_phase_gate * [0, 1; 1, 0] * U1_phase_gate;
U1_Zphi   = RZ(U1_phi);
U1_Ytheta = RY(U1_theta);
U1_Zlamb  = RZ(U1_lamb);

U1_reconstr_frm_gates = U1_global_phase * U1_Zphi * U1_Ytheta * U1_Zlamb;
U1;

U2_phase_gate = phasegate(U2_phase);
U2_global_phase = [0, 1; 1, 0] * U2_phase_gate * [0, 1; 1, 0] * U2_phase_gate;
U2_Zphi   = RZ(U2_phi);
U2_Ytheta = RY(U2_theta);
U2_Zlamb  = RZ(U2_lamb);

U2_reconstr_frm_gates = U2_global_phase * U2_Zphi * U2_Ytheta * U2_Zlamb;
U2;

V1_phase_gate = phasegate(V1_phase);
V1_global_phase = [0, 1; 1, 0] * V1_phase_gate * [0, 1; 1, 0] * V1_phase_gate;
V1_Zphi   = RZ(V1_phi);
V1_Ytheta = RY(V1_theta);
V1_Zlamb  = RZ(V1_lamb);

V1_reconstr_frm_gates = V1_global_phase * V1_Zphi * V1_Ytheta * V1_Zlamb;
V1;

V2_phase_gate = phasegate(V2_phase);
V2_global_phase = [0, 1; 1, 0] * V2_phase_gate * [0, 1; 1, 0] * V2_phase_gate;
V2_Zphi   = RZ(V2_phi);
V2_Ytheta = RY(V2_theta);
V2_Zlamb  = RZ(V2_lamb);

V2_reconstr_frm_gates = V2_global_phase * V2_Zphi * V2_Ytheta * V2_Zlamb;
V2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECK IF THE ANGLES CAN RECONSTRUCT THE 4x4 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XI = kron([0 1; 1 0],eye(2));

%%% Reconstructing U1, U2, V1, V2, MID
recon_left_U1 = eye(4);
recon_left_U1(3:end,3:end) = U1_reconstr_frm_gates;
recon_left_U2 = eye(4);
recon_left_U2(3:end,3:end) = U2_reconstr_frm_gates;
recon_left = XI * recon_left_U1 * XI * recon_left_U2;
recon_mid = [COSINES -SINES; SINES COSINES];
recon_right_V1 = eye(4);
recon_right_V1(3:end,3:end) = V1_reconstr_frm_gates;
recon_right_V2 = eye(4);
recon_right_V2(3:end,3:end) = V2_reconstr_frm_gates;
recon_right = XI * recon_right_V1 * XI * recon_right_V2;
Unitary_reconstructed_frm_phases = recon_left * recon_mid * recon_right; % OK!


fprintf('------------------------------------------------------------------\n\n')

fprintf('CHECK: Reconstructing the unitary from the angles of the elementary quantum gates \n\n')
disp(real(Unitary_reconstructed_frm_phases))

fprintf('Original unitary (before decomposition into elementary quantum gates) \n\n')
disp(Unitary)












%%
Phasegate = phasegate(pi/2);

function Zgate = RZ(angle)
    Zgate = [exp(-1i*angle/2) 0; 0 exp(1i*angle/2)];
end

function Ygate = RY(angle)
    Ygate = [cos(angle/2) -sin(angle/2); sin(angle/2) cos(angle/2)];
end

function Phasegate = phasegate(angle)
    Phasegate = [1 0; 0 exp(1i*angle)];
end

function gate = mtplx4(angle1, angle2)
    gate = [cos(angle1/2), 0, sin(angle1/2), 0; 0, cos(angle2/2), 0, sin(angle2/2); -sin(angle1/2), 0, cos(angle1/2), 0; 0, -sin(angle2/2), 0, cos(angle2/2)];
end

function gate = mtplx8(angle1 ,angle2, angle3, angle4)
    gate = [cos(angle1/2), 0, 0, 0, sin(angle1/2), 0, 0, 0;             0, cos(angle2/2), 0, 0, 0, sin(angle2/2), 0, 0;             0, 0, cos(angle3/2), 0, 0, 0, sin(angle3/2), 0;             0, 0, 0, cos(angle4/2), 0, 0, 0, sin(angle4/2);             -sin(angle1/2), 0, 0, 0, cos(angle1/2), 0, 0, 0;             0, -sin(angle2/2), 0, 0, 0, cos(angle2/2), 0, 0;             0, 0, -sin(angle3/2), 0, 0, 0, cos(angle3/2), 0;             0, 0, 0, -sin(angle4/2), 0, 0, 0, cos(angle4/2)];
end


function output = ZYZ_decomposition(U)
%--------------------------------------------------------------------------
% To find the ZYZ rotation gates of (not special) unitary matrix mat.
% output = ZYZ_decomposition(mat)
% i.e. mat = [0,1;1,0] * phasegate(gamma) * [0,1;1,0] * phasegate(gamma) * RZ(phi) * RY(theta) * RZ(lambda);
% 
% Written by Matthew Ho on 2020-07-07, 2210 hrs.
%--------------------------------------------------------------------------

% Convert to SU(2)
coeff = 1/sqrt(det(U));
coeff = det(U)^(-0.5);

SU = coeff * U;

% Get angles
gamma = -angle(coeff); % <-- coefficient for converting to SU(2)
theta = 2 * atan2(abs(SU(2,1)),abs(SU(1,1)));
phiplamb = 2 * angle(SU(2,2));
phimlamb = 2 * angle(SU(2,1));
phi = (phiplamb + phimlamb)/2;
lamb = (phiplamb - phimlamb)/2;

su_mat_angles = [gamma, phi, theta, lamb];

% phase gate
phase_gate = [1 0 ; 0 exp(1j*gamma)];  % <-- coefficient for converting to SU(2)

% ZYZ Decomposition su_mat
RZGate_lamb = [exp(-1i*lamb/2) 0; 0 exp(1i*lamb/2)];
RYGate_theta = [cos(theta/2) -sin(theta/2); sin(theta/2) cos(theta/2)];
RZGate_phi = [exp(-1i*phi/2) 0; 0 exp(1i*phi/2)];

%[theta, phi, lamb, coeff, phase]
output{1,1} = sprintf('Angles: \n(1,1)gamma - for phase gate \n(2,1)phi - for RZ(phi) \n(3,1)theta - for RY(theta) \n(4,1)lamb - for RZ(lamb)');

output{2,1}(1,1) = gamma;
output{2,1}(2,1) = phi;
output{2,1}(3,1) = theta;
output{2,1}(4,1) = lamb;

%[RZGate1_mat, RYGate_mat, RZGate2_mat];

output{1,2} = sprintf('phase');
output{1,3} = sprintf('RZ(phi)');
output{1,4} = sprintf('RY(theta)');
output{1,5} = sprintf('RZ(lamb)');
output{1,6} = sprintf('coeff to convert to SU');

output{2,2} = phase_gate;
output{2,3} = RZGate_phi;
output{2,4} = RYGate_theta;
output{2,5} = RZGate_lamb;
output{2,6} = coeff;

end




