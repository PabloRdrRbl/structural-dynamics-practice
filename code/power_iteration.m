%% PART II: Power Iteration Method
%
% Author: Pablo R. Robles (2017)
% Comment: Modification from Prof. Rixen code.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Finite Element assembly of a beam-bar 3D structure
%
%    Project part of the course Structural Dynamics, TUM MW 2136
%    D.R. 05.06.2017
%
%   This file contains the assembly of a general 3D mesh consisting of beam
%   and bar elements. The nodal, material and connectivity matrices are
%   defined in the 'Nodes', 'Mat' and 'Elements' matrices repsectively. The
%   example given here is a hanguar as described in a joined document.
%   The beam element matrices are given. 
%
%%%%%%%%%%%%%%%%%%%%%

 clear

%%%%%%%%%% Geometric parameters %%%%%%%%%%%%%
% 1.387418246946373e+10
% define here all the geometric quantities you need (diameters, lengths, thinkness, etc.)

L=80;
l1=4;
l2=5;
l3=5;
l4=4;
l5=6;
l6=4;
l7=10;
l8=3;
l9=9;
l10=10;

h=25;

H=25;

D1=0.20; d1=0.195;    % beams group 1 (lower arch)
D2=0.20; d2=0.195;    % beams group 2 (mid arch)
D3=0.25; d3=0.244;    % beams group 3 (high arch and span)
D4=0.05; d4=0.046;    % bars in arch
D5=0.40; d5=0.380;    % column beams
D6=0.08; d6=0.075;    % cross bars
D7=0.25; d7=0.24;    % cross beams


%  Material properties 
%
% define here the material constant. 

 E=2.1e11;
 nu=0.3;
 rho=7500;

%%%%%%% Computing characteristics %%%%%%%%%%%%%
%
% Here you compute the quantities needed for the 'Mat' array (see below). 
%

 G = E/(2*(1+nu));

 A1  = pi*(D1^2-d1^2)/4;
 A2  = pi*(D2^2-d2^2)/4;
 A3  = pi*(D3^2-d3^2)/4;
 A4  = pi*(D4^2-d4^2)/4;
 A5  = pi*(D5^2-d5^2)/4;
 A6  = pi*(D6^2-d6^2)/4;
 A7  = pi*(D7^2-d7^2)/4;
 
 I1  = pi*(D1^4-d1^4)/64; J1  = pi*(D1^4-d1^4)/32;
 I2  = pi*(D2^4-d2^4)/64; J2  = pi*(D2^4-d2^4)/32;
 I3  = pi*(D3^4-d3^4)/64; J3  = pi*(D3^4-d3^4)/32;
 I5  = pi*(D5^4-d5^4)/64; J5  = pi*(D5^4-d5^4)/32;
 I7  = pi*(D7^4-d7^4)/64; J7  = pi*(D7^4-d7^4)/32;

%%%%%%%  Nodes    %%%%%%%%
%  x    ,   y   ,   z  
%
% Create here the 'Nodes' matrix. It contains the nodal coordinates x,y,z.
%

Nodes(1,:)  = [ 0               0       H ];    % arch
Nodes(2,:)  = [ 0               0       H+l1];  
Nodes(3,:)  = [ l1              0       H+l1];
Nodes(4,:)  = [ l3              0       H+l1+l2];      
Nodes(5,:)  = [ l3+l1           0       H+l1+l2];
Nodes(6,:)  = [ l3+l5           0       H+l1+l2+l4];
Nodes(7,:)  = [ l3+l5+l1        0       H+l1+l2+l4];
Nodes(8,:)  = [ l3+l5+l7        0       H+l1+l2+l4+l6];
Nodes(9,:)  = [ l3+l5+l7+l9     0       H+l1+l2+l4+l6+l8];   
Nodes(10,:) = [ l3+l5+l7+l9+l10 0       H+l1+l2+l4+l6+l8]; 
for i=21:29,
    Nodes(i,:) = Nodes(i-20,:); 
    Nodes(i,1) = L-Nodes(i,1);
end
Nodes(30,:)   = [0     0    0]; % ground nodes
Nodes(31,:)   = [L     0    0];
for i=1:31,
   Nodes(i+100,:) = Nodes(i,:) + [0 h 0];
end
Nodes(135,:) = [0 h/2 H];
Nodes(136,:) = [L h/2 H];

%%%%%%% Dummy node for beams %%%%%%%%%
%
% this node is use to orient the beam cross section in 3D (see lecture
% notes).
% Since in this example the cross section is axisymmetric any orientation
% is fine. This we choose a dummy node randonly, just making sure it is not
% colinear with the 2 nodes of a beam.

 Nodes(200,:) = [0 40 0]; 

%%%%%%%  Material and section properties    %%%%%
%
%  E    ,   G   , rho ,   A   ,  Ix  ,   Iy   , Iz
%  
% Matrix 'Mat' contains the material and section properties definition. Here below an example with
% 7 different materials. Materials 1, 2, 3, 5 and 7 are for beam elements, materials 4 and 6 are for bar elements.

 Mat(1,1:7)=  [ E   G   rho   A1  J1    I1     I1];  % beams group 1 (lower arch)
 Mat(2,:)  =  [ E   G   rho   A2  J2    I2     I2];  % beams group 2 (mid arch)
 Mat(3,:)  =  [ E   G   rho   A3  J3    I3     I3];  % beams group 3 (high arch and span)
 Mat(4,:)  =  [ E   G   rho   A4  0      0      0];  % bars in arch
 Mat(5,:)  =  [ E   G   rho   A5  J5    I5     I5];  % column beams
 Mat(6,:)  =  [ E   G   rho   A6  0      0      0];  % cross bars
 Mat(7,:)  =  [ E   G   rho   A7  J7    I7     I7];  % cross beams

%%%%%%% Elements   %%%%%%%%%
%
%  type(1=bar,2=beam), mat, node 1 , node 2, node 3(for beam)
%
%  The Matrix 'Elements' contains the connectivity (localization) informations for all the elements of the mesh,
%  together with the material and geometric related properties. 
%  For example, here below: element no.1 is a beam (type=2), is made of
%  material 7, connects node 1 and node 135 and has the dummy node 200
%  (necessary for a beam)

Elements(1,:)  = [2   1       1  2  200];  % arch
Elements(2,:)  = [2   1       1  3  200];
Elements(3,:)  = [2   1       2  3  200];
Elements(4,:)  = [2   2       2  4  200];  
Elements(5,:)  = [2   2       3  5  200]; 
Elements(6,:)  = [2   2       4  6  200]; 
Elements(7,:)  = [2   2       5  7  200]; 
Elements(8,:)  = [2   2       6  8  200];  
Elements(9,:)  = [2   2       7  8  200];
Elements(10,:) = [2   3       8  9  200];
Elements(11,:) = [1   4       3  4  200];
Elements(12,:) = [1   4       4  5  200];
Elements(13,:) = [1   4       5  6  200];
Elements(14,:) = [1   4       6  7  200];
Elements(15,:) = [2   3       9  10 200];   % middle span
Elements(16,:) = [2   3       10 29 200];  
% create right side of arch
for i=1:14,
    Elements(16+i,:) = Elements(i,:);
    Elements(16+i,3:4) = Elements(16+i,3:4) + [20 20];
end
% vertical columns
Elements(31,:) = [2   5       31 21  200];
Elements(32,:) = [2   5       30  1  200];
% duplicate the arch at y=h
for i=1:32,
    Elements(i+32,:) = Elements(i,:)+[0 0 100 100 0];
end
% add beams for crane rails
Elements(65,:) = [2   7       1   135 200];
Elements(66,:) = [2   7       135 101 200];
Elements(67,:) = [2   7       21  136 200];
Elements(68,:) = [2   7       136 121 200];
% add cross bar and beams
Elements(69,:) = [1   6       9  129  200];
Elements(70,:) = [1   6       29 109  200];
Elements(71,:) = [1   6       8  108  200];
Elements(72,:) = [1   6       28 128  200];
Elements(73,:) = [1   6       6  104  200];
Elements(74,:) = [1   6       26 124  200];
Elements(75,:) = [1   6       4  106  200];
Elements(76,:) = [1   6       24 126  200];

% plot the mesh

figure(1),plotmesh(Nodes,Elements)  

%%%%%%%%% build elementary matrices and assemble  %%%%%%

% build index table for dof: locnod(i,:) = list of degrees of freedom number attached to node i
%
Nodes_active = unique(Elements(:,3:4));   % find which node really are used in elements
                                          % the degrees of freedom fixed on the ground will
                                          % be removed later.
Ndof = size(Nodes_active,1)*6;            % 6 dof per node
locnod(Nodes_active,1:6) = reshape([1:Ndof],6,Ndof/6)';         
Nele   = size(Elements,1);

K = sparse(Ndof,Ndof);
M = sparse(Ndof,Ndof);

for el = 1:Nele,

   type   = Elements(el,1);
   matel  = Mat(Elements(el,2),:);
   node1  = Nodes(Elements(el,3),:);  % coordinates of node 1
   node2  = Nodes(Elements(el,4),:);  % coordinates of node 2
   node3  = Nodes(Elements(el,5),:);  % coordinates of node 3 (for beams)

   % Build element matrices
   if type==1,
      [Kel,Mel] = bar(node1,node2,matel);
         dof=[locnod(Elements(el,3),[1 2 3]) ...
              locnod(Elements(el,4),[1 2 3]) ] ; % only the translational
   elseif type==2,
      [Kel,Mel] = beam(node1,node2,node3,matel);  
         dof=[locnod(Elements(el,3),:) ...
              locnod(Elements(el,4),:) ] ; 
   end
   
   %% Assemble
   K(dof,dof) = K(dof,dof)+Kel;
   M(dof,dof) = M(dof,dof)+Mel;
      
   clear matel type Kel Mel dof
end

%% Apply boundary conditions
% fix dofs of nodes 30 31 130 131
% a. find the degrees of freedom that are fixed
dof_fix = [];
for n=[30 31 130 131], 
    dof_fix =   [dof_fix  locnod(n,:) ];
end
Ndof_fix = size(dof_fix,2);
% b. remove these fixed dof from the list of dof and eliminate them in the matrices
dof_rem = setdiff([1:Ndof],dof_fix);% remaining degrees of freedom
Ndof_rem=Ndof-Ndof_fix;
K=K(dof_rem,dof_rem);   % from here on, K and M are related to the dof_remaining
M=M(dof_rem,dof_rem);

%% Compute eigensolutions: Power iteration method

% Make K and M full
K = full(K);
M = full(M);

% Check K is not singular using rank
[m, n] = size(K);
if rank(K) < m 
    error('K is not singular')
end

% Dynamic flexibility matrix or iteration matrix
D = inv(K) * M;

% Tolerance for convergence check using Rayleight quotients
eps = 1e-5;

% Arbitrary starting vector
z0 = rand(n, 1);
z = z0;

% Rayleight quotient
rayquo = (z' * K * z) ./ (z' * M * z);

% Previous Rayleight quotient
rayquo_0 = 0;

% Convergence condition
condition = 1;

iter = 0;

while condition
    
    % New iterate z
    z = D * z;
    
    % Normalized iterate
    z = z / norm(z);
    
    % Rayleight quotient
    rayquo_0 = rayquo;
    rayquo = (z' * K * z) / (z' * M * z);
    
    % Convergence condition. If greater than the tolerance keep iterating
    condition = abs(rayquo - rayquo_0) > eps;
    
    iter = iter + 1;
end

Omega2 = rayquo;

%% Compute eigensolutions: Matlab eig

[X,Omega2] = eig(full(K),full(M));
[Omega2,k] = sort(real(diag(Omega2)));
X = real(X(:,k));
clear k
omega2m = Omega2(1:neig);
sqrt(Omega2)/(2*pi)

figure(2)
Xplot=zeros(Ndof,Ndof_rem); % express the modes in the non-fixed numbering used in locnod
Xplot(dof_rem,:)=X;
plotmesh(Nodes,Elements);
plotdis(Nodes,Elements,locnod,Xplot(:,1),5);