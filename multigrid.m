function [phi, err, resid] = multigrid(phi,f, rho, h, nIterMax, alpha, smoothCounts, cycleType)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MULTIGRID: Calculates the variable coefficient poisson equation using V,W 
% or F cycle in 2D with nuemnann boundary condition of d(phi(x,y))/d(n) 
% equal to zero for M by M mesh points for equal dx and dy mesh size                                                                                                              
% PDE being solved:
% 
%                                      yMMNNNNNNNNMd`     o                                     
%                                .o:    sMm.     :y`   - .oo:    :o-                      `/+so`
% -oooooooooooooooooo-         `yMd:     +MN-   /s    ym ms:N+   -dMh`                   sMmyo+`
% `hMMMNmmmmmmmmmmNMh`        `mMo        :NM/`o+     M+ M+ sN     /Mm.                 sMy`    
%   yMMm.        `hs          dMo          -mMh/      yh/Ny/ds      +Mm`  `........ `oyyMMhyyo  
%    oMMN:      .do    `oo`  /MN            .o-        -/Nh/-        mM+  oddddddds `++mMy++/:  
%     /NMN+    -d/     .hh.  yMy                         h:          sMd               NM.      
%      :NMMo  /d:            hMs      oooooooooooooooooooooooooo     oMd  .-------.   :Md       
%       .mMMyod.             oMd                .:/:`                yMy  +hhhhhhh/   yM+       
%        `hMMh`              -MM.              sms+ym-              `NM:             .MN`       
%         `//`                oMm.            :M/   Mo             .hMs              -y/        
%                              +NNo.          +My::sm.           .oNNo                          
%                               `/y+          oM/oo/`            /y+`                           
%                                             yd                                                
%
%   • phi: 2D array of cell centered solution variable including ghost cells, containing 
%   upon call the initial guess of the solution. Size = (M+2) x (N+2)
%
%   • f: 2D array of cell centered PDE right hand side including ghost cells (even 
%   though these are not necessarily defined). Size = (M+2) x (N+2)
%
%   • rho: 2D array of cell centered PDE variable coefficents including ghost cells 
%   (even though these are not necessarily defined). Size = (M+2) x (N+2)
%
%   • h: mesh spacing
%
%   • nIterMax: integer number of maximum V-cycle iterations to be performed
%
%   • alpha: convergence threshold for residual
%
%   • smoothCounts: Number of smoothing operations in pre and post smoothing
%
%   • cycleType: 1 for V-Cycle, 2 for F-Cycle and 3 for W-Cycle
%
% OutPuts:
%   • phi: 2D array of cell centered solution variable including ghost cells, containing 
%   the solution
%
%   • err: error code. Set to 0 if solution converged and no errors occurred, set to 1 
%   if the solution did not converge.
%
%   • resid: array of residuals for each iteration
%
% Note:
%   • Ensure prime factors of M and N is mostly comprised of 2. Any non-two 
%   factors will slow down convergence. 
%
% Author: Adil Ansari, MS
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global counts;
counts = smoothCounts;
resid = zeros(1,nIterMax+1); % define residual matrix
resid(1) = max(max(abs(residual(phi,f,rho,h)))); % Compute initial resid
err = 1;
cycles = {@VCycle,@FCycle,@WCycle}; % store cycle functions in cell array
for k = 2:nIterMax+1
    phi = cycles{cycleType}(phi,f,rho,h);
    resid(k) = max(max(abs(residual(phi,f,rho,h))));
    if resid(k)<alpha
        err = 0;
        break
    end
end
end

function phi = FCycle(phi,f,rho,h)
    % Recursive F-Cycle Multigrid
    global counts;
    % Pre-Smoothing
	v = smooth(phi,f,rho,h,counts);
	r = residual(v,f,rho,h);
	rhs = restrict(r);
    rho_2h = restrict(rho);
    
    M = size(rhs,1)-2;    N = size(rhs,2)-2;
    disp(log2(M))
    eps = zeros(M+2,N+2);
    % stop recursion when encountering odd sized grid 
    if mod(M,2) && mod(N,2) 
        eps = smooth(eps,rhs,rho_2h,2*h,M*N);
    else
        % "We need to go deeper!" - Leonardo Dicaprio
        eps = FCycle(eps,rhs,rho_2h,2*h);        
    end
	v = v + prolong(eps);
    % fracture iteration into a w
    v = smooth(v,f,rho,h,counts);
	r = residual(v,f,rho,h);
	rhs = restrict(r);
    rho_2h = restrict(rho);
    
    M = size(rhs,1)-2;    N = size(rhs,2)-2;

    % stop recursion when encountering odd sized grid 
    eps = zeros(M+2,N+2);
    if mod(M,2) && mod(N,2) 
        eps = smooth(eps,rhs,rho_2h,2*h,M*N);
    else
        % "We need to go deeper again!" - Leonardo Dicaprio
        eps = VCycle(eps,rhs,rho_2h,2*h);        
    end
	v = v + prolong(eps);
    % Post-Smoothing
	phi = smooth(v,f,rho,h,counts);
    disp(log2(M))
end

function phi = VCycle(phi,f,rho,h)
    % Recursive V-Cycle Multigrid
    global counts;
    % Pre-Smoothing
	v = smooth(phi,f,rho,h,counts);
	r = residual(v,f,rho,h);
	rhs = restrict(r);
    rho_2h = restrict(rho);
    
    M = size(rhs,1)-2;    N = size(rhs,2)-2;
    eps = zeros(M+2,N+2);
    % stop recursion when encountering odd sized grid 
    if mod(M,2) && mod(N,2) 
        eps = smooth(eps,rhs,rho_2h,2*h,M*N*counts);
    else
        % "We need to go deeper!" - Leonardo Dicaprio
        eps = VCycle(eps,rhs,rho_2h,2*h);        
    end
	v = v + prolong(eps);
    % Post-Smoothing
	phi = smooth(v,f,rho,h,counts);    
end

function phi = WCycle(phi,f,rho,h)
    % Recursive W-Cycle Multigrid
    global counts;
    % Pre-Smoothing
	phi = smooth(phi,f,rho,h,counts);
	r = residual(phi,f,rho,h);
	rhs = restrict(r);
    rho_2h = restrict(rho);
    
    M = size(rhs,1)-2;    N = size(rhs,2)-2;
    eps = zeros(M+2,N+2);
    % stop recursion when encountering odd sized grid 
    if mod(M,2) && mod(N,2) 
        eps = smooth(eps,rhs,rho_2h,2*h,M*N);
    else
        % "We need to go deeper!" - Leonardo Dicaprio
        eps = WCycle(eps,rhs,rho_2h,2*h);        
    end
	phi = phi + prolong(eps);
    % fracture iteration into 2 parts
    phi = smooth(phi,f,rho,h,counts);
	r = residual(phi,f,rho,h);
	rhs = restrict(r);
    rho_2h = restrict(rho);
    
    M = size(rhs,1)-2;    N = size(rhs,2)-2;
    eps = zeros(M+2,N+2);
    % stop recursion when encountering odd sized grid 
    if mod(M,2) && mod(N,2) 
        eps = zeros(M+2,N+2);
        eps = smooth(eps,rhs,rho_2h,2*h,M*N);
    else
        % "We need to go deeper again!" - Leonardo Dicaprio
        eps = WCycle(zeros(M+2,N+2),rhs,rho_2h,2*h);        
    end
	phi = phi + prolong(eps);
    
    % Post-Smoothing
	phi = smooth(phi,f,rho,h,counts);
end

function phi = smooth(phi,f,rho,h,n)
% Smoothing operator

M = size(phi,1)-2;
N = size(phi,2)-2;

% Start Gauss Seidel Iterations
for k = 1:n
    for i = (2:M+1)
        for j = (2:N+1)
            im = 1/(h*h)*(...
                1./(rho(i+1,j)+rho(i,j))+...
                1./(rho(i,j)+rho(i-1,j))+...
                1./(rho(i,j+1)+rho(i,j))+...
                1./(rho(i,j-1)+rho(i,j)));
            pg = 1/(h*h)*(...
                phi(i+1,j)./(rho(i+1,j)+rho(i,j))+...
                phi(i-1,j)./(rho(i,j)+rho(i-1,j))+...
                phi(i,j+1)./(rho(i,j+1)+rho(i,j))+...
                phi(i,j-1)./(rho(i,j-1)+rho(i,j)))-f(i,j);
            phi(i,j) = 1.0./im.*pg;
        end
    end
    
    % update ghost points
    phi = BC(M,N,phi);
end
end

function resid = residual(phi,f,rho,h)
%calculates the residuals
[M,N] = size(phi);
M = M - 2; N = N - 2;

resid = zeros(M+2,N+2);
i = 2:M+1;
j = 2:N+1;
pg = 1/(h*h)*(...
    (phi(i+1,j)-phi(i,j))./(rho(i+1,j)+rho(i,j))-...
    (phi(i,j)-phi(i-1,j))./(rho(i,j)+rho(i-1,j))+...
    (phi(i,j+1)-phi(i,j))./(rho(i,j+1)+rho(i,j))-...
    (phi(i,j)-phi(i,j-1))./(rho(i,j-1)+rho(i,j)));
resid(i,j) = f(i,j) - pg;

% update ghost points
resid = BC(M,N,resid);
end

function A = restrict(B)
% make mesh coarser
M = (size(B,1)-2)/2;
N = (size(B,2)-2)/2;
i = 2:M+1;
j = 2:N+1;
A = zeros(M+2,N+2);
A(i,j) = 0.25*(B(2*i-1,2*j-1) + ...
    B(2*i-2,2*j-1) + ...
    B(2*i-1,2*j-2) + ...
    B(2*i-2,2*j-2));

% update ghost points
A = BC(M,N,A);
end

function B = prolong(A)
% create finer mesh with bilinear interpolation
M = size(A,1)-2;
N = size(A,2)-2;

B = zeros(2*M+2,2*N+2);
i = 2:M+1;
j = 2:N+1;
B(2*i-1, 2*j-1) = 0.5625*A(i,j) + 0.1875*(A(i,j+1)+A(i+1,j)) + 0.0625*A(i+1,j+1);
B(2*i-2, 2*j-1) = 0.5625*A(i,j) + 0.1875*(A(i,j-1)+A(i+1,j)) + 0.0625*A(i-1,j+1);
B(2*i-1, 2*j-2) = 0.5625*A(i,j) + 0.1875*(A(i,j+1)+A(i-1,j)) + 0.0625*A(i+1,j-1);
B(2*i-2, 2*j-2) = 0.5625*A(i,j) + 0.1875*(A(i,j-1)+A(i-1,j)) + 0.0625*A(i-1,j-1);

% update ghost points
B = BC(2*M,2*N,B);
end

function A = BC(M,N,A)
A(1,:) =  A(2,:);
A(:,1) =  A(:,2);
A(M+2,:) =  A(M+1,:);
A(:,N+2) =  A(:,N+1);
end