fprintf("Code Started!\n");
colormap(parula(1024));
M = 128*2;
N = 128*2;
dim = [-1 1 -1 1]*2;
x = linspace(dim(1),dim(2),M+2);
y = linspace(dim(3),dim(4),N+2);
h = x(2)-x(1);

[X,Y] = meshgrid(x,y); X = X'; Y = Y';

rng(28)
Nc = 100;
rad = rand(1,Nc)/7;
cx = 1.5*rand(1,Nc)-0.75;
cy = 1.5*rand(1,Nc)-0.75;

levelset = -sqrt((X+1).^2+(Y).^2)+0.8;
for k = 2:length(cx)
    levelset = max(levelset, -sqrt((X-cx(k)).^2+(Y-cy(k)).^2)+rad(k));
end

Hs = @(x,h) double(x>=h/2)+(x>-h/2&x<h/2).*(0.5+0.5*sin(pi*x/2/(h/2)));
psi = Hs(levelset,h);
rho = 15*psi+(1-psi);
f = rho*0;
f(:,2:end-1) = (rho(:,3:end) - rho(:,1:end-2))/h;
f = f+.01*sin(2*X)+.01*atan(Y); 
f(5,5) = f(5,5)-sum(f(:));
% imagesc(x,y,f')
phi = zeros(size(rho));
iters = 40;

%%%%%%%%%%%%%%%%%%%%%%%% Main Usage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
[phi,  resid] = multigrid_c(phi*0, f, rho, h, iters, 1e-10, 10, 3);
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1)
imagesc(x,y,f')
title('Known: f(x,y)')
xlabel('x-axis');
ylabel('y-axis');
set(gca,'ydir','normal')
colorbar
axis equal
axis(dim)
subplot(2,2,2)
imagesc(x,y,rho')
title('Known: \rho(x,y)')
xlabel('x-axis');
ylabel('y-axis');
set(gca,'ydir','normal')
colorbar
axis equal
axis(dim)
subplot(2,2,3)
imagesc(x,y,phi')
title('Solution: \phi(x,y)')
xlabel('x-axis');
ylabel('y-axis');
set(gca,'ydir','normal')
colorbar
axis equal
axis(dim)
subplot(2,2,4)
semilogy(0:iters,resid,'o-');
title('Residual')
xlabel('Iterations');
ylabel('|L_{\infty}|');xlim([0 sum(resid~=0)-1]);
grid on
grid minor
subplot(2,2,1)
text(0.90,2.5,'Solve: $\nabla\cdot\left(\frac{\nabla \phi(x,y)}{\rho(x,y)} \right )=f(x,y)$','interpreter','latex','fontsize',15)
