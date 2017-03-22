%CompareBoundaries models and plots diffusion processes generated with
% MMMilstein, and applies a range of boundary conditions to these processes
%
% Lachlan Bridges
% 14/01/17

seed=10113;
N=100;
r = 3;  %number of realisations to plot
reps=10;

T=[0,30];
Q=[-3,3;1,-1];
x0=0;
mu = [5,-3];
sigma = [1,10];
drift = @(y,t,J) mu(J);
diff = @(y,t,J) sigma(J);
dbdy = @(y,t,J) 0;

%boundary conditions
lb=-20;
ub=[];
ljump = 10;
ujump = -10;

%cell of boundary condition functions and there names
boundtype=cell(2);
boundtype{1,1}='absorbing';
boundtype{1,2}='instantaneousjump';
boundtype{1,3}='reflected';
%boundtype{1,4}='reflected';
%boundtype{1,5}='sticky';
boundtype{2,1}=@(x,lb,ljump,ub,ujump) BMAbsorbing(x,lb,ub);
boundtype{2,2}=@(x,lb,ljump,ub,ujump) BMInstantJump(x,lb,ljump,ub,ujump);
boundtype{2,3}=@(x,lb,ljump,ub,ujump) BMReflect(x,lb);
%boundtype{2,4}=@(x,lb,ljump,ub,ujump) BM2Reflect(x,lb,ub);
%boundtype{2,5}=@(x,lb,ljump,ub,ujump) BMSticky(x,lb);

%generating original diffusion process
[t, x, ~, ~, ~] = MMMilstein(Q,reps,drift,diff,dbdy,x0,T,N,seed);

%used in titles and saving of plot
bounds=[lb;ub];
if numel(bounds)==1
    str = sprintf('%i',bounds(1));
elseif numel(bounds)==2
    str = sprintf('%i and %i',bounds(1),bounds(2));
end

%iterates through each boundary condition
for i=1:size(boundtype,2)
    figure;
    y = boundtype{2,i}(x,lb,ljump,ub,ujump); % applying boundary condition
    plot(t(1,:),x(1,:),'b-'), hold on
    plot(t(1,:),y(1,:),'r-')
    xlabel('t','FontSize',10)
    ylabel('X_t','FontSize',10)
    title(sprintf('No boundary condition (blue)\n vs\n %s boundary condition (red)',boundtype{1,i}),...
        'FontSize', 8)
    %saveas(gcf,sprintf('../figures/bound_%s_mu_%1.0f_%1.0f_sig_%1.0f_%1.0f_lbound_%i_ubound_%i.png'...
    %    ,boundtype{1,i},mu(1),mu(2),sigma(1),sigma(2),lb,ub));
    %hold off
    %close all
    % plotting, labelling and saving figures
    figure;
    ax1 = subplot(2,2,1);
    plotmm(t(1:r,:),x(1:r,:))
    title(sprintf('No boundary conditions\n (%i realisations)',r),...
        'FontSize', 8)
    ylabel('X_t','FontSize',10)
    ax2 = subplot(2,2,2);
    [tm,xm] = meanmm(t,x,10000);
    plot(tm,xm,'b-')
    title(sprintf('No boundary conditions\n (average of %i realisations)'...
        ,N),'FontSize', 8)
    ax3 = subplot(2,2,3);
    plotmm(t(1:r,:),y(1:r,:))
    title(sprintf('%s boundary at %s\n (%i realisations)',...
        boundtype{1,i}, str, r),'FontSize', 8)
    ylabel('Y_t','FontSize',10)
    xlabel('t','FontSize',10)
    ax4 = subplot(2,2,4);
    [tm,xm] = meanmm(t,y,10000);
    plot(tm,xm,'r-')
    title(sprintf('%s boundary at %s\n (average of %i realisations)',...
        boundtype{1,i},str,N),'FontSize', 8)
    xlabel('t','FontSize',10)
    linkaxes([ax2,ax3,ax1, ax4],'xy');
    %saveas(gcf,sprintf('../figures/bound_%s_comp_mu_%1.0f_%1.0f_sig_%1.0f_%1.0f_lbound_%i_ubound_%i.png'...
    %    ,boundtype{1,i},mu(1),mu(2),sigma(1),sigma(2),lb,ub));
    %close all
end