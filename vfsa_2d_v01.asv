%%%%%%%%%%%%%%%% 
close all
clear 
clc
funcy=3; %%%% 1 for sign function 2D with global minimum at(0,0) where f(x, y)=0.0
%%%%%%%%%%%%% 2 for Banana(Rosen)function with global minimum at(1,1) where f(x, y)=0.0
%%%%%%%%%%%%% 3 for Griewank function with global minimum at(0,0) where f(x, y)=0.0
reproduce=0; % 0 for random, 1 for reproducable results
nruns=1; %%%%% Number of runs(cases) 
n=100; % number of iteration per run
xx=zeros(n,nruns); %%% tracer for x position at different iteration 
yy=zeros(n,nruns); %%% tracer for y position at different iteration 
err_plot=zeros(n,nruns); %%% tracer for error at different iteration per run
xx_best=zeros(nruns,1); yy_best=zeros(nruns,1); 
error_tresh=0.001;
%%%% range of two varaibles for optimization
xmax=10; xmin=-10; dx=0.05;
ymax=10; ymin=-10; dy=0.05;
nx=round(abs(xmax-xmin)./dx);
ny=round(abs(ymax-ymin)./dy);


%%%% Parameter definition for basic VFSA
    nmov=3;
    temp0=1;
    decay=0.999;
    t0(1)=1;
    t0(2)=1;   


if reproduce==1;
    rand('seed',31415927)
    randn('seed',3111113)
end


%%%% Main loop
for jrun=1:nruns 
        xmod=monte_carlo_sample(xmin,dx,nx,xmax);%%% Initial guess
        ymod=monte_carlo_sample(ymin,dy,ny,ymax);
% %         emod=(1-func1(xmod,ymod)).^2;
        emod=(1-funcy_gen2d(xmod,ymod,funcy)).^2;
        jtemp=1;
        xx(jtemp,jrun)=xmod;
        yy(jtemp,jrun)=ymod;  
        err_plot(jtemp,jrun)=funcy_gen2d(xmod,ymod,funcy);
        while(jtemp<=n-1)
            temp(jtemp)=temp0.*exp(-decay.*(jtemp-1).^0.5);
            tmp1=t0(1).*exp(-decay.*(jtemp-1).^0.5); tmp2=t0(2).*exp(-decay.*(jtemp-1).^0.5);
            for jmov=1:nmov
                xtrial=walk(xmod,dx,xmin,xmax,tmp1);
                ytrial=walk(ymod,dy,ymin,ymax,tmp2);    
                xtrial=round((xtrial-xmin)./dx).*dx+xmin;
                ytrial=round((ytrial-ymin)./dy).*dy+ymin;
                etrial=funcy_gen2d(xtrial,ytrial,funcy);
                if etrial< emod
                    %%%% hist_updat
                    emod=etrial;
                    xmod=xtrial;
                    ymod=ytrial;
                else
                    arg=(etrial-emod)./temp(jtemp);
                    if arg>1.e6
                        pde=0.001;
                    else
                        pde=exp(-arg);
                    end
                    if pde>rand
                        %%%% hist_updat
                        emod=etrial;
                        xmod=xtrial;
                        ymod=ytrial;
                    end
                end
            end %%%% end move
            err_plot(jtemp+1,jrun)=emod;    
            xx(jtemp+1,jrun)=xmod;
            yy(jtemp+1,jrun)=ymod;
            %%%%%%%%%%%%% Exit from temp loop if emod is so small
            if emod<=error_tresh
                jtemp=n;                
            end            
            jtemp=jtemp+1;       
        end
        xx_best(jrun)=xmod;
        yy_best(jrun)=ymod;      
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plotting
lw=2;
fs=12;
%%%% range of two varaibles for plotting only 
xmax=10; xmin=-10; dx=0.05;
ymax=10; ymin=-10; dy=0.05;
nx=round(abs(xmax-xmin)./dx);
ny=round(abs(ymax-ymin)./dy);
x=xmin:dx:xmax;
y=ymin:dy:ymax;
%%%%% generate function to display
myfun=zeros(ny,nx);
for i=1:nx
    for j=1:ny
        myfun(j,i)=funcy_gen2d(x(i),y(j),funcy);
    end 
end

figure(1)
imagesc(x,y,myfun);
colorbar;xlabel('X Axes','FontSize',fs);ylabel('Y Axes','FontSize',fs);set(gca,'FontSize',fs)
title('Objective Function','FontSize',fs);
for jrun=1:nruns
    hold on
    plot(xx(:,jrun),yy(:,jrun),'--w','LineWidth',lw); %%% all iteration
    hold on
    plot(xx(1,jrun),yy(1,jrun),'dw'); %%% first iteration for jrun(initial guess)
    hold on
    plot(xx(n,jrun),yy(n,jrun),'*k');%%% last iteration for jrun(solution)
    colorbar;xlabel('X Axes','FontSize',fs);ylabel('Y Axes','FontSize',fs);set(gca,'FontSize',fs)
end


figure(2)
plot(err_plot,'r','LineWidth',lw)
xlabel('Iteration number','FontSize',fs);
ylabel('Objective function','FontSize',fs);
set(gca,'FontSize',fs)

figure(3)        
    subplot(1,2,1)
    plot(xx,'b','LineWidth',lw)
    xlabel('Iteration number','FontSize',fs);
    ylabel('First model parameter(x)','FontSize',fs);
    ylim([xmin,xmax]);
    grid;
    set(gca,'FontSize',fs)
    subplot(1,2,2)
    plot(yy,'b','LineWidth',lw)
    xlabel('Iteration number','FontSize',fs);
    ylabel('Second model parameter(y)','FontSize',fs);
    ylim([ymin,ymax]);
    grid;
    set(gca,'FontSize',fs)

if nruns>1
    figure(4)
    subplot(3,1,1)
    imagesc(err_plot');colorbar;
    ylabel('Number of runs','FontSize',fs);
    title('Objective Function','FontSize',fs);
    subplot(3,1,2)
    imagesc(xx');colorbar;
    ylabel('Number of runs','FontSize',fs);
    title('First model parameter(x)','FontSize',fs);
    subplot(3,1,3)
    imagesc(yy');colorbar;
    xlabel('Iteration number','FontSize',fs);
    ylabel('Number of runs','FontSize',fs);
    title('Second model parameter(y)','FontSize',fs);
end




        
    
    
    
    
    
    
    
    