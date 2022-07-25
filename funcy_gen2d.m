function func=funcy_gen2d(x,y,funcy)

if funcy==1
    func1=sign(sin(x)./x).*abs(sin(x)./x).^0.25.*sign(sin(y)./y).*abs(sin(y)./y).^0.25;
    func=(1-func1).^2;
elseif funcy==2
    func=(1-x).^2 + 100.*(y-x.^2).^2;
elseif funcy==3
    func=1+(x.^2+y.^2)./4000-cos(x).*cos(y./2);    
else
    sprintf('Other functions not included yet')
    return

end


% % close all
% % clear
% % clc
% % %%%% 2D plotting
% % funcy=3;
% % xmax=10; xmin=-10; dx=0.05;
% % ymax=10; ymin=-10; dy=0.05;
% % nx=round(abs(xmax-xmin)./dx);
% % ny=round(abs(ymax-ymin)./dy);
% % x=xmin:dx:xmax;
% % y=ymin:dy:ymax;
% % %%%%% generate function to display
% % myfun=zeros(ny,nx);
% % for i=1:nx;
% %     for j=1:ny;
% %         myfun(j,i)=funcy_gen2d(x(i),y(j),funcy); %%%% Sign function in VFSA book        
% %     end 
% % end
% % 
% % imagesc(x,y,myfun);


