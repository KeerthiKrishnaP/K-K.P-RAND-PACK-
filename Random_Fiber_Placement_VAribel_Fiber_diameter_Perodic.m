%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Created by: Keerthi Krishna PARVATHANENI
%% Date:::18-09-2019
%% NO Licence Needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc
r=0.014/2;
UL=0.20;LL=0.0; Rf=0.019;P=0;
bx=[LL UL UL LL];
by=[LL LL UL UL];
AP=((UL+0.007)*(UL+0.007));
B = convhull(bx,by);
polyarea(bx(B),by(B));
plot(bx,by,'.')
axis equal
hold on
fill ( bx(B),by(B), 'K','facealpha', 0.5 ); 
%%Main_Code
NumTrail=0;

while P<0.40
    NumTrail=NumTrail+20;
L3x=LL;U3x=UL;
M=[1 1];
B3Cx=(L3x-U3x).*rand(NumTrail,1)+U3x;
L3y=LL;U3y=UL;
B3Cy=(L3y-U3y).*rand(NumTrail,1)+U3y;
%%%Variation of fiber radius
Rmin=0.010; Rmax=0.014;
R=(Rmin-Rmax)*rand(NumTrail,1)+Rmax;
%%%%%%%%%ù
Data3=[B3Cx, B3Cy,R];
%%%%%%%ù
AL=Data3(1,1:3);
AL(1,4)=(pi*Data3(1,3)^2)/4;
for i=2:NumTrail
    tempx=Data3(i,1);
    tempy=Data3(i,2);
    tempR=Data3(i,3);
    tempDx=zeros(1);
    Ndx=zeros(1);
    Edge_element=0;
for j=1:M(1)
    tempDx(j)=sqrt((tempx-AL(j,1))^2+(tempy-AL(j,2))^2);
end
Pass=0;
    for k=1:length(tempDx)
        T2=0;T3=0;T4=0;
        T_1=0;
    if tempDx(k)>(R(i)/2+AL(k,3)/2+0.2*(R(i)/2+AL(k,3)/2))
        Pass=Pass+1;
    end
    if Pass==length(tempDx)
        d4=tempx;
        d1=tempy;
        d2=UL-tempx;
        d3=UL-tempy;
     %corber element check
     if (d1<tempR && d4<tempR)
     Edge_element=1;
         N(1,:)=[tempx, tempy,tempR];
         N(2,:)=[UL+d4,tempy,tempR];
         N(3,:)=[UL+d4,UL+d1,tempR];
         N(4,:)=[tempx,UL+d1,tempR];
         for n=2:4
             T2=0;T3=0;T4=0;
            for q=1:M(1)
              Ndx=sqrt((N(n,1)-AL(q,1))^2+(N(n,2)-AL(q,2))^2);
              if Ndx>(R(i)/2+tempR/2+0.2*(R(i)/2+tempR/2))
                  if n==2
                      T2=T2+1;
                  elseif n==3
                      T3=T3+1;
                  elseif n==4
                      T4=T4+1;
                  end
              end
            end
         end
     %
     elseif (d1<tempR && d2<tempR)
         Edge_element=1;
         N(1,:)=[tempx, tempy,tempR];
         N(2,:)=[UL+d2,UL+d1,tempR];
         N(3,:)=[-d2,UL+d1,tempR];
         N(4,:)=[-d2,tempy,tempR];
         for n=2:4
             T2=0;T3=0;T4=0;
            for q=1:M(1)
              Ndx=sqrt((N(n,1)-AL(q,1))^2+(N(n,2)-AL(q,2))^2);
              if Ndx>(R(i)/2+tempR/2+0.2*(R(i)/2+tempR/2))
                  if n==2
                      T2=T2+1;
                  elseif n==3
                      T3=T3+1;
                  elseif n==4
                      T4=T4+1;
                  end
              end
            end
         end
     elseif (d2<tempR && d3<tempR)
         Edge_element=1;
         N(1,:)=[tempx, tempy,tempR];
         N(2,:)=[-d2,tempy,tempR];
         N(3,:)=[-d2,-d3,tempR];
         N(4,:)=[tempx,-d3,tempR];
         for n=2:4
             T2=0;T3=0;T4=0;
            for q=1:M(1)
              Ndx=sqrt((N(n,1)-AL(q,1))^2+(N(n,2)-AL(q,2))^2);
              if Ndx>(R(i)/2+tempR/2+0.2*(R(i)/2+tempR/2))
                  if n==2
                      T2=T2+1;
                  elseif n==3
                      T3=T3+1;
                  elseif n==4
                      T4=T4+1;
                  end
              end
            end
         end
     elseif (d3<tempR && d4<tempR)
         Edge_element=1;
         N(1,:)=[tempx, tempy,tempR];
         N(2,:)=[tempx,-d3,tempR];
         N(3,:)=[UL+d4,-d3,tempR];
         N(4,:)=[UL+d4,tempy,tempR];
         for n=2:4
             T2=0;T3=0;T4=0;
            for q=1:M(1)
              Ndx=sqrt((N(n,1)-AL(q,1))^2+(N(n,2)-AL(q,2))^2);
              if Ndx>(R(i)/2+tempR/2+0.2*(R(i)/2+tempR/2))
                  if n==2
                      T2=T2+1;
                  elseif n==3
                      T3=T3+1;
                  elseif n==4
                      T4=T4+1;
                  end
              end
            end
         end
     elseif d1<tempR
          T_1=0;
         Edge_element=1;
         E(1,:)=[tempx, tempy, tempR];
         E(2,:)=[tempx,UL+d1,tempR];
         for q=1:M(1)
              Edx=sqrt((E(2,1)-AL(q,1))^2+(E(2,2)-AL(q,2))^2);
              if Edx>(R(i)/2+tempR/2+0.2*(R(i)/2+tempR/2))
                  T_1=T_1+1;
              end
         end
     elseif d2<tempR
          T_1=0;
         Edge_element=1;
         E(1,:)=[tempx, tempy, tempR];
         E(2,:)=[-d2,tempy,tempR];
         for q=1:M(1)
              Edx=sqrt((E(2,1)-AL(q,1))^2+(E(2,2)-AL(q,2))^2);
              if Edx>(R(i)/2+tempR/2+0.2*(R(i)/2+tempR/2))
                  T_1=T_1+1;
              end
         end
     elseif d3<tempR
          T_1=0;
         Edge_element=1;
         E(1,:)=[tempx, tempy, tempR];
         E(2,:)=[tempx,-d3,tempR];
         for q=1:M(1)
              Edx=sqrt((E(2,1)-AL(q,1))^2+(E(2,2)-AL(q,2))^2);
              if Edx>(R(i)/2+tempR/2+0.2*(R(i)/2+tempR/2))
                  T_1=T_1+1;
              end
         end         
     elseif d4<tempR
          T_1=0;
         Edge_element=1;
         E(1,:)=[tempx, tempy, tempR];
         E(2,:)=[UL+d4,tempy,tempR];
         for q=1:M(1)
              Edx=sqrt((E(2,1)-AL(q,1))^2+(E(2,2)-AL(q,2))^2);
              if Edx>(R(i)/2+tempR/2+0.2*(R(i)/2+tempR/2))
                  T_1=T_1+1;
              end
         end         
     end
     if Edge_element==1
     if (T2==length(AL) && T3==length(AL) && T4==length(AL))
         for m=1:4
           AL(end+1,1:3)=N(m,:);
           AL(end,4)=(pi*tempR^2)/4;
         end
     elseif (T_1==length(AL))
           AL(end+1,1:3)=E(1,:);
           AL(end,4)=(pi*tempR^2)/4;
           AL(end+1,1:3)=E(2,:);
           AL(end,4)=(pi*tempR^2)/4;    
     else
         continue
     end
     else
    AL(end+1,1:3)=[tempx,tempy,tempR];
    AL(end,4)=(pi*tempR^2)/4;
    end
    else 
        continue
    end
    end
M=size(AL); 
Af=sum(AL(:,4));
end
P=Af/AP
end
  
%norm = normpdf(AL(:,3),pd[1],pd[2]);
for k=1:M(1)
     r=AL(k,3)/2;
     x=AL(k,1);
     y=AL(k,2);
 th = 0:pi/50:2*pi;
 xunit = r * cos(th) + x;
 yunit = r * sin(th) + y;   
 h = plot(xunit, yunit,'w');
end
hold off
% AL(:,3)=AL(:,3)/2;
%%Write the information of fiber centers
fileID = fopen('Center of fibers Varible FD.txt','wt');
fprintf(fileID,'%6s,%12s,%6s\n','Xc','Yc','R');
fprintf(fileID,'\n');
for ii=1:M(1)
fprintf(fileID,'%12.8f,%12.8f,%12.8f\n',AL(ii,1:3));
end
fclose(fileID);
% pd=fitdist(AL(2:end,3),'Normal');
% m = mean(pd);
% y = pdf(pd,AL(2:end,3));
% hold on
% histogram(AL(2:end,3),'Normalization','pdf')
% %scatter(AL(2:end,3),y)
% hold off