clc
clear all
close all

%load the constants 
[~,~,data]=xlsread('mechanical spreadsheet.xlsx','Constants'); 
global p Cd         %global
drod=data{3,3};     %diameter of the rod
dpiston=data{4,3};  %Effective diameter of the piston

%dependent variables
Arod=pi*drod^2/4;        %Rod area
Apiston=pi*dpiston^2/4;     %Piston area
Aannular=Apiston-Arod; %Annular area
L=dpiston/4;            %shim length

%load the inputs
[~,~,data]=xlsread('mechanical spreadsheet.xlsx','Inputs'); 

p=data{22,3};       %density of lubricant
Cd= data{23,3};     %coeficient of drag
E=data{24,3};       %young Modulus


V=cell2mat(data(3:32,8));         %velocities
Force=cell2mat(data(3:32,9));    %Actual force


dfvalveU=data{3,5};  %the upper limit of diameter of the fixed valve
dVvalveU=data{4,5};  %the upper limit of diameter of the variable valve
tshimU=data{5,5};  %the upper limit of  shim thickness
preloadU=data{6,5};  %the upper limit of diameter of the shim thickness
bshimU=data{7,5};     %the upper limit of the width of shim 

dfvalveL=data{3,4};  %the lower limit of diameter of the fixed valve
dVvalveL=data{4,4};  %the lower limit of diameter of the variable valve
tshimL=data{5,4};    %the  lower limit of the shim thickness
preloadL=data{6,4};  %the upper limit of diameter of the shim thickness
bshimL=data{7,4};     %the  lower limit of the width of shim 

dfvalves=linspace(dfvalveL,dfvalveU,data{3,6});    %valve diameter
dVvalves=linspace(dVvalveL,dVvalveU,data{4,6});    %valve diameter
tshim=linspace(tshimL,tshimU,data{5,6});          %shim thickness
preloads=linspace(preloadU,preloadL,data{6,6});    %valve diameter
bshim=linspace(bshimU,bshimL,data{7,6});    %width of shim


results=[]; %result

for i=1:length(dfvalves)
    for j=1:length(dVvalves)
        for k=1:length(bshim)
              for n=1:length(tshim)
                        
                dfvalve=dfvalves(i);    %diameter of the constantly opened valve
                dVvalve=dVvalves(j);    %maximum diameter of the hole with valve
                
                t=tshim(n);        %shim thickness
                b=bshim(k);        %shim width

          
                Afvalve=pi*dfvalve^2/4; %area of the constantly opened valve
                AVvalve=pi*dVvalve^2/4; %maximum area of the hold with valve

                %fluid velocities
                Vf=V*Apiston/Afvalve;
             
                %compute the change in pressure
                dP=p*Vf.^2/2;
                                
                %compute the moment of inertia
                I=L*t^3/12;
                a=abs(L-b);

                %compute the variable area
                Avalve=(b*dP*AVvalve*L^3)/(3*E*I);            
            
                %total area of the valve
                Atotal=Afvalve + Avalve;
                          
                %compute the flow rate when valve is fully open
                Qopen= Cd*Atotal.*V;
            
                 P=dP.*Atotal;            
                 error=mean((P-Force).^2);           %error 
            
                 wd=max((P*L^3)/(3*E*I));      %Maximum  displacement
             
                flag=b*wd >max(Avalve);       %flag            
            
                if flag
                else                
                    d=[dfvalve,dVvalve,t,b,error];
                    results=vertcat(results,d);            
                end

            end
         end 
     end
  end


Fixed_valve_diameter=results(:,1)*1000;  %constant diameter of the valve
Maximum_valve_diameter=results(:,2)*1000;  %maximum diameter of the hold with valve
shim_thickness=results(:,3)*1000;  %shim thickness
shim_width=results(:,4)*1000;      %shim width
error=results(:,5);                %error

[~,idx]=min(error); %find the minimum error


%create Table
T=table(Fixed_valve_diameter,Maximum_valve_diameter,shim_thickness,shim_width,error)
writetable(T,'summary1.xlsx');

optimized_fixed_diameter=Fixed_valve_diameter(idx)*1e-3;
optimized_variable_diameter=Maximum_valve_diameter(idx)*1e-3;
optimized_shim_thickness=shim_thickness(idx)*1e-3;
optimized_shim_width=shim_width(idx)*1e-3;

optimized_fixed_Area=pi*optimized_fixed_diameter^2/4;
optimized_variable_Area=pi*optimized_variable_diameter^2/4;

%fluid velocities
Vf=V*Apiston/optimized_fixed_Area;

%compute the change in pressure
dP=p*Vf.^2/2;

%compute the moment of inertia
I=L*optimized_shim_thickness^3/12;

%compute the variable area
Avalve=(optimized_shim_width*dP*optimized_variable_Area*L^3)/(3*E*I); 

%total area of the valve
Atotal=optimized_fixed_Area + Avalve;

%actual force
Factual=dP.*Atotal;


figure(1)
plot(V,Factual,'LineWidth',1.5)
hold on
plot(V,Force,'LineWidth',1.5)

xlabel('Piston velocity [m/s]')
ylabel('Force [N]')
legend('Actual Force','Desired Force','lcn','northeast')

Velocity_piston=V; %velocity of piston
Velocity_fluid=Vf; %velocity of fluid
change_pressure=dP; %change in pressure
variable_area=Avalve; %Area of variable valve
Total_area=optimized_fixed_Area + Avalve; %Total area of the valve
System_Force=Factual; %actual Force
Desired_Force=Force;  %Desired Force
Error=(System_Force-Desired_Force).^2; %error

T2=table(Velocity_piston,Velocity_fluid,change_pressure,variable_area,Total_area,...
    System_Force,Desired_Force,Error)
writetable(T2,'summary2.xlsx');