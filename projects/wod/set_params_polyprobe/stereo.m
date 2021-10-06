function [stereotaxie]=stereo(X_bl,Y_bl,X_PA,Y_PA,X_PB,Y_PB)

stereotaxie=struct;

%deviation probe A
stereotaxie.x_pa= (X_PA/X_bl)*(((X_bl^2)+(Y_bl^2))^0.5);
stereotaxie.y_pa= cosd(atand(Y_bl/X_bl))*(Y_PA - Y_bl*(X_PA/X_bl));

%deviation probe B
stereotaxie.x_pb= (X_PB/X_bl)*(((X_bl^2)+(Y_bl^2))^0.5);
stereotaxie.y_pb= cosd(atand(Y_bl/X_bl))*(Y_PB - Y_bl*(X_PB/X_bl));

end