clear all
clc
close all

%%Author: Alejandra Vesga  and José David Sanabria Gómez
%%date: May 20th 2016

%Very simple code to calculate the cross-rock distances for each muon direction, considering in detail the topography around each volcano that is available in the global digital elevation model of the Earth at NASA Shuttle Radar Topography Mission (SRTM)


%%parameters

conv=111.1949;                                                            % km por degree, %Average radius of the earth: 6371km.
conv=conv*1000;                                                           % in meters
dx=3/3600*conv;

X=readhgt(4:5,-76:-75,'merge','interp','crop',[4.47,4.505,-75.4,-75.37]); % Download topographic data from Nasa


lat_t=X.lat;                                                              % in geographic coordinates.
lon_t=X.lon;
z_t=X.z;
z_t=double(z_t);


dato= [-75.381092, 4.492298,207,237,1, 2000,6]                            %localization observation point latitude, longitude

%%let's work

for n=1:1

close all
punto = num2str(dato(n,5));
figure(1)
imagesc(lon_t,lat_t,z_t)
colorbar
colormap jet
title(['Location Geographical Coordinates Observation Point',num2str(dato(n,5))])
xlabel('Longitude')
ylabel('Latitude')
set(gca,'YDir','normal')


G0=[dato(n,1),dato(n,2)];                                                  %point observation in geographic coordinates.
zg= interp2(lon_t,lat_t,z_t,G0(1),G0(2));

zg=double(zg);
G0_0=[G0(1) G0(2) zg];                                                     %point observation 3D
hold on
plot3(G0(1),G0(2),zg,'o','MarkerEdgeColor','r','MarkerFaceColor','b','MarkerSize',3);
hold on
str = 'Relieve2DGeoCOLCIENCIAS';
nombre = strcat(str,punto,'.pdf');
saveas(gcf,nombre)


str=' Redefiniendo sistema de coordenadas...';
disp(str)


lat_t_0=lat_t - min(lat_t);                                                %The reference system is redefined.
lon_t_0=lon_t - min(lon_t);

figure(2)
imagesc(lon_t_0,lat_t_0,z_t)
colorbar
colormap jet
title(['Location Local Coordinates Observation Point',num2str(dato(n,5))])
xlabel('Longitude')
ylabel('Latitude')
set(gca,'YDir','normal')

P0_lon=G0(1)-min(lon_t);
P0_lat=G0(2)-min(lat_t);
P0_0=[P0_lon, P0_lat, zg];                                                 %point observation 3D is redefined
hold on
plot3(P0_lon,P0_lat,zg,'o','MarkerEdgeColor','r','MarkerFaceColor','y','MarkerSize',3);
str = 'LocalesCOLCIENCIAS';
nombre = strcat(str,punto,'.pdf');
saveas(gcf,nombre)
%%

%to meters
lon_t_m=lon_t_0*conv;                                                      %The reference system is redefined
lat_t_m=lat_t_0*conv;                                                      
P0_m=[P0_lon, P0_lat, zg/conv]*conv; 

% interpolation...

lon_t_m_int=min(lon_t_m):dx/2:max(lon_t_m);
lat_t_m_int=min(lat_t_m):dx/2:max(lat_t_m);
[XX,YY] = meshgrid(lon_t_m_int,lat_t_m_int);
zint = interp2(lon_t_m,lat_t_m,z_t,XX,YY);
zint = double(zint);

figure(3)
surface(lon_t_m_int,lat_t_m_int,zint)
colorbar
colormap jet
zlim([1500 4000])
shading interp
xlabel('X (m)')
ylabel('Y (m)')
zlabel('altitud (m)')
az=45;
el=45;
view(az,el)
set(gca,'YDir','normal')
str = 'metrosCOLCIENCIAS';
nombre = strcat(str,punto,'.pdf');
nombre1 = strcat(str,punto,'.fig');
hold on
plot3(P0_m(1),P0_m(2),zg,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3);
hold on
str=' Generando grilla de referencia...';
disp(str)


%% Calculating distances
num=0;
polar=dato(n,3):1:dato(n,4);
cenit_inv=dato(n,7):1:30;
npolar = length(polar);
ncenit = length(cenit_inv);
nlong  = npolar*ncenit;
nh=dato(n,6);                                                               %Number of steps on the line
h=1;                                                                        %Step size on the line...

long = zeros(nlong,3);                                                      %Arrangement with the distance traveled for each pair of angles

for k=polar;                                                                %For the polar angle measured from the x-axis
    for j=cenit_inv;                                                        %para cenit inverso medido desde la horizontal (horizonte)
        num = num+1;                                                        %contador para almacenar las longitudes recorridas..
        
        VV  = [cosd(k), sind(k), tand(j)];                                  
        V   =  VV/norm(VV);                                                 %Unit vector that determines the ray.
        
        len = 0;
        clear P topo
        for i=1:nh
            P  = P0_m + i*h*V;                                              %The many points on the straight
            i_lat = floor(P(1) /(dx/2))+1;
            i_lon = floor(P(2) /(dx/2))+1;
            topo = zint(i_lon,i_lat);
            if (topo > P(3))                                                %If the topography is above the line...
               len = len + h;                                               %The accumulated length traveled within the volcano...
                
            end
        end
        long(num,1)= k;                                                     %The vector of lengths as a function of the polar angle and the inverse centit
        long(num,2)= j;                                                     
        long(num,3)=len;                                                    
        hold on
        
        if(mod(k,1)==0 )%& j==2)mod(j,2)==0)
            hold on
            line([P0_m(1), P(1)],[P0_m(2), P(2)],[P0_m(3), P(3)],'Marker','.','LineStyle','-','Color',[0 0 1],'LineWidth',1)
           
        end
    end
           
end

%%
saveas(gcf,nombre)



%%


dato1=reshape(long(:,3),length(cenit_inv),[]);
figure(4)
imagesc(polar,cenit_inv,fliplr(dato1))
set(gca,'xDir','reverse','YDir','normal')
colorbar
colormap jet
axis ([min(polar) max(polar) min(cenit_inv) max(cenit_inv)])
title(['Length L(m) of the rays crossing the lava dome in point',num2str(dato(n,5))])
xlabel('theta (degrees)')
ylabel('elevation(degrees)')
str = 'DistanciaCOLCIENCIAS';
nombre = strcat(str,punto,'.pdf');
print(nombre,'-depsc')
saveas(gcf,nombre)

str = 'distancias';
nombre1 = strcat(str,punto,'.txt');
save (nombre1,'long','-ASCII','-TABS')


end



