function [ times_rs, xpRs, ypRs, alts, v_nom] = preprocess_latlon_data( plane_t )

        times=[plane_t.time-plane_t.time(1)]/60;
        times_rs=times;
        lons=interp1(times,plane_t.ll(1,:)',times_rs);
        lats=interp1(times,plane_t.ll(2,:)',times_rs);
        alts=interp1(times,plane_t.alt(1,:)',times_rs);

        %         figure(1)
        %         plot(lons,lats), hold on
        %         plot(lons(1),lats(1),'o','Color',[1 0 0]), hold off
        %         figure(2)
        %         subplot(3,1,1)
        %         plot(times_rs,lons,'-o')
        %         subplot(3,1,2)
        %         plot(times_rs,lats,'-o')
        %         subplot(3,1,3)
        %         plot(times_rs,alts,'-o')

        %------------------------------------------------------------------------
        %Properties of Earth, the oblate spheroid thing
        num_pts=length(lons);
        lat_center=median(lats);
        lon_center=median(lons);
        a=3443.917;   %Radius at equator
        b=3432.370;   %Radius at the poles
        R_num=(a^2*cos(lats)).^2+(b^2*sin(lats)).^2;
        R_den=(a*cos(lats)).^2+(b*sin(lats)).^2;
        RadEarth=median(sqrt(R_num./R_den));   %Radius at each longitude, latitude pair

        cosc=sind(lat_center)*sind(lats)+cosd(lat_center)*cosd(lats).*cosd(lons-lon_center);
        xp=RadEarth*cosd(lats).*sind(lons-lon_center)./cosc;
        yp=RadEarth*cosd(lat_center)*sind(lats)-sind(lat_center)*cosd(lats).*cosd(lons-lon_center)./cosc;
        %------------------------------------------------------------------------
        %Rotate Frame
        xp=xp-xp(1);
        yp=yp-yp(1);
        dx=xp(end)-xp(1);
        dy=yp(end)-yp(1);
        Rang=  atan2(dy,dx);
        R=[cos(Rang) sin(Rang);-sin(Rang) cos(Rang)];
        xy=R*[xp; yp];
        xpR=xy(1,:)';
        ypR=xy(2,:)';
    
        %------------------------------------------------------------------------
    
        xpRs=smooth(times_rs,xpR,num_pts/40,'lowess');
        ypRs=smooth(times_rs,ypR,num_pts/40,'lowess');

        cdist=cumsum((diff(xpRs).^2+diff(ypRs).^2).^.5);
        dtime=times(end)-times(1);
        v_nom=cdist(end)/dtime;


        
end
