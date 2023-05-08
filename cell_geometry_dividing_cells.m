byhandR=1.447538462;
files=dir("*.tif");
j=1;
clear DATA
for i=1:length(files)
    aux=imread(strcat(files(i).folder,'/',files(i).name));
    aux=bwareafilt(aux==255,[500,1.5*10^7])
    [L, num] = bwlabel(aux);
    dx=1/9.2308;
    for k = 1 : num
        DATA.Cell(:,:,j)=ismember(L, k);
        Border=conv2(double(DATA.Cell(:,:,j)),[0 1 0; 1 -4 1; 0 1 0],'same')<0;
        [rbor, cbor]=find(Border==1);
        
        % HOW TO FIND SYMMETRY AXIS BY PCA
        [xx yy]=find(DATA.Cell(:,:,j)==1);
        pp=pca([xx, yy]);
        p(1)=pp(2,1)/pp(1,1);
        baric=[mean(xx) mean(yy)];
        p(2)=baric(2)-p(1)*baric(1);            
            
        x1=0; %point 1
        x2=0; %point 2
        d1=1e10; d2=1e10;
        for t=1:1:length(rbor)
            if F_DistancePL([rbor(t) cbor(t)],p)<d1
                d1=F_DistancePL([rbor(t) cbor(t)],p);
                x1=t;
            end
        end
        for t=1:1:length(rbor)
            if F_DistancePL([rbor(t) cbor(t)],p)<d2 && norm([rbor(x1) cbor(x1)]-[rbor(t) cbor(t)])>50
                d2=F_DistancePL([rbor(t) cbor(t)],p);
                x2=t;
            end
        end
        DATA.x1(j,:)=[rbor(x1) cbor(x1)]; %extreme point 1 of the symmetry axis
        DATA.x2(j,:)=[rbor(x2) cbor(x2)]; %extreme point 2 of the symmetry axis
                
        %CELL PROFILE ALONG THE SEMGENT x1...x2 (rotation axis)
        Fradius=[]; 
        Xradius=[]; 
        for b=1:length(rbor)
            [Fradius(b) Xradius(b)]=F_Distance3P(DATA.x1(j,:),DATA.x2(j,:),[rbor(b) cbor(b)]);
        end
        fp=find(Fradius>=0);
        fn=find(Fradius<=0);
        Xradiusp=Xradius(fp)*dx;
        [Xradiusp,Ip]=sort(Xradiusp,'ascend');
        Xradiusn=Xradius(fn)*dx;
        [Xradiusn,In]=sort(Xradiusn,'ascend');
        Fradiusp=Fradius(fp(Ip))*dx;
        Fradiusn=Fradius(fn(In))*dx;
        
        
        Fradiusp=[0 mean(Fradiusp(1:3)) (Fradiusp(1:end-4)+Fradiusp(2:end-3)+Fradiusp(3:end-2)+Fradiusp(4:end-1)+Fradiusp(5:end))/5 mean(Fradiusp(end-2:end)) 0];
        Fradiusn=[0 mean(Fradiusn(1:3)) (Fradiusn(1:end-4)+Fradiusn(2:end-3)+Fradiusn(3:end-2)+Fradiusn(4:end-1)+Fradiusn(5:end))/5 mean(Fradiusn(end-2:end)) 0];
        DFradiusp=[(Fradiusp(2)-Fradiusp(1))/(Xradiusp(2)-Xradiusp(1)), (Fradiusp(3:end)-Fradiusp(1:end-2))./(Xradiusp(3:end)-Xradiusp(1:end-2)), (Fradiusp(end)-Fradiusp(end-1))/(Xradiusp(end)-Xradiusp(end-1))];
        DFradiusn=[(Fradiusn(2)-Fradiusn(1))/(Xradiusn(2)-Xradiusn(1)), (Fradiusn(3:end)-Fradiusn(1:end-2))./(Xradiusn(3:end)-Xradiusn(1:end-2)), (Fradiusn(end)-Fradiusn(end-1))/(Xradiusn(end)-Xradiusn(end-1))];
        DFradiusp(find(abs(DFradiusp)>10))=10*sign(DFradiusp(find(abs(DFradiusp)>10)));
        DFradiusn(find(abs(DFradiusn)>10))=10*sign(DFradiusn(find(abs(DFradiusn)>10)));
                  
        alfa=0.7;
        DATA.GEO(1,j)=mean(Fradiusp(find(Fradiusp>alfa*max(Fradiusp))))/2+mean(-Fradiusn(find(-Fradiusn>alfa*max(-Fradiusn))))/2; %radius
        DATA.GEO(2,j)=norm(DATA.x1(j,:)-DATA.x2(j,:))*dx; %length
        DATA.GEO(3,j)=pi*trapz(Xradiusp,Fradiusp.*sqrt(1+DFradiusp.^2))-pi*trapz(Xradiusn,Fradiusn.*sqrt(1+DFradiusn.^2)); %area rotation
        DATA.GEO(4,j)=0.5*pi*trapz(Xradiusp,Fradiusp.^2)+0.5*pi*trapz(Xradiusn,Fradiusn.^2); %volume rotation
        DATA.GEO(5,j)=2*pi*DATA.GEO(1,j)*DATA.GEO(2,j); %area=2piRL
        DATA.GEO(6,j)=pi*DATA.GEO(1,j)^2*DATA.GEO(2,j)-2/3*pi*DATA.GEO(1,j)^3; %volume=piR^2L-2/3piR^3
        j=j+1;
    end
end
DATA.GEO(7,:)=2*pi*mean(DATA.GEO(1,:))*DATA.GEO(2,:); %area=2piRL
DATA.GEO(8,:)=pi*mean(DATA.GEO(1,:))^2*DATA.GEO(2,:)-2/3*pi*mean(DATA.GEO(1,:))^3; %volume=piR^2L-2/3piR^3
DATA.GEO(9,:)=2*pi*byhandR*DATA.GEO(2,:); %area=2piRL
DATA.GEO(10,:)=pi*byhandR^2*DATA.GEO(2,:)-2/3*pi*byhandR; %volume=piR^2L-2/3piR^3
labels=["Mean width","Length","area(rot)","vol(rot)","area","vol","area_mr","vol_mr","areabyhandR","volbyhandR"];
for i=1:6
subplot(4,2,i)
histogram(DATA.GEO(i,:),10)
title(labels(i))
end
csvwrite("data.csv",DATA.GEO')
            %%FUNCTIONS---------------------------------------------------------------

function  d = F_DistancePL(x0,p)    %distance of x0 from line y=ax+b
    a=p(1); b=p(2);
    d=abs(a*x0(1)-x0(2)+b)/sqrt(a^2+1);
end

function  [d, d1] = F_Distance3P(x1,x2,x3)     %distance of x3 from the x1-x2 line
    if norm(x2-x3)~=0 && norm(x1-x3)~=0
        v12=double(x2(:)-x1(:));
        v13=double(x3(:)-x1(:));
        d1=real((v13'*v12)/norm(v12));
        d=real(norm(v13)*sqrt(1-((v13'*v12)/(norm(v12)*norm(v13)))^2));
        if x1(1)~=x2(1)
            p=polyfit([x1(1), x2(1)],[x1(2), x2(2)],1);
            d=d*sign(x3(2)-x3(1)*p(1)-p(2));        
        else       
            d=d*sign(x3(1)-x1(1));
        end
    else
        d=0;
        if norm(x2-x3)==0
            d1=norm(x2-x1);
        else
            d1=0;
        end
    end
end

function [mo ma]=F_SumClustering (x,y,a)  %grouping pixels according to their position along the symmetry axis
    [o a0]=hist(x,[(3*a(1)-a(2))/2, a, (3*a(end)-a(end-1))/2]);
    o=o(2:end-1);
    q=find(x<=(a(1)/2+a(2)/2));
    ma(1)=mean(x(q));
    mo(1)=sum(y(q));
    for j=2:length(a)-1
        x1=(a(j)+a(j-1))/2;
        x2=(a(j)+a(j+1))/2;
        q=find(((x>x1).*(x<=x2))==1);
        ma(j)=mean(x(q));
        mo(j)=sum(y(q));
    end
    q=find(x>(a(end)/2+a(end-1)/2));
    ma(length(a))=mean(x(q));
    mo(length(a))=sum(y(q));

    t=find(isnan(mo+ma)==0);
    mo=mo(t);
    ma=ma(t);
end

function [m s n c]=F_mean_variance_cake (X,f) %fitting of the Cdr2 profile along the symmetry axis
    X=X(:);
    f=f(:);
    f=(f-0.5*(mean(f(1:20))+mean(f(end-19:end))));
    f=f.*(f>0);
%     f=f(:)/sum(f);
    m=sum(X.*f)/sum(f);
    s=sum((X.^2).*f)/sum(f)-m^2;
    p0=[m sqrt(s) 1 0];
    
    F=@(p)F_error(X,f,p);
    options=optimset('MaxFunEvals',1e5);
    [popt fval]=fminsearch(F,p0,options);
    m=popt(1);
    s=popt(2)^2;
    n=abs(popt(3));
    c=abs(popt(4));  
end

function e=F_error (X,Y,p)
    f=abs(p(3))*pdf('normal',X,p(1),abs(p(2)));
    R=(max(X)-min(X))/2;
    xc=(max(X)+min(X))/2;
    c=abs(p(4))*real(sqrt(R^2-(X-xc).^2));
    e=norm(Y-f-c);
end