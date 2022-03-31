function [Pos , Kvsp , Kasp , Kjsp] = Cubic_Bspline(tt , signal)

spfactor=3;
k=4;


for i = 1:size(signal,2)
    Kpos=signal(:,i);
    s=size(Kpos);
    sp1x=spap2(fix(s(1)/spfactor),k,tt,Kpos);
    Kjsp(:,i)=fnval(fnder(sp1x,3),tt(:,1));
    Kasp(:,i)=fnval(fnder(sp1x,2),tt(:,1));
    Kvsp(:,i)=fnval(fnder(sp1x,1),tt(:,1));
    Pos(:,i)=fnval(fnder(sp1x,0),tt(:,1));
end

% sp1x=spap2(fix(s(1)/spfactor),k,Kpos(:,2)/1000,Kpos(:,30));
% Kasp(:,30)=fnval(fnder(sp1x,2),Kasp(:,2)/1000);
% Kvsp(:,30)=fnval(fnder(sp1x,1),Kasp(:,2)/1000);
% sp1z=spap2(fix(s(1)/spfactor),k,Kpos(:,2)/1000,Kpos(:,32));
% Kasp(:,32)=fnval(fnder(sp1z,2),Kasp(:,2)/1000);
% Kvsp(:,32)=fnval(fnder(sp1z,1),Kasp(:,2)/1000);
