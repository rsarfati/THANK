% load DataTHANK;
% y = DataTHANK;
% T = length(y);

[G1,C,impact,eu,SDX,zmat,NY,NX]=modelTHANK(postmode');

hor=[0:4];
CORR=zeros(7,7,length(hor));
CORRmod=zeros(7,7,length(hor));
countj=0;
for j=hor
    countj=countj+1;
    CORR(:,:,countj)=corr(y(j+1:T,:),y(1:T-j,:));
    Vobs=(zmat*dlyapgio(G1,impact*SDX*SDX'*impact')*zmat');
    CORRmod(:,:,countj)=zmat*[(G1^j)*dlyapgio(G1,impact*SDX*SDX'*impact')]*zmat'./(sqrt(diag(Vobs))*sqrt(diag(Vobs))');
end

figure(2)
countplot=0;
for i=1:7
    for j=1:7
        countplot=countplot+1;
        subplot(7,7,countplot); plot(hor,squeeze(CORR(i,j,:)),'LineWidth',2.5); hold on;
        plot(hor,squeeze(CORRmod(i,j,:)));
        %axis tight
        axis([0 4 -1 1])
        
        if i==1 & j==1; ylabel('dY_{t}','FontWeight','bold');
        elseif i==2 & j==1; ylabel('dC_{t}','FontWeight','bold');
        elseif i==3 & j==1; ylabel('dI_{t}','FontWeight','bold');
        elseif i==4 & j==1; ylabel('h_{t}','FontWeight','bold');
        elseif i==5 & j==1; ylabel('dW_{t}','FontWeight','bold');
        elseif i==6 & j==1; ylabel('p_{t}','FontWeight','bold');
        elseif i==7 & j==1; ylabel('R_{t}','FontWeight','bold');
        end
        
        if i==1 & j==1; title('dY_{t-k}');
        elseif i==1 & j==2; title('dC_{t-k}');
        elseif i==1 & j==3; title('dI_{t-k}');
        elseif i==1 & j==4; title('h_{t-k}');
        elseif i==1 & j==5; title('dW_{t-k}');
        elseif i==1 & j==6; title('p_{t-k}');
        elseif i==1 & j==7; title('R_{t-k}');
        end
    end
end
sgtitle('cross-correlogram')
