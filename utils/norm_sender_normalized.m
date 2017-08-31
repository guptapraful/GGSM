% try to fix the range of z field
function [subband size_band] = norm_sender_normalized_gen(pyro,pind,Nsc,Nor,parent,neighbor,blSzX,blSzY,nbins)
%function [shape] = norm_sender_normalized_gen(pyro,pind,Nsc,Nor,parent,neighbor,blSzX,blSzY,nbins)

guardband = 16;
pyro = real(pyro);
Nband = size(pind,1)-1;
band = [1 3 6 8 9 11];
zRange = 15;

p = 1;
for scale=1:Nsc
    for orien=1:Nor
        nband = (scale-1)*Nor+orien+1; % except the ll
        %      if(prod(band-nband-1) ~=0)
        %            continue;
        %      end
        aux = pyrBand(pyro, pind, nband);
        [Nsy,Nsx] = size(aux);
        
        prnt = parent & (nband < Nband-Nor);   % has the subband a parent?
        BL = zeros(size(aux,1),size(aux,2),1 + prnt);
        BL(:,:,1) = aux;
        if prnt,
            auxp = pyrBand(pyro, pind, nband+Nor);
            %    if nband>Nor+1,     % resample 2x2 the parent if not in the high-pass oriented subbands.
            % 	   auxp = real(imenlarge2(auxp)); % this was uncommented
            auxp = real(imresize(auxp,2)); %
            %   end
            %  fprintf('parent band and size is %d %d %d \n',nband+Nor,Nsy,Nsx)
            BL(:,:,2) = auxp(1:Nsy,1:Nsx);
        end
        y=BL;
        [nv,nh,nb] = size(y);
        block = [blSzX blSzY];
        
        nblv = nv-block(1)+1;	% Discard the outer coefficients
        nblh = nh-block(2)+1;   % for the reference (centrral) coefficients (to avoid boundary effects)
        nexp = nblv*nblh;			% number of coefficients considered
        N = prod(block) + prnt; % size of the neighborhood
        
        Ly = (block(1)-1)/2;		% block(1) and block(2) must be odd!
        Lx = (block(2)-1)/2;
        if (Ly~=floor(Ly))|(Lx~=floor(Lx)),
            error('Spatial dimensions of neighborhood must be odd!');
        end
        Y = zeros(nexp,N);		% It will be the observed signal (rearranged in nexp neighborhoods)
        % Rearrange observed samples in 'nexp' neighborhoods
        n = 0;
        for ny=-Ly:Ly,	% spatial neighbors
            for nx=-Lx:Lx,
                n = n + 1;
                foo = shift(y(:,:,1),[ny nx]);
                foo = foo(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
                Y(:,n) = (foo(:));
            end
        end
        
        if prnt,	% parent
            n = n + 1;
            foo = y(:,:,2);
            foo = foo(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
            Y(:,n) = (foo(:));
        end
        
        %      including neighbor
        if neighbor,
            for neib=1:Nor
                if neib == orien
                    continue;
                end
                n=n+1;
                nband1 = (scale-1)*Nor+neib+1; % except the ll
                aux1 = pyrBand(pyro, pind, nband1);
                aux1 = aux1(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
                Y(:,n) = (aux1(:));
            end
        end
        
        C_x = innerProd(Y)/nexp;
        % C_x is positive definete covariance matrix
        [Q,L] = eig(C_x);
        % correct possible negative eigenvalues, without changing the overall variance
        L = diag(diag(L).*(diag(L)>0))*sum(diag(L))/(sum(diag(L).*(diag(L)>0))+(sum(diag(L).*(diag(L)>0))==0));
        C_x = Q*L*Q';
        
        o_c = aux(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
        o_c = (o_c(:));
        
        o_c = o_c - mean(o_c);
        
        o_c_p = auxp(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
        o_c_p = (o_c_p(:));
        
        o_c_p = o_c_p - mean(o_c_p);
        
        o_c_n = aux1;
        o_c_n = (o_c_n(:));
        
        o_c_n = o_c_n - mean(o_c_n);
        
        [hy rangeo] = hist(o_c,nbins);
        hy=hy/sum(hy);
        
        %         window = ones(blSzX, blSzY);
        %         window = window/sum(sum(window));
        %         [alpha1,ss]   = estimateggdparam(aux1(:));
        %         z = ((alpha1+1)*filter2(window, (abs(aux1)).^(alpha1),  'same')).^(1/(alpha1));
        %         z = z(:);
        
        tempY = (Y*pinv(C_x)).*Y/N;
        z1 = sqrt(sum(tempY,2));
        
        Y_old = Y;
        [sigma_y, s_y_new] = estimate_params_mv_ggd_mm(Y);
        z_new = (sum((Y*pinv(sigma_y)).*Y,2)).*(((s_y_new)./n).^(1./(s_y_new)));
        z_old = ones(size(z_new));
        s_y_old = 0;
        j = 1;
        max_iter = 50;
        while j <= max_iter & abs(s_y_new - s_y_old) > 1e-4 %norm(z_new - z_old)/numel(z_old) > 1e-5
            
            z_old = z_new;
            s_y_old = s_y_new;
            %Y = Y_old./repmat(sqrt(z_new), 1, n);
            o_c_norm = o_c./sqrt(z_new);
            o_c_norm = reshape(o_c_norm, nblv, nblh);
            o_c_p_norm = o_c_p./sqrt(z_new);
            o_c_p_norm = reshape(o_c_p_norm, nblv, nblh);
            n = 0;
            for ny=-Ly:Ly,	% spatial neighbors
                for nx=-Lx:Lx,
                    n = n + 1;
                    foo = shift(o_c_norm,[ny nx]);
                    Y(:,n) = (foo(:));
                end
            end
            if prnt,	% parent
                n = n + 1;
                foo = o_c_p_norm;
                Y(:,n) = (foo(:));
            end
            
            %      including neighbor
            if neighbor,
                for neib=1:Nor
                    if neib == orien
                        continue;
                    end
                    n=n+1;
                    nband1 = (scale-1)*Nor+neib+1; % except the ll
                    aux1 = pyrBand(pyro, pind, nband1);
                    aux1 = aux1(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
                    Y(:,n) = (aux1(:))./sqrt(z_new);
                end
            end
            
            [sigma_y, s_y_new] = estimate_params_mv_ggd_mm(Y);
            y_z = sum((Y_old*pinv(sigma_y)).*Y_old,2);
            z_new = (y_z.*(((s_y_new)./n).^(1./s_y_new)));
            norm_z(j) = norm(z_new - z_old)/numel(z_old);
            j= j+1;
            
        end
        z = sqrt(z_new);
        
        ind = find(z~=0);
        
        g_c =o_c(ind)./z(ind);
        
        
        %% consider the guardband
        
        g_c=reshape(g_c,nblv,nblh);
        gb = guardband/(2^(scale-1));
        g_c=g_c(gb+1:end-gb, gb+1:end-gb);
        size_band(p,:) = size(g_c);
        g_c= (g_c(:));
        g_c = g_c - mean(g_c);
        
        subband{p} = g_c;
        p = p+1;
    end
end
return;

