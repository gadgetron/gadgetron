

function rho = simple_cg_nc_test(rhs, csm, co, w, D)
  if (isempty(D)), 
    D = ones(size(rhs));
  end
  
  r = D .* rhs;
  
  rr_0 = r(:)' * r(:);

  rr_1 = 0;
  rr = 0;
  
  rr_last = 1e10;
  
  rho = zeros(size(rhs));
  p = rho;
  
  for it=1:5,
      rr_1 = rr;
      rr = r(:)' * r(:);
      
      if (it == 1),
          p = r;
      else
          beta = rr/rr_1;
          p = p .* beta;
          p = p + r;
      end
      
      p1 = D .* p;
      %tmp_pp = p(:)' * p(:)

      q = mult_MH_M(p1, csm, co, w);
      %tmp_qq = q(:)' * q(:)

      q = D .* q;
      
      alpha = rr/(p(:)'* q(:));
      
      rho = rho + alpha * p;
      r = r - alpha*q;
      
      rr_rr0 = rr/rr_0
      
  end
  
  rho = D .* rho; 

return


function b = mult_MH_M(x, csm, co, w)
    coils = size(csm,length(size(csm)));
    data = itok(repmat(x,[1 1 size(csm,length(size(csm)))]) .* csm ,[1,2]);
    
    for c=1:coils,
        tmp = grid_data_bck(permute(co,[2 1]) .* size(csm,1),data(:,:,c),[0 0]);
        data_out(:,c) = tmp(:);
    end
    
    b = zeros(size(csm));
    for c=1:coils,
        tmp = grid_data(permute(co,[2 1]) .* size(csm,1), data_out(:,c), w(:), [size(csm,1) size(csm,2)], [0 0]);
        b(:,:,c) = tmp;
    end
    b = sum(ktoi(b,[1,2]) .* conj(csm),3);
return