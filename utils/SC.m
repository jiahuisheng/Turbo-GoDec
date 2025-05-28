function  X_ica=SC(X_ica,j,no_lines,no_rows,no_bands)
% sparse X_ica
      T=X_ica;
      [~,idx]=sort(abs(T(:)),'descend');
      card=j*no_lines*no_rows;
      Num=no_lines*no_rows*no_bands;
      X_ica(idx(card+1:Num))=0;