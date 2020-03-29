function d = dist_cal(p_XZ,q_YZ,distype,dimtype)
    if dimtype == 1 %dimtype = 2 => X; dimtype = 1 => Z;
        dim1 = 1; dim2 = 2;
    else
        dim1 = 2; dim2 = 1;
    end
    pX_labswi = sum(p_XZ,dim1);
    pY_labswi = sum(q_YZ,dim1);
    pXY_mean = (pX_labswi + pY_labswi)/2;
   if distype == 1 %Jesen-Shannon Divergence between X and Y;
       d = 0.5*sum(pX_labswi.*log(pX_labswi./pXY_mean),dim2) + 0.5*sum(pY_labswi.*log(pY_labswi./pXY_mean),dim2);
   elseif distype == 2 %KL divergence over Y and X
       d = sum(pX_labswi.*log(pX_labswi./pY_labswi));
   else %KL divergence between YZ and XZ
       d = sum(sum(p_XZ.*log(p_XZ./q_YZ)));       
   end
  
        