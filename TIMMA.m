% The fomer version using maximum averaging to predict the M value. In this version
% also implement the minimum averaging to get the upper boundary


% Tao Xu
function [M,err] = TIMMA(ds, S, y_vals, a, LOO)  
                                            %ds = number of drugs
                                            %S = kinase profiles of ds drugs
                                            %a = number of kinases per S
                                            %y_vals = specificity
                                            %LOO = whether to use LOO
                                            
%  tic()

% version 4 of obtaining decimal gray code
[rows,cols,G_dec] = graycode(a);

% We build an inhibition matrix for each drug using multi-dimensional arrays
% IM_d: cells which are directly observed
% IM_i: cells which are indirectly observed, e.g the subsets of which are
% directly observed in S
% IM_o: cells which are overly observed, e.g the supersets of which are
% directly observed in S
IM_d = NaN(rows,cols,ds);
IM_i = NaN(rows,cols,ds);
IM_o = NaN(rows,cols,ds);

% index of the drugs in G
S_r = zeros(1,ds); 
for i = 1:ds 
    
    dec_val = bin2dec(char(S(i,:)+48));
    tmp = G_dec==dec_val;
	% change logical from  to double
    tmp = +tmp;   % equivalent of S(i,:)
   
    [super_set,sub_set]=binary_set(S(i,:));
	%superset of S(i,:)
    tmp2=+ismember(G_dec,super_set); 
	% absolute subset of S(i,:)
    tmp3=+ismember(G_dec,sub_set); 

    
    S_r(i) = find(tmp);

    tmp(tmp==0)=NaN;
    tmp2(tmp2==0)=NaN;
    tmp3(tmp3==0)=NaN;
  

    IM_d(:,:,i) = tmp.*y_vals(i);
    IM_i(:,:,i) = tmp2.*y_vals(i);
    IM_o(:,:,i) = tmp3.*y_vals(i);

end

% the brutal averaging on direct observed cells
M_d = nansum(IM_d,3)./sum(~isnan(IM_d),3); 

[M_i, ds_hat] = max(IM_i,[],3);
[M_o, min_o] = min(IM_o,[],3);

% pinpoint the cells which needs maximization averaging
% e.g. the cells which are supersets of observed drugs
cel = find(isnan(M_d)&~isnan(M_i));

%cel = find(isnan(M_d));
if (isempty(cel))
else
    for i=cel'
		% find the index for the cell
        [i_c,j_c] = ind2sub([rows,cols],i); 
		% the drug sets which are subsets of the cell
        tmp5 = squeeze(~isnan(IM_i(i_c,j_c,:))); 
		% the drug which achieves max sensitivity
        index = ds_hat(i); 
		 % the corresponding gray code G{i_g,j_g}
        [i_g,j_g] = find(IM_d(:,:,index)>=0);
    	%[i_g,j_g] = find(~isnan(IM_d(:,:,index)));
        % find the supersets of S(index,:) in S that has smaller sensitivity
        tmp4 = squeeze(IM_o(i_g,j_g,:) < M_i(i));
        tmp6 = find(tmp5 & tmp4);
        k = 1;
        if ~isempty(tmp6)
			% max-averaging
            for j=tmp6' 
                M_i(i) = (M_i(i)*k+y_vals(j))/(k+1);
                k = k+1;
            end
        end
    end
end

% pinpoint the cells which needs minimization averaging
% e.g. the cells which are subsets of observed drugs

cel2 = find(isnan(M_d)&~isnan(M_o));

if (isempty(cel2))
    
else
    for i=cel2'
        [i_c,j_c] = ind2sub([rows,cols],i); % find the index for the cell
        tmp5 = squeeze(~isnan(IM_o(i_c,j_c,:))); % the drug sets which are supersets of the cell
        index = min_o(i); % the drug which achieves max sensitivity
        [i_g,j_g] = find(IM_d(:,:,index)>=0); % the corresponding gray code G{i_g,j_g}
        % find the subsets of S(index,:) in S that has bigger sensitivity
        tmp4 = squeeze(IM_i(i_g,j_g,:) > M_o(i));
        tmp6 = find(tmp5 & tmp4);
        k = 1;
        if ~isempty(tmp6)
            for j=tmp6' % min-averaging
%                 disp('------- min average running now ------');
                M_o(i) = (M_o(i)*k+y_vals(j))/(k+1);
                k = k+1;
            end
        end
    end
end



M = M_d;
M(cel) = M_i(cel);
M(cel2) = M_o(cel2);
% cels that not only have lower boundery and also have upper boundary 
average_index = intersect(cel,cel2);
M(average_index) = (M_i(average_index)+M_o(average_index))/2;





err = NaN(1,ds);
if LOO==0
    % test error
    for i=1:ds
        err(i)=abs(M(S_r(i))-y_vals(i));
    end
else
    for i = 1:ds
        % remove drug i
       IM_d_LOO = IM_d;
       IM_i_LOO = IM_i;
       IM_o_LOO = IM_o;
       y_vals_LOO = y_vals;
       
       IM_d_LOO(:,:,i) = [];
       IM_i_LOO(:,:,i) = [];
       IM_o_LOO(:,:,i) = [];
       y_vals_LOO(i) = [];
       
       M_d_LOO = nansum(IM_d_LOO,3)./sum(~isnan(IM_d_LOO),3);
       M_LOO = M_d_LOO;
       [M_i_LOO,ds_hat] = max(IM_i_LOO,[],3);
       [M_o_LOO,min_hat] = min (IM_o_LOO,[],3);
       
       cel = find(isnan(M_d_LOO)&~isnan(M_i_LOO));
       cel2 = find(isnan(M_d_LOO)&~isnan(M_o_LOO));
       
       j_max = find(cel==S_r(i)); 
       j_min = find(cel2==S_r(i)); 
       [avg_ind max_ind min_ind]= intersect(j_max,j_min);
       j_max(max_ind)=[];
       j_min (min_ind)=[];
     
       if (j_max)
           [i_c,j_c] = ind2sub([rows,cols],cel(j_max));
			
           tmp5 = squeeze(~isnan(IM_i_LOO(i_c,j_c,:))); 
           index = ds_hat(cel(j_max)); % the drug which achieves max sensitivity
           [i_g,j_g] = find(IM_d_LOO(:,:,index)>=0); % the corresponding gray code G{i_g,j_g}

           % find the supersets of S(index,:) in S that has smaller
           % sensitivity than maximal
           tmp4 = squeeze(IM_o_LOO(i_g,j_g,:) < M_i_LOO(cel(j_max)));
           tmp6 = find(tmp5 & tmp4);
           k = 1;
           if ~isempty(tmp6)
               for m=tmp6' % max-averaging
                   M_i_LOO(cel(j_max)) = (M_i_LOO(cel(j_max))*k+y_vals_LOO(m))/(k+1);
                   k = k+1;
               end
           end
            err(i) =  abs(M_i_LOO(S_r(i))-y_vals(i));
       
       
       else if j_min
           [i_c,j_c] = ind2sub([rows,cols],cel2(j_min));
           tmp5 = squeeze(~isnan(IM_o_LOO(i_c,j_c,:))); % the drug sets which are superset of the cell
           index = min_hat(cel2(j_min)); % the drug which achieves min sensitivity
           [i_g,j_g] = find(IM_d_LOO(:,:,index)>=0); % the corresponding gray code G{i_g,j_g}

           % find the subsets of S(index,:) in S that has bigger
           % sensitivity than minimal
           tmp4 = squeeze(IM_i_LOO(i_g,j_g,:) > M_o_LOO(cel2(j_min)));
           tmp6 = find(tmp5 & tmp4);
           k = 1;
           if ~isempty(tmp6)
               for m=tmp6' % min-averaging
                   M_o_LOO(cel2(j_min)) = (M_o_LOO(cel2(j_min))*k+y_vals_LOO(m))/(k+1);
                   k = k+1;
               end
           end
             err(i) =  abs((M_o_LOO(S_r(i)))-y_vals(i));
           else if avg_ind
				[i_c,j_c] = ind2sub([rows,cols],cel(avg_ind));
                % the drug sets which are subsets of the cell
				tmp5_i = squeeze(~isnan(IM_i_LOO(i_c,j_c,:))); 
				% the drug sets which are superset of the cell
				tmp5_o = squeeze(~isnan(IM_o_LOO(i_c,j_c,:))); 
				% the drug which achieves max sensitivity
                ind_i = ds_hat(cel(avg_ind)); 
                ind_o = min_hat(cel2(avg_ind));
 				% the corresponding gray code G{i_g,j_g}
                [i_g,j_g] = find(IM_d_LOO(:,:,ind_i)>=0);
                [i_g_o,j_g_o] = find(IM_d_LOO(:,:,ind_o)>=0);
                % find the supersets of S(index,:) in S that has smaller
                % sensitivity than maximal
                tmp4_i = squeeze(IM_o_LOO(i_g,j_g,:) < M_i_LOO(cel(avg_ind)));
                tmp4_o = squeeze(IM_i_LOO(i_g_o,j_g_o,:) > M_o_LOO(cel2(avg_ind)));
                tmp6_i = find(tmp5_i & tmp4_i);
                tmp6_o = find(tmp5_o & tmp4_o);
                k = 1;
                if ~isempty(tmp6_i)
                   	for m=tmp6_i' % max-averaging
                    	M_i_LOO(cel(avg_ind)) = (M_i_LOO(cel(avg_ind))*k+y_vals_LOO(m))/(k+1);
                        k = k+1;
                    end
                end
                k = 1;
                if ~isempty(tmp6_o)
                   for m=tmp6_o' % min-averaging
                       M_o_LOO(cel(avg_ind)) = (M_o_LOO(cel(avg_ind))*k+y_vals_LOO(m))/(k+1);
                       k = k+1;
                   end
                end
                   err(i) = abs((M_i_LOO(S_r(i))+M_o_LOO(S_r(i)))/2-y_vals(i));
               else
                   err(i) = abs(M_LOO(S_r(i))-y_vals(i));
            	end
         	end
       end       
    end
end

%   toc()
