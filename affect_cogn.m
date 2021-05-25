%% reproducibility script of cognition-affect project
% Lina Schaare, May 2021

%% load utils
dir = []; % add data directory
P = [dir 'utils/surfstat/']; 
addpath(genpath(P));
addpath([P, 'surfstat_tutorial/surfstat']);
addpath([P, 'FreeSurfer5.3']);
addpath(genpath([dir 'utils/BrainSpace-0.1.1/matlab/']))
addpath(genpath([dir 'utils/BrainSpace-0.1.1/gifti-master/'])) 
addpath(genpath([dir 'utils/cifti-matlab-master']))

addpath([dir 'utils/cbrewer/'])
addpath([dir 'utils/misc/'])

DPATH = [dir 'data/HCP/'];                          
MASKPATH = [P, 'FreeSurfer5.3/fsaverage5/label/'];
RPATH = [dir 'affect_cognition/results/FIG/'];

%% load HCP data
HCP=readtable([DPATH, 'unrestricted_hlschaare_5_21_2020_5_11_56.csv']);
HCP_r=readtable([DPATH, 'RESTRICTED_hlschaare_5_21_2020_5_14_38.csv']);

% behavioural measures of interest
list_of_vars = [HCP.CogTotalComp_Unadj, HCP.CogFluidComp_Unadj, HCP.CogCrystalComp_Unadj,...
    HCP.LifeSatisf_Unadj,HCP.MeanPurp_Unadj,HCP.PosAffect_Unadj,...
    HCP.AngAffect_Unadj,HCP.FearAffect_Unadj,HCP.Sadness_Unadj];

%% load brain measures
for standard200_7 = 1
    
load_new = 0; % 0 = load pre-saved data, 1 = load all data from HCP repository
    
if load_new == 0
    
    HCP200_CT = csvread([dir 'data/HCP_fromSofie/CT_200_7_nevi.csv']);
    area200 = csvread([dir 'data/HCP_fromSofie/SA_200_7_nevi.csv']);
    
    % annotate parcels
    for parcels200 = 1
        [vertices, label, colortablel] = ...
            fs_read_annotation([MASKPATH 'lh.Schaefer2018_200Parcels_7Networks_order.annot']);
        parcel_left = label;
        label_left = label;
        for i = 1:size(colortablel.table, 1)
            mycode = colortablel.table(i,5);
            parcel_left(find(parcel_left == mycode)) = i;
        end
        
        [vertices, label, colortabler] = ...
            fs_read_annotation([MASKPATH 'rh.Schaefer2018_200Parcels_7Networks_order.annot']);
        parcel_right = label;
        label_right = label;
        for i = 1:size(colortabler.table, 1)
            mycode = colortabler.table(i,5);
            parcel_right(find(parcel_right == mycode)) = i;
        end
        
        parcels200 = [parcel_left; parcel_right+1000];
        parcels200 = parcels200';
        
        names200 = [colortablel.struct_names(2:end);colortabler.struct_names(2:end)];
        
    end
    
else        
     % load cortical thickness  
%      HCP200 = ft_read_cifti([dir 'utils/misc/Schaefer2018_200Parcels_7Networks_order.dlabel.nii'],'mapname','array');
%      
%      CTX_fs32k = zeros(length(HCP_r.Subject),64984);
%      
%     for i = 1:length(HCP_r.Subject)
%         i
%         try
%             L = gifti(['/Volumes/BnB3/BnB1/Raw_Data_nonBIDS/HCP/', num2str(HCP_r.Subject(i)),'/MNINonLinear/fsaverage_LR32k/', num2str(HCP_r.Subject(i)),'.L.thickness.32k_fs_LR.shape.gii'])
%             R = gifti(['/Volumes/BnB3/BnB1/Raw_Data_nonBIDS/HCP/', num2str(HCP_r.Subject(i)),'/MNINonLinear/fsaverage_LR32k/', num2str(HCP_r.Subject(i)),'.R.thickness.32k_fs_LR.shape.gii'])
%             
%             CTX_fs32k(i,:) = [L.cdata;R.cdata];
%         catch
%         end
%     end 
%     
%     % parcellate CT data
%     HCP200_CT= zeros(length(HCP_r.Subject),200);
%     for i = 1:200
%         HCP200_CT(:,i) = trimmean(CTX_fs32k(:,find(HCP200.dlabel==i))',10)';
%     end
%     
%     % annotate parcels
%     for parcels200 = 1
%         [vertices, label, colortablel] = ...
%             fs_read_annotation([MASKPATH 'lh.Schaefer2018_200Parcels_7Networks_order.annot']);
%         parcel_left = label;
%         label_left = label;
%         for i = 1:size(colortablel.table, 1)
%             mycode = colortablel.table(i,5);
%             parcel_left(find(parcel_left == mycode)) = i;
%         end
%         
%         [vertices, label, colortabler] = ...
%             fs_read_annotation([MASKPATH 'rh.Schaefer2018_200Parcels_7Networks_order.annot']);
%         parcel_right = label;
%         label_right = label;
%         for i = 1:size(colortabler.table, 1)
%             mycode = colortabler.table(i,5);
%             parcel_right(find(parcel_right == mycode)) = i;
%         end
%         
%         parcels200 = [parcel_left; parcel_right+1000];
%         parcels200 = parcels200';
%         
%         names200 = [colortablel.struct_names(2:end);colortabler.struct_names(2:end)]
%         
%     end
%      
%     % load surface area and parcellate
%     % use the unsmoothed data for the parcels
%     isthere_ct = zeros(size(HCP_r.Subject));
%     namesct_left = strcat('/Volumes/BnB_TEMP/Sofie/Genetics/2018.CIVIT_FS/FS_surf/', num2str(HCP_r.Subject), '_lh2areaj_fsaverage_1.mgh')
%     namesct_right = strcat('/Volumes/BnB_TEMP/Sofie/Genetics/2018.CIVIT_FS/FS_surf/', num2str(HCP_r.Subject), '_rh2areaj_fsaverage_1.mgh')
%     
%     AREA = zeros(length(HCP_r.Subject),size(SW.coord,2));
%     for i = 1:length(HCP_r.Subject)
%         try
%             AREA(i,1:10242)     = SurfStatReadData1(namesct_left(i,:));
%             AREA(i,10243:20484) = SurfStatReadData1(namesct_right(i,:));
%             isthere_ct(i) = 1;
%         catch
%             disp([namesct_left(i,:) ' not there'])
%         end
%         sum(isthere_ct)/length(isthere_ct)
%     end
%     
%     area200 = []
%     for i = 1:100
%         area200(i,:) = sum(AREA(:,find(parcels200==i+1)),2);
%     end
%     
%     for i = 1:100
%         area200(i+100,:) = sum(AREA(:,find(parcels200==i+1001)),2);
%     end

    % save data
    %csvwrite([DPATH 'CT_200_7.csv'], HCP200_CT);
    %csvwrite([DPATH 'SA_200_7.csv'], area200');
    %fid = fopen(['DPATH labels_200_7.csv','w')
    %fprintf(fid,'%s\n',names200{:,1})
    %fclose(fid)   
    
end
    
    % load subcortical volumes
    subcort_fs(:,1) = table(HCP.FS_L_AccumbensArea_Vol);
    subcort_fs(:,2) = table(HCP.FS_R_AccumbensArea_Vol);
    subcort_fs(:,3) = table(HCP.FS_L_Amygdala_Vol);
    subcort_fs(:,4) = table(HCP.FS_R_Amygdala_Vol);
    subcort_fs(:,5) = table(HCP.FS_L_Caudate_Vol);
    subcort_fs(:,6) = table(HCP.FS_R_Caudate_Vol);    
    subcort_fs(:,7) = table(HCP.FS_L_Hippo_Vol);
    subcort_fs(:,8) = table(HCP.FS_R_Hippo_Vol);
    subcort_fs(:,9) = table(HCP.FS_L_Pallidum_Vol);
    subcort_fs(:,10) = table(HCP.FS_R_Pallidum_Vol);    
    subcort_fs(:,11) = table(HCP.FS_L_Putamen_Vol);
    subcort_fs(:,12) = table(HCP.FS_R_Putamen_Vol);
    subcort_fs(:,13) = table(HCP.FS_L_ThalamusProper_Vol);
    subcort_fs(:,14) = table(HCP.FS_R_ThalamusProper_Vol);
    subcort_fs(:,15) = table(HCP.FS_L_VentDC_Vol);
    subcort_fs(:,16) = table(HCP.FS_R_VentDC_Vol);

    for i = 1:size(subcort_fs, 2)
        subcort_fs{:,i} = subcort_fs{:,i}./1000;
    end
    
    
    subcort_fs.Properties.VariableNames = {'accumb_l','accumb_r',...
        'amy_l','amy_r','caud_l','caud_r','hipp_l','hipp_r','pall_l',...
        'pall_r','put_l','put_r','thal_l','thal_r','ventDC_l','ventDC_r'};
    
end

%% find missing data and outliers
for out = 1
    
    % find participants with missing data
    studykeep = csvread([dir 'affect_cognition/data/studykeep_indices.csv']);
    [m1, m2] = find(isnan(list_of_vars(studykeep,:)));
    studykeep(m1) = [];
    %studykeep   = mintersect(find(mean(HCP200_CT,2)>0), find(HCP.FS_IntraCranial_Vol>0), find(~isnan(mean(list_of_vars,2))));
    
    outlier = ones(3,1206);
    for c = 1:1206
        r_a = corrcoef(area200(c,:),mean(area200));
        r_c = corrcoef(HCP200_CT(c,:),mean(HCP200_CT));
        if r_a(2) < 0.8
            outlier(1,c) = c;
            outlier(2,c) = r_a(2);
        elseif r_c(2) < 0.8
            outlier(1,c) = c;
            outlier(3,c) = r_c(2);
        else
            continue;
        end
    end
    %outlier in sa
    out_sa = find(outlier(2,:)<1)
    %outlier in ctx
    out_ct = find(outlier(3,:)<1)
    
    % exclude outliers
    studykeep(find(ismember(studykeep, [out_ct, out_sa]))) = [];

end

%% Phenotypic correlation, heritability and genetic correlation of cognition and affect 
for figure1 = 1
    
    % NIH positive affect composite  
    nihPA       = zeros(1,1206);
    nihPA(studykeep) = mean([HCP.LifeSatisf_Unadj(studykeep), HCP.MeanPurp_Unadj(studykeep), HCP.PosAffect_Unadj(studykeep)],2);
    
    % NIH negative affect composite  
    nihNA       = zeros(1,1206);
    nihNA(studykeep) = mean([HCP.AngAffect_Unadj(studykeep), HCP.AngHostil_Unadj(studykeep), HCP.Sadness_Unadj(studykeep),...
        HCP.FearAffect_Unadj(studykeep),HCP.PercStress_Unadj(studykeep)],2);
    
    % mean affect
    aff            = zeros(1,1206);
    aff(studykeep) = mean([nihPA(studykeep)', -nihNA(studykeep)'],2);
    
    
    
    list_of_vars = [HCP.CogTotalComp_Unadj, HCP.CogFluidComp_Unadj, HCP.CogCrystalComp_Unadj, nihPA', nihNA', aff'];    
    comptitlevar = {'Total Cognition', 'Fluid', 'Crystallized', 'Positive Affect', 'Negative Affect', 'Mean Affect'};
    
    
    T = table(mean(list_of_vars(studykeep,:))',std(list_of_vars(studykeep,:))',...
        min(list_of_vars(studykeep,:))',max(list_of_vars(studykeep,:))',...
         'VariableNames', {'Mean','SD','Min','Max'}, 'RowNames',comptitlevar')
    writetable(T,[RPATH 'Table1.csv']);         
    
   
    f = figure,
    subplot(1,6,1),
    hist(list_of_vars(studykeep,1))
    subplot(1,6,2),
    hist(list_of_vars(studykeep,2))
    subplot(1,6,3),
    hist(list_of_vars(studykeep,3))
    subplot(1,6,4),
    hist((list_of_vars(studykeep,4)))
    subplot(1,6,5),
    hist(list_of_vars(studykeep,5))
    subplot(1,6,6),
    hist(list_of_vars(studykeep,6))
    exportfigbo(f,[RPATH, 'F1A.personality.hist.png'],'png', 6)
    close(f)
    
    for i = 1:size(list_of_vars,2)
        vary   = list_of_vars(:,i);
        keep   = studykeep;
        vark   = zscore(vary(keep));
        agek   = zscore(HCP_r.Age_in_Yrs(keep));
        sexk   = cellstr(HCP.Gender(keep));
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark);
        slm  = SurfStatLinMod(zscore(list_of_vars(keep,:)),M);
        slm  = SurfStatT(slm, vark);
        t_var(i,:) = slm.t;
        beta_var(i,:) = slm.coef(end,:);
        pp_var(i,:) = 1- tcdf(slm.t,slm.df);
        pn_var(i,:) = 1-pp_var(i,:);
        
    end
    t_var(eye(6)==1)=1;
    beta_var(eye(6)==1)=1;
    pp_var(eye(6)==1)=1;
    pn_var(eye(6)==1)=1;
       
    
    fdr_bh([pp_var(1,2),pp_var(1,3),pp_var(1,4),pp_var(1,5),pp_var(1,6),...
        pp_var(2,3),pp_var(2,4),pp_var(2,5),pp_var(2,6),...
        pp_var(3,4),pp_var(3,5),pp_var(3,6),...
        pp_var(4,5),pp_var(4,6),...
        pp_var(5,6)],0.025) 
    ppn = squareform(ans);
    fdr_bh([pn_var(1,2),pn_var(1,3),pn_var(1,4),pn_var(1,5),pn_var(1,6),...
        pn_var(2,3),pn_var(2,4),pn_var(2,5),pn_var(2,6),...
        pn_var(3,4),pn_var(3,5),pn_var(3,6),...
        pn_var(4,5),pn_var(4,6),...
        pn_var(5,6)],0.025) 
   pnn = squareform(ans);
   
   ppp = ppn+pnn;
    
   ppp(eye(6)==1) =1;
   for explo = 1
        psignr = ppp == 1;
        %rleft = t_var.*psignr;
        rleft = beta_var.*psignr;
        mat = rleft;
        f = figure;
        imagesc(mat)
        colorbar;
        caxis([-1,1]);
        textStrings = num2str(mat(:), '%0.2f');       % Create strings from the matrix values
        textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
        textStrings = replace(textStrings,'0.00','-');
        [x, y] = meshgrid(1:6);  % Create x and y coordinates for the strings
        hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
            'HorizontalAlignment', 'center');
        midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
        textColors = repmat(abs(mat(:)) > 0.5, 1, 3);  % Choose white or black for the text color of the strings so they can be easily seen over the background color
        
        set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors
        set(gca, 'XTick', 1:(length(list_of_vars)));
        set(gca, 'XTickLabel', comptitlevar)
        xtickangle(-45)
        colormap(flipud(cbrewer('div','RdBu',99)));
        colorbar;
        set(gca, 'YTick', 1:(length(list_of_vars)));
        set(gca, 'YTickLabel', comptitlevar)
        caxis([-1,1]);
        title('Phenotypic correlation: betas, FDRq<0.05')
        
        print(f,[RPATH 'Corr_comp_explo_n'],'-dpng')
        close(f)
   end

   
   %plot total cognition and mean affect
   f=figure
   scatter(list_of_vars(:,1), list_of_vars(:,6),'.','k')
   l = lsline;
   l.Color = 'k'; 
   l.LineWidth = 2;
   xlabel('Total Cognition')
   ylabel('Mean Affect')
   print(f,[RPATH 'Scatter_cog_affect'],'-dpng')
   close(f)   
   
   
   
   %heritability:
    F1b = readtable([dir 'affect_cognition/solar_subc/cog_affect_heritability.csv']);

    %rearrange order of traits, so it matches with order above
    ord = [6 3 2 5 4 1];
    F1b = F1b(ord,:);   
    
    f = figure,
    bar(F1b.H2r)
    set(gca, 'XTickLabel', comptitlevar)
    xtickangle(-45)
    print(f,[RPATH 'Heri_1B'],'-dpng')
    close(f)
    
    
    %genetic correlation
    F1c = readtable([dir 'affect_cognition/solar_subc/cog_aff_gencorr.csv']);      

    gc_bv = eye(6);  
    gc_bv(gc_bv==0)=F1c.rG;    
    gc_bv_p = eye(6); 
    gc_bv_p(gc_bv_p==0)=F1c.p;
    
    %rearrange order of traits, so it matches with order above
    ord = [6 3 2 5 4 1];
    gc_bv = gc_bv(ord,ord);
    gc_bv_p = gc_bv_p(ord,ord);
    
    %FDR correction   
    pp_var = gc_bv_p;
       fdr_bh([pp_var(1,2),pp_var(1,3),pp_var(1,4),pp_var(1,5),pp_var(1,6),...
        pp_var(2,3),pp_var(2,4),pp_var(2,5),pp_var(2,6),...
        pp_var(3,4),pp_var(3,5),pp_var(3,6),...
        pp_var(4,5),pp_var(4,6),...
        pp_var(5,6)],0.025) 
    ppn = squareform(ans);
   
    
    for explo = 1
        psignr = ppn;
        rleft = gc_bv.*psignr;
        mat = rleft;
        f = figure;
        imagesc(mat)
        colorbar;
        caxis([-1,1]);
        textStrings = num2str(mat(:), '%0.2f');       % Create strings from the matrix values
        textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
        textStrings = replace(textStrings,'0.00','-');
        [x, y] = meshgrid(1:6);  % Create x and y coordinates for the strings
        hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
            'HorizontalAlignment', 'center');
        midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
        textColors = repmat(mat(:) > 0.5, 1, 3);  % Choose white or black for the text color of the strings so they can be easily seen over the background color
        
        set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors
        set(gca, 'XTick', 1:(length(list_of_vars)));
        set(gca, 'XTickLabel', comptitlevar)
        xtickangle(-45)
        colormap(flipud(cbrewer('div','RdBu',99)));
        colorbar;
        set(gca, 'YTick', 1:(length(list_of_vars)));
        set(gca, 'YTickLabel', comptitlevar)
        caxis([-1,1]);
        title('Genetic correlation: rhoG, FDRq<0.05')
        
        print(f,[RPATH 'F1C'],'-dpng')
        close(f)
   end

    
   
end

%% Phenotypic correlation of cognition and brain structure
%cognition: total cognition, fluid, crystallized

% cortical thickness
for ctx = 1
    % total
    for k = 1
        total  = HCP.CogTotalComp_Unadj;        
        keep   = studykeep;
        vark   = (total(keep));
        agek   = HCP_r.Age_in_Yrs(keep);
        sexk   = cellstr(HCP.Gender(keep));
        icv    = HCP.FS_IntraCranial_Vol(keep);
        gb     = mean(HCP200_CT(keep,:),2);
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(gb) + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark)
        slm  = SurfStatLinMod(HCP200_CT(keep,:),M)
        slm  = SurfStatT(slm, vark)
        pp   = 1 - tcdf(slm.t, slm.df)
        pn   = 1 - tcdf(-slm.t, slm.df)
        p= zeros(size(pp));
        p = pp<pn;
        p_all= zeros(size(p));
        p_all(p==1) = pp(p==1);
        p_all(p==0) = pn(p==0);        
        
        [h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pp,0.025);
        [h2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pn,0.025);
        
        h = h1+h2;

        var_ctx_t(1,:) = slm.t;
        var_ctx_p(1,:) = p_all;
        var_ctx_FDR(1,:) = h;
        
        if max(h) == 1
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure;
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(cbrewer('div','RdBu',11)));
            SurfStatColLim([-3 3])
            exportfigbo(f,[RPATH, 'F2pheno.total.ctx.png'],'png', 10)
            close(f)
            
        end
        
    end
    % fluid
    for k = 1
        fluid  = HCP.CogFluidComp_Unadj;
        crys   = HCP.CogCrystalComp_Unadj;
        keep   = studykeep;
        vark   = (fluid(keep));
        agek   = HCP_r.Age_in_Yrs(keep);
        sexk   = cellstr(HCP.Gender(keep));
        icv    = HCP.FS_IntraCranial_Vol(keep);
        gb     = mean(HCP200_CT(keep,:),2);
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(gb) + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark)
        slm  = SurfStatLinMod(HCP200_CT(keep,:),M)
        slm  = SurfStatT(slm, vark)
        pp   = 1 - tcdf(slm.t, slm.df)
        pn   = 1 - tcdf(-slm.t, slm.df)
        p= zeros(size(pp));
        p = pp<pn;
        p_all= zeros(size(p));
        p_all(p==1) = pp(p==1);
        p_all(p==0) = pn(p==0);        
        
        [h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pp,0.025);
        [h2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pn,0.025);
        
        h = h1+h2;
        var_ctx_t(2,:) = slm.t;
        var_ctx_p(2,:) = p_all;
        var_ctx_FDR(2,:) = h;
        
        if max(h) == 1
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(cbrewer('div','RdBu',11)));            
            SurfStatColLim([-3 3])
            exportfigbo(f,[RPATH, 'F2pheno.fluid.ctx.png'],'png', 10)
            close(f)
            
        end
        
    end
    % crystallized
    for k = 1
        fluid  = HCP.CogFluidComp_Unadj; 
        crys   = HCP.CogCrystalComp_Unadj;
        keep   = studykeep;
        vark   = (crys(keep));
        agek   = HCP_r.Age_in_Yrs(keep);
        sexk   = cellstr(HCP.Gender(keep));
        icv    = HCP.FS_IntraCranial_Vol(keep);
        gb     = mean(HCP200_CT(keep,:),2);
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(gb) + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark)
        slm  = SurfStatLinMod(HCP200_CT(keep,:),M)
        slm  = SurfStatT(slm, vark)
        pp   = 1 - tcdf(slm.t, slm.df)
        pn   = 1 - tcdf(-slm.t, slm.df)
        p= zeros(size(pp));
        p = pp<pn;
        p_all= zeros(size(p));
        p_all(p==1) = pp(p==1);
        p_all(p==0) = pn(p==0);        
        
        [h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pp,0.025);
        [h2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pn,0.025);
        
        h = h1+h2;
        
        var_ctx_t(3,:) = slm.t;
        var_ctx_p(3,:) = p_all;
        var_ctx_FDR(3,:) = h;
   
        if max(h) == 1 
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure;
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(cbrewer('div','RdBu',11)));            
            SurfStatColLim([-3 3])
            exportfigbo(f,[RPATH, 'F2pheno.crys.ctx.png'],'png', 10)
            close(f)
            
        end
        
    end
end

% surface area 
for area = 1
    % total
    for k = 1
        total  = HCP.CogTotalComp_Unadj;        
        keep   = studykeep;
        vark   = (total(keep));
        agek   = HCP_r.Age_in_Yrs(keep);
        sexk   = cellstr(HCP.Gender(keep));
        icv    = HCP.FS_IntraCranial_Vol(keep); 
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(icv) + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark)
        slm  = SurfStatLinMod(area200(:,keep)',M)
        slm  = SurfStatT(slm, vark)
        pp   = 1 - tcdf(slm.t, slm.df)
        pn   = 1 - tcdf(-slm.t, slm.df)
        p= zeros(size(pp));
        p = pp<pn;
        p_all= zeros(size(p));
        p_all(p==1) = pp(p==1);
        p_all(p==0) = pn(p==0);        
        
        [h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pp,0.025);
        [h2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pn,0.025);
        
        h = h1+h2;
        
        var_area_t(1,:) = slm.t;
        var_area_p(1,:) = p_all;
        var_area_FDR(1,:) = h;
       
        if max(h) == 1
          
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure;
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(cbrewer('div','RdBu',11)));
            SurfStatColLim([-3 3])
            exportfigbo(f,[RPATH, 'F2pheno.total.area.png'],'png', 10)
            close(f)
            
        end
        
    end
    % crystallized
    for k = 1
        fluid  = HCP.CogFluidComp_Unadj;
        crys   = HCP.CogCrystalComp_Unadj;
        keep   = studykeep;
        vark   = (crys(keep));
        agek   = HCP_r.Age_in_Yrs(keep);
        sexk   = cellstr(HCP.Gender(keep));
        icv    = HCP.FS_IntraCranial_Vol(keep);        
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(icv) + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark) 
        slm  = SurfStatLinMod(area200(:,keep)',M)
        slm  = SurfStatT(slm, vark)
        pp   = 1 - tcdf(slm.t, slm.df)
        pn   = 1 - tcdf(-slm.t, slm.df)
        p= zeros(size(pp));
        p = pp<pn;
        p_all= zeros(size(p));
        p_all(p==1) = pp(p==1);
        p_all(p==0) = pn(p==0);        
        
        [h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pp,0.025);
        [h2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pn,0.025);
        
        h = h1+h2;
        
        var_area_t(3,:) = slm.t;
        var_area_p(3,:) = p_all;
        var_area_FDR(3,:) = h;
       
        if max(h) == 1
          
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(cbrewer('div','RdBu',11)));
            SurfStatColLim([-3 3])
            exportfigbo(f,[RPATH, 'F2pheno.crys.area.png'],'png', 10)
            close(f)
            
        end
        
    end
    % fluid 
    for k = 1
        fluid  = HCP.CogFluidComp_Unadj; 
        crys   = HCP.CogCrystalComp_Unadj;
        keep   = studykeep;
        vark   = (fluid(keep));
        agek   = HCP_r.Age_in_Yrs(keep);
        sexk   = cellstr(HCP.Gender(keep));
        icv    = HCP.FS_IntraCranial_Vol(keep); 
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(icv) + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark) 
        slm  = SurfStatLinMod(area200(:,keep)',M)
        slm  = SurfStatT(slm, vark)
        pp   = 1 - tcdf(slm.t, slm.df)
        pn   = 1 - tcdf(-slm.t, slm.df)
        p= zeros(size(pp));
        p = pp<pn;
        p_all= zeros(size(p));
        p_all(p==1) = pp(p==1);
        p_all(p==0) = pn(p==0);        
        
        [h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pp,0.025);
        [h2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pn,0.025);
        
        h = h1+h2;
        
        var_area_t(2,:) = slm.t;
        var_area_p(2,:) = p_all;
        var_area_FDR(2,:) = h;
     
        if max(h) == 1
           
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(cbrewer('div','RdBu',11)));
            SurfStatColLim([-3 3])
            exportfigbo(f,[RPATH, 'F2pheno.fluid.area.png'],'png', 10)
            close(f)
            
        end
        
    end
end

% subcortical voumes
for subcort = 1
    % total
    for k = 1
        total  = HCP.CogTotalComp_Unadj;        
        keep   = studykeep;
        vark   = (total(keep));
        agek   = HCP_r.Age_in_Yrs(keep);
        sexk   = cellstr(HCP.Gender(keep));
        icv    = HCP.FS_IntraCranial_Vol(keep); 
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(icv) + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark)
        slm  = SurfStatLinMod(subcort_fs{keep,:},M)
        slm  = SurfStatT(slm, vark)
        pp   = 1 - tcdf(slm.t, slm.df)
        pn   = 1 - tcdf(-slm.t, slm.df)
        p= zeros(size(pp));
        p = pp<pn;
        p_all= zeros(size(p));
        p_all(p==1) = pp(p==1);
        p_all(p==0) = pn(p==0);        
        
        [h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pp,0.025);
        [h2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pn,0.025);
        
        h = h1+h2;
        
        var_sub_t(1,:) = slm.t;
        var_sub_p(1,:) = p_all;
        var_sub_FDR(1,:) = h;
       
        if max(h) == 1                    
            
            f = figure                      
            bar(slm.t)
            set(gca, 'XTick', 1:(length(subcort_lina.Properties.VariableNames)));
            set(gca, 'XTickLabel', subcort_lina.Properties.VariableNames)
            xtickangle(-45)
            exportfigbo(f,[RPATH, 'F2pheno.total.sub.png'],'png', 10)
            close(f)       
            
        end
        
    end
    % crystallized
    for k = 1
        fluid  = HCP.CogFluidComp_Unadj;
        crys   = HCP.CogCrystalComp_Unadj;
        keep   = studykeep;
        vark   = (crys(keep));
        agek   = HCP_r.Age_in_Yrs(keep);
        sexk   = cellstr(HCP.Gender(keep));
        icv    = HCP.FS_IntraCranial_Vol(keep);        
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(icv) + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark)
        slm  = SurfStatLinMod(subcort_fs{keep,:},M)
        slm  = SurfStatT(slm, vark)
        pp   = 1 - tcdf(slm.t, slm.df)
        pn   = 1 - tcdf(-slm.t, slm.df)
        p= zeros(size(pp));
        p = pp<pn;
        p_all= zeros(size(p));
        p_all(p==1) = pp(p==1);
        p_all(p==0) = pn(p==0);        
        
        [h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pp,0.025);
        [h2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pn,0.025);
        
        h = h1+h2;
        
        var_sub_t(3,:) = slm.t;
        var_sub_p(3,:) = p_all;
        var_sub_FDR(3,:) = h;
       
        if max(h) == 1
           
            p_val = h;
            
            f = figure,
            bar(slm.t)
            set(gca, 'XTick', 1:(length(subcort_lina.Properties.VariableNames)));
            set(gca, 'XTickLabel', subcort_lina.Properties.VariableNames)
            xtickangle(-45)
            exportfigbo(f,[RPATH, 'F2pheno.crys.sub.png'],'png', 10)
            close(f)        
            
        end
    end
    % fluid
    for k = 1
        fluid  = HCP.CogFluidComp_Unadj; 
        crys   = HCP.CogCrystalComp_Unadj;
        keep   = studykeep;
        vark   = (fluid(keep));
        agek   = HCP_r.Age_in_Yrs(keep);
        sexk   = cellstr(HCP.Gender(keep));
        icv    = HCP.FS_IntraCranial_Vol(keep); 
        sex_num = grp2idx(sexk);

        
        M    = 1 + term(icv) + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark)
        slm  = SurfStatLinMod(subcort_fs{keep,:},M)
        slm  = SurfStatT(slm, vark)
        pp   = 1 - tcdf(slm.t, slm.df)
        pn   = 1 - tcdf(-slm.t, slm.df)
        p= zeros(size(pp));
        p = pp<pn;
        p_all= zeros(size(p));
        p_all(p==1) = pp(p==1);
        p_all(p==0) = pn(p==0);        
        
        [h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pp,0.025);
        [h2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pn,0.025);
        
        h = h1+h2;
        
        var_sub_t(2,:) = slm.t;
        var_sub_p(2,:) = p_all;
        var_sub_FDR(2,:) = h;
       
        if max(h) == 1
           
            p_val = h;
            
            f = figure,
            bar(slm.t)
            set(gca, 'XTick', 1:(length(subcort_lina.Properties.VariableNames)));
            set(gca, 'XTickLabel', subcort_lina.Properties.VariableNames)
            xtickangle(-45)
            exportfigbo(f,[RPATH, 'F2pheno.fluid.sub.png'],'png', 10)
            close(f)                       
            
        end
        
    end
end

%% Phenotypic correlation of affect and brain structure
% affect: mean affect, positive affect, negative affect

% cortical thickness
for ctx = 1
    for j = 4:6
        vary  =  list_of_vars(:,j);
        keep   = studykeep;
        vark   = (vary(keep));
        agek   = HCP_r.Age_in_Yrs(keep);
        sexk   = cellstr(HCP.Gender(keep));
        icv    = HCP.FS_IntraCranial_Vol(keep);
        gb     = mean(HCP200_CT(keep,:),2);
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(gb) + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark)
        slm  = SurfStatLinMod(HCP200_CT(keep,:),M)
        slm  = SurfStatT(slm, vark)
        pp   = 1 - tcdf(slm.t, slm.df)
        pn   = 1 - tcdf(-slm.t, slm.df)
        p= zeros(size(pp));
        p = pp<pn;
        p_all= zeros(size(p));
        p_all(p==1) = pp(p==1);
        p_all(p==0) = pn(p==0);
        
        [h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pp,0.025);
        [h2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pn,0.025);
        
        h = h1+h2;
        
        var_ctx_FDR(j,:) = h;
        var_ctx_p(j,:) = p_all;
        var_ctx_t(j,:) = slm.t; 
        if max(h) == 1
          
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(cbrewer('div','RdBu',11)));
            SurfStatColLim([-3 3])
            exportfigbo(f,[RPATH, 'F3pheno.', comptitlevar{j},'.ctx.png'],'png', 10)
            close(f)
            
        end
        
    end

end

% surface area
for area = 1
    for j = 4:6
        vary  =  list_of_vars(:,j);
        keep   = studykeep;
        vark   = (vary(keep));
        agek   = HCP_r.Age_in_Yrs(keep);
        sexk   = cellstr(HCP.Gender(keep));
        icv    = HCP.FS_IntraCranial_Vol(keep);
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(icv) + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark)
        slm  = SurfStatLinMod(area200(:,keep)',M)
        slm  = SurfStatT(slm, vark)
        pp   = 1 - tcdf(slm.t, slm.df)
        pn   = 1 - tcdf(-slm.t, slm.df)
        p= zeros(size(pp));
        p = pp<pn;
        p_all= zeros(size(p));
        p_all(p==1) = pp(p==1);
        p_all(p==0) = pn(p==0);        
        
        [h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pp,0.025);
        [h2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pn,0.025);
        
        h = h1+h2;
        
        var_area_FDR(j,:) = h;
        var_area_p(j,:) = p_all;
        var_area_t(j,:) = slm.t; 
        if max(h) == 1
           
            p_val = h;
            
            heri_ct = zeros(1,20484);
            for i = 1:100
                heri_ct(:,find(parcels200==i+1)) = (slm.t(i)).*(p_val(i));
            end
            for i = 1:100
                heri_ct(:,find(parcels200==i+1001)) = (slm.t(i+100)).*(p_val(i+100));
            end
            
            ctx2 = cbrewer('div','Spectral',10);
            f = figure,
            BoSurfStatViewData(heri_ct,SN,'')
            colormap(flipud(cbrewer('div','RdBu',11)));            
            SurfStatColLim([-3 3])
            exportfigbo(f,[RPATH, 'F3pheno.', comptitlevar{j},'.area.png'],'png', 10)
            close(f)
            
        end
        
    end
end

% subcortical
for subcort = 1
    for j = 4:6
        vary  =  list_of_vars(:,j);
        keep   = studykeep;
        vark   = (vary(keep));
        agek   = HCP_r.Age_in_Yrs(keep);
        sexk   = cellstr(HCP.Gender(keep));
        icv    = HCP.FS_IntraCranial_Vol(keep);
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(icv) + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark)
        slm  = SurfStatLinMod(subcort_fs{keep,:},M)
        slm  = SurfStatT(slm, vark)
        pp   = 1 - tcdf(slm.t, slm.df)
        pn   = 1 - tcdf(-slm.t, slm.df)
        p= zeros(size(pp));
        p = pp<pn;
        p_all= zeros(size(p));
        p_all(p==1) = pp(p==1);
        p_all(p==0) = pn(p==0);        
        
        [h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pp,0.025);
        [h2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pn,0.025);
        
        h = h1+h2;
        
        var_sub_t(j,:) = slm.t;
        var_sub_p(j,:) = p_all;
        var_sub_FDR(j,:) = h;
       
        if max(h) == 1
            
            f = figure,
            bar(slm.t)
            set(gca, 'XTick', 1:(length(subcort_lina.Properties.VariableNames)));
            set(gca, 'XTickLabel', subcort_lina.Properties.VariableNames)
            xtickangle(-45)
            exportfigbo(f,[RPATH, 'F3pheno.', comptitlevar{j},'.sub.png'],'png', 10)
            close(f) 
            
        end
        
    end
    
end

%% Phenotypic correlation tables
T = table(pheno.id);
writetable(T, [RPATH 'included_ids.csv'], 'WriteRowNames', false)

try
    
for ctx_post_hoc = 1
    %% all in total cognition
    i =1
     rn = names200(var_ctx_FDR(i,:)==1)
    T= table(var_ctx_t(1,find(var_ctx_FDR(i,:)==1))',var_ctx_p(1,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(2,find(var_ctx_FDR(i,:)==1))',var_ctx_p(2,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(3,find(var_ctx_FDR(i,:)==1))',var_ctx_p(3,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(4,find(var_ctx_FDR(i,:)==1))',var_ctx_p(4,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(5,find(var_ctx_FDR(i,:)==1))',var_ctx_p(5,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(6,find(var_ctx_FDR(i,:)==1))',var_ctx_p(6,find(var_ctx_FDR(i,:)==1))',...
        'RowNames',rn);
     T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar
    
     writetable(T, [RPATH 'total_ctx.csv'], 'WriteRowNames', true)
     
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'total_ctx'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);     
     
    
    %% all in fluid
     i =2
     rn = names200(var_ctx_FDR(i,:)==1)
    T= table(var_ctx_t(1,find(var_ctx_FDR(i,:)==1))',var_ctx_p(1,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(2,find(var_ctx_FDR(i,:)==1))',var_ctx_p(2,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(3,find(var_ctx_FDR(i,:)==1))',var_ctx_p(3,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(4,find(var_ctx_FDR(i,:)==1))',var_ctx_p(4,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(5,find(var_ctx_FDR(i,:)==1))',var_ctx_p(5,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(6,find(var_ctx_FDR(i,:)==1))',var_ctx_p(6,find(var_ctx_FDR(i,:)==1))',...
        'RowNames',rn)
   T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'fluid_ctx.csv'], 'WriteRowNames', true)
     
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'fluid_ctx'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);     
     
     
    %% all in crystallized
    i =3
     rn = names200(var_ctx_FDR(i,:)==1)
    T= table(var_ctx_t(1,find(var_ctx_FDR(i,:)==1))',var_ctx_p(1,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(2,find(var_ctx_FDR(i,:)==1))',var_ctx_p(2,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(3,find(var_ctx_FDR(i,:)==1))',var_ctx_p(3,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(4,find(var_ctx_FDR(i,:)==1))',var_ctx_p(4,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(5,find(var_ctx_FDR(i,:)==1))',var_ctx_p(5,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(6,find(var_ctx_FDR(i,:)==1))',var_ctx_p(6,find(var_ctx_FDR(i,:)==1))',...
        'RowNames',rn)
   T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'crystal_ctx.csv'], 'WriteRowNames', true)
     
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'crystal_ctx'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);     
          
     
    %% all in positive affect
    i =4
     rn = names200(var_ctx_FDR(i,:)==1)
    T= table(var_ctx_t(1,find(var_ctx_FDR(i,:)==1))',var_ctx_p(1,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(2,find(var_ctx_FDR(i,:)==1))',var_ctx_p(2,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(3,find(var_ctx_FDR(i,:)==1))',var_ctx_p(3,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(4,find(var_ctx_FDR(i,:)==1))',var_ctx_p(4,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(5,find(var_ctx_FDR(i,:)==1))',var_ctx_p(5,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(6,find(var_ctx_FDR(i,:)==1))',var_ctx_p(6,find(var_ctx_FDR(i,:)==1))',...
        'RowNames',rn)
   T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'nihPA_ctx.csv'], 'WriteRowNames', true)

    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'nihPA_ctx'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);     

    %% all in negative affect 
    i =5 %(!!!) 
     rn = names200(var_ctx_FDR(i,:)==1)
    T = table(var_ctx_t(1,find(var_ctx_FDR(i,:)==1))',var_ctx_p(1,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(2,find(var_ctx_FDR(i,:)==1))',var_ctx_p(2,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(3,find(var_ctx_FDR(i,:)==1))',var_ctx_p(3,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(4,find(var_ctx_FDR(i,:)==1))',var_ctx_p(4,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(5,find(var_ctx_FDR(i,:)==1))',var_ctx_p(5,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(6,find(var_ctx_FDR(i,:)==1))',var_ctx_p(6,find(var_ctx_FDR(i,:)==1))',...
        'RowNames',rn)
    T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'nihNA_ctx.csv'], 'WriteRowNames', true)
 
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'nihNA_ctx'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);      
     
     %% all in affect
     i =6 %(!!!) 
     rn = names200(var_ctx_FDR(i,:)==1)
    T = table(var_ctx_t(1,find(var_ctx_FDR(i,:)==1))',var_ctx_p(1,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(2,find(var_ctx_FDR(i,:)==1))',var_ctx_p(2,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(3,find(var_ctx_FDR(i,:)==1))',var_ctx_p(3,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(4,find(var_ctx_FDR(i,:)==1))',var_ctx_p(4,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(5,find(var_ctx_FDR(i,:)==1))',var_ctx_p(5,find(var_ctx_FDR(i,:)==1))',...
        var_ctx_t(6,find(var_ctx_FDR(i,:)==1))',var_ctx_p(6,find(var_ctx_FDR(i,:)==1))',...
        'RowNames',rn)
    T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'aff_ctx.csv'], 'WriteRowNames', true)
     
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'aff_ctx'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);      

 
end
catch
    disp('There was an empty table somewhere');
end

try
for area_post_hoc = 1
    %% all in total
     rn = names200(var_area_FDR(1,:)==1)
     T = table(var_area_t(1,find(var_area_FDR(1,:)==1))',var_area_p(1,find(var_area_FDR(1,:)==1))',...
        var_area_t(2,find(var_area_FDR(1,:)==1))',var_area_p(2,find(var_area_FDR(1,:)==1))',...
        var_area_t(3,find(var_area_FDR(1,:)==1))',var_area_p(3,find(var_area_FDR(1,:)==1))',...
        var_area_t(4,find(var_area_FDR(1,:)==1))',var_area_p(4,find(var_area_FDR(1,:)==1))',...
        var_area_t(5,find(var_area_FDR(1,:)==1))',var_area_p(5,find(var_area_FDR(1,:)==1))',...
        var_area_t(6,find(var_area_FDR(1,:)==1))',var_area_p(6,find(var_area_FDR(1,:)==1))',...
        'RowNames',rn)
    T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
   
      writetable(T, [RPATH 'total_area.csv'], 'WriteRowNames', true)
      
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'total_area'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);       
 
    %% all in fluid IQ
    i =2
       rn = names200(var_area_FDR(2,:)==1)
    T = table(var_area_t(1,find(var_area_FDR(i,:)==1))',var_area_p(1,find(var_area_FDR(i,:)==1))',...
        var_area_t(2,find(var_area_FDR(i,:)==1))',var_area_p(2,find(var_area_FDR(i,:)==1))',...
        var_area_t(3,find(var_area_FDR(i,:)==1))',var_area_p(3,find(var_area_FDR(i,:)==1))',...
        var_area_t(4,find(var_area_FDR(i,:)==1))',var_area_p(4,find(var_area_FDR(i,:)==1))',...
        var_area_t(5,find(var_area_FDR(i,:)==1))',var_area_p(5,find(var_area_FDR(i,:)==1))',...
        var_area_t(6,find(var_area_FDR(i,:)==1))',var_area_p(6,find(var_area_FDR(i,:)==1))',...
        'RowNames',rn)
     T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'fluid_area.csv'], 'WriteRowNames', true)
     
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'fluid_area'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);      
 
    %% all in crystallized
    i =3
       rn = names200(var_area_FDR(3,:)==1)
    T = table(var_area_t(1,find(var_area_FDR(i,:)==1))',var_area_p(1,find(var_area_FDR(i,:)==1))',...
        var_area_t(2,find(var_area_FDR(i,:)==1))',var_area_p(2,find(var_area_FDR(i,:)==1))',...
        var_area_t(3,find(var_area_FDR(i,:)==1))',var_area_p(3,find(var_area_FDR(i,:)==1))',...
        var_area_t(4,find(var_area_FDR(i,:)==1))',var_area_p(4,find(var_area_FDR(i,:)==1))',...
        var_area_t(5,find(var_area_FDR(i,:)==1))',var_area_p(5,find(var_area_FDR(i,:)==1))',...
        var_area_t(6,find(var_area_FDR(i,:)==1))',var_area_p(6,find(var_area_FDR(i,:)==1))',...
        'RowNames',rn)
    T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
    writetable(T, [RPATH 'crystal_area.csv'], 'WriteRowNames', true)
    
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'crystal_area'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);     

    %% all in positive affect
    i =4
     rn = names200(var_area_FDR(i,:)==1)
    T= table(var_ctx_t(1,find(var_area_FDR(i,:)==1))',var_ctx_p(1,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(2,find(var_area_FDR(i,:)==1))',var_ctx_p(2,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(3,find(var_area_FDR(i,:)==1))',var_ctx_p(3,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(4,find(var_area_FDR(i,:)==1))',var_ctx_p(4,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(5,find(var_area_FDR(i,:)==1))',var_ctx_p(5,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(6,find(var_area_FDR(i,:)==1))',var_ctx_p(6,find(var_area_FDR(i,:)==1))',...
        'RowNames',rn)
   T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'nihPA_area.csv'], 'WriteRowNames', true)

    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'nihPA_area'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);     

    %% all in negative affect
    i =5
     rn = names200(var_area_FDR(i,:)==1)
    T= table(var_ctx_t(1,find(var_area_FDR(i,:)==1))',var_ctx_p(1,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(2,find(var_area_FDR(i,:)==1))',var_ctx_p(2,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(3,find(var_area_FDR(i,:)==1))',var_ctx_p(3,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(4,find(var_area_FDR(i,:)==1))',var_ctx_p(4,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(5,find(var_area_FDR(i,:)==1))',var_ctx_p(5,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(6,find(var_area_FDR(i,:)==1))',var_ctx_p(6,find(var_area_FDR(i,:)==1))',...
        'RowNames',rn)
   T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'nihNA_area.csv'], 'WriteRowNames', true)

    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'nihNA_area'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d); 

    %% all in affect
    i =6
     rn = names200(var_area_FDR(i,:)==1)
    T= table(var_ctx_t(1,find(var_area_FDR(i,:)==1))',var_ctx_p(1,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(2,find(var_area_FDR(i,:)==1))',var_ctx_p(2,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(3,find(var_area_FDR(i,:)==1))',var_ctx_p(3,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(4,find(var_area_FDR(i,:)==1))',var_ctx_p(4,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(5,find(var_area_FDR(i,:)==1))',var_ctx_p(5,find(var_area_FDR(i,:)==1))',...
        var_ctx_t(6,find(var_area_FDR(i,:)==1))',var_ctx_p(6,find(var_area_FDR(i,:)==1))',...
        'RowNames',rn)
   T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'aff_area.csv'], 'WriteRowNames', true)

    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'aff_area'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d); 

end
catch
    disp('There was an empty table somewhere');
end

try
for subc_post_hoc = 1
    %% all in total
    i =1
     rn = subcort_lina.Properties.VariableNames(var_sub_FDR(i,:)==1);
    T= table(var_sub_t(1,find(var_sub_FDR(i,:)==1))',var_sub_p(1,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(2,find(var_sub_FDR(i,:)==1))',var_sub_p(2,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(3,find(var_sub_FDR(i,:)==1))',var_sub_p(3,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(4,find(var_sub_FDR(i,:)==1))',var_sub_p(4,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(5,find(var_sub_FDR(i,:)==1))',var_sub_p(5,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(6,find(var_sub_FDR(i,:)==1))',var_sub_p(6,find(var_sub_FDR(i,:)==1))',...
        'RowNames',rn)
     T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'total_subc.csv'], 'WriteRowNames', true)
     
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'total_subc'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d); 

     
    
    %% all in fluid
     i =2
     rn = subcort_lina.Properties.VariableNames(var_sub_FDR(i,:)==1);
    T= table(var_sub_t(1,find(var_sub_FDR(i,:)==1))',var_sub_p(1,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(2,find(var_sub_FDR(i,:)==1))',var_sub_p(2,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(3,find(var_sub_FDR(i,:)==1))',var_sub_p(3,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(4,find(var_sub_FDR(i,:)==1))',var_sub_p(4,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(5,find(var_sub_FDR(i,:)==1))',var_sub_p(5,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(6,find(var_sub_FDR(i,:)==1))',var_sub_p(6,find(var_sub_FDR(i,:)==1))',...
        'RowNames',rn)
   T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'fluid_subc.csv'], 'WriteRowNames', true)
     
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'fluid_subc'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);      
     
    %% all in crystallized
    i =3
     rn = subcort_lina.Properties.VariableNames(var_sub_FDR(i,:)==1);
    T= table(var_sub_t(1,find(var_sub_FDR(i,:)==1))',var_sub_p(1,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(2,find(var_sub_FDR(i,:)==1))',var_sub_p(2,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(3,find(var_sub_FDR(i,:)==1))',var_sub_p(3,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(4,find(var_sub_FDR(i,:)==1))',var_sub_p(4,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(5,find(var_sub_FDR(i,:)==1))',var_sub_p(5,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(6,find(var_sub_FDR(i,:)==1))',var_sub_p(6,find(var_sub_FDR(i,:)==1))',...
        'RowNames',rn)
   T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'crystal_subc.csv'], 'WriteRowNames', true)
     
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'crystal_subc'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);       
     
    %% all in positive affect
    i =4
     rn = subcort_lina.Properties.VariableNames(var_sub_FDR(i,:)==1);
    T= table(var_sub_t(1,find(var_sub_FDR(i,:)==1))',var_sub_p(1,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(2,find(var_sub_FDR(i,:)==1))',var_sub_p(2,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(3,find(var_sub_FDR(i,:)==1))',var_sub_p(3,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(4,find(var_sub_FDR(i,:)==1))',var_sub_p(4,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(5,find(var_sub_FDR(i,:)==1))',var_sub_p(5,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(6,find(var_sub_FDR(i,:)==1))',var_sub_p(6,find(var_sub_FDR(i,:)==1))',...
        'RowNames',rn)
   T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'nihPA_subc.csv'], 'WriteRowNames', true)
     
     
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'nihPA_subc'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);     
     
    %% all in negative affect 
    i = 5 %(!!!) 
     rn = subcort_lina.Properties.VariableNames(var_sub_FDR(i,:)==1);
    T= table(var_sub_t(1,find(var_sub_FDR(i,:)==1))',var_sub_p(1,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(2,find(var_sub_FDR(i,:)==1))',var_sub_p(2,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(3,find(var_sub_FDR(i,:)==1))',var_sub_p(3,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(4,find(var_sub_FDR(i,:)==1))',var_sub_p(4,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(5,find(var_sub_FDR(i,:)==1))',var_sub_p(5,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(6,find(var_sub_FDR(i,:)==1))',var_sub_p(6,find(var_sub_FDR(i,:)==1))',...
        'RowNames',rn)
    T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'nihNA_subc.csv'], 'WriteRowNames', true)
     
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'nihNA_subc'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);  

     %% all in affect
     i =6 %(!!!) 
     rn = subcort_lina.Properties.VariableNames(var_sub_FDR(i,:)==1);
    T= table(var_sub_t(1,find(var_sub_FDR(i,:)==1))',var_sub_p(1,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(2,find(var_sub_FDR(i,:)==1))',var_sub_p(2,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(3,find(var_sub_FDR(i,:)==1))',var_sub_p(3,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(4,find(var_sub_FDR(i,:)==1))',var_sub_p(4,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(5,find(var_sub_FDR(i,:)==1))',var_sub_p(5,find(var_sub_FDR(i,:)==1))',...
        var_sub_t(6,find(var_sub_FDR(i,:)==1))',var_sub_p(6,find(var_sub_FDR(i,:)==1))',...
        'RowNames',rn)
    T.Properties.VariableNames([1 3 5 7 9 11]) = comptitlevar;
    
     writetable(T, [RPATH 'aff_subc.csv'], 'WriteRowNames', true)

     
    mltableObj = MATLABTable(T);
    th = mltableObj.Header;
    thentry11 = entry(th,1,1);
    thentry11.Children(1).Children(1).Content = 'Variables';
    mltableObj.Style = t.Style;    
    mltableObj.RowNamesRule = true;
    mltableObj.Header.Style = [mltableObj.Header.Style {Italic(true)}];
    mltableObj.TableEntriesHAlign = 'right';
    mltableObj.Header.TableEntriesHAlign = 'center';
    
    d = Document([RPATH 'aff_subc'],'docx');
    append(d,mltableObj);       
    close(d);
%     rptview(d);  

end
catch
    disp('There was an empty table somewhere');
end

%% Brain structure heritability and genetic correlation analyes
% heritability
for load_her_CTX = 1
    ctx_her = readtable([dir '/affect_cognition/solar_CT/ct_heritability.csv']);
    ctx_her.name = cellfun(@(x) x(1:end-6), ctx_her.name, 'un', 0);
    
    idx = zeros(size(names200));
    for n = 1:length(names200)
        idx(n,1) = find(endsWith(ctx_her.name, names200{n}));
    end
    ctx_her = ctx_her(idx,:);
    writetable(ctx_her,[RPATH 'ctx_heritability.csv'], 'WriteRowNames', true);
    
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(ctx_her{:,4},0.05)
    max(h)
    
    mean(ctx_her{:,3})
    std(ctx_her{:,3})
     
    
    heri_ct = zeros(1,20484);
    for i = 1:100
        heri_ct(:,find(parcels200==i+1)) = ctx_her{i,3};
    end
    for i = 1:100
        heri_ct(:,find(parcels200==i+1001)) = ctx_her{i+100,3};
    end
    
    ctx2 = cbrewer('seq','RdPu',10);
    f = figure;
    BoSurfStatViewData(heri_ct,SN,'h2')
    colormap(ctx2)    
    exportfigbo(f,[RPATH 'F2.her.CTX.png'],'png', 10)
    close(f)
    
end

for load_her_area = 1
    area_her = readtable([dir '/affect_cognition/solar_SA/sa_heritability.csv']);
    area_her.name = cellfun(@(x) x(1:end-6), area_her.name, 'un', 0);
    
    idx = zeros(size(names200));
    for n = 1:length(names200)
        idx(n,1) = find(endsWith(area_her.name, names200{n}));
    end
    area_her = area_her(idx,:);
    writetable(area_her,[RPATH 'area_heritability.csv'], 'WriteRowNames', true);
    
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(area_her{:,4},0.05)
    max(h)
    
    mean(area_her{:,3})
    std(area_her{:,3})
     
    
    heri_ct = zeros(1,20484);
    for i = 1:100
        heri_ct(:,find(parcels200==i+1)) = area_her{i,3};
    end
    for i = 1:100
        heri_ct(:,find(parcels200==i+1001)) = area_her{i+100,3};
    end
    
        ctx2 = cbrewer('seq','OrRd',10);
        f = figure;
        BoSurfStatViewData(heri_ct,SN,'h2')
        colormap(ctx2)
        exportfigbo(f,[RPATH 'F2.her.area.png'],'png', 10)
        close(f)
end

for load_her_sub = 1
    sub_her = readtable(['affect_cognition/solar_subc/subc_heritability.csv']);
    writetable(sub_her,[RPATH 'sub_heritability.csv'], 'WriteRowNames', true);
    
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(sub_her{:,4},0.05)
    max(h)
    
    mean(sub_her{:,3})
    std(sub_her{:,3})
    
    f = figure;
    bar(sub_her{:,3})
    set(gca, 'XTick', 1:(length(subcort_lina.Properties.VariableNames)));
    set(gca, 'XTickLabel', subcort_lina.Properties.VariableNames)
    xtickangle(-45)
    exportfigbo(f,[RPATH, 'f2.her.sub.png'],'png', 10)
    close(f)
    
end

%genetic correlation
for load_coher_CTX = 1
    CTX_CH = readtable([dir 'affect_cognition/solar_CT/ctXcog_aff_gencorr.csv']);
    
    %rearrange order of traits, so it matches with order above
    CTX_CH = CTX_CH([1001:1200 401:600 201:400 801:1000 601:800 1:200],:);
       
    
        for t = 1:6
            if t == 1;
                trait = CTX_CH(1:t*200,:);
            else
                trait = CTX_CH(((t-1)*200)+1:t*200,:);
            end
            
            %rearrange order of parcels, so it matches with order above
            idx = zeros(size(names200));
            for n = 1:length(names200)
                idx(n,1) = find(startsWith(trait.name, [names200{n} '_']));
            end
            trait = trait(idx,:);
            
            % write table of overlap between phenotypic and genetic corr.
            rn = names200(var_ctx_FDR(t,:)==1);
            T = table(trait{find(var_ctx_FDR(t,:)==1),4:7},...
             'RowNames',rn);
            %T.Properties.VariableNames([1 2 3 4]) = {'rE' 'rp' 'rG' 'p'};
            writetable(T, [RPATH, 'gencorr_', trait{1,1}{1}, '_ctx.csv'], 'WriteRowNames', true)           
            
            [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(trait{:,7},0.05);
            max(h);
       
            
            heri_ct = zeros(1,20484);
            for x = 1:100
                heri_ct(:,find(parcels200==x+1)) = trait{x,6}*h(x);
            end
            for x = 1:100
                heri_ct(:,find(parcels200==x+1001)) = trait{x+100,6}.*h(x+100);
            end
            
            ctx2 = flipud(cbrewer('div','RdBu',11));
            f = figure
            BoSurfStatViewData(heri_ct,SN,[trait{1,1},' FDR'])
            colormap(ctx2)
            SurfStatColLim([-0.5 0.5])
            exportfigbo(f,[RPATH, 'gencorr', trait{1,1}{1}, 'wholebrainctxFDR.png'],'png', 10)
            close(f)
            
            heri_ct = zeros(1,20484);
            for x = 1:100
                heri_ct(:,find(parcels200==x+1)) = trait{x,6};
            end
            for x = 1:100
                heri_ct(:,find(parcels200==x+1001)) = trait{x+100,6};
            end
            
            ctx2 = flipud(cbrewer('div','RdBu',11));
            f = figure
            BoSurfStatViewData(heri_ct,SN,trait{1,1})
            colormap((ctx2))
            SurfStatColLim([-0.5 0.5])
            exportfigbo(f,[RPATH, 'gencorr', trait{1,1}{1}, 'wholebrainctx.png'],'png', 10)
            close(f)
            
            % plot only with significant phenotypic and genetic correlations
            idx_rE = var_ctx_FDR(t,:)'.*h;            
            heri_ct = zeros(1,20484);
            for x = 1:100
                heri_ct(:,find(parcels200==x+1)) = trait{x,6}*idx_rE(x);
            end
            for x = 1:100
                heri_ct(:,find(parcels200==x+1001)) = trait{x+100,6}.*idx_rE(x+100);
            end
            
            ctx2 = flipud(cbrewer('div','RdBu',11));
            f = figure;
            BoSurfStatViewData(heri_ct,SN,[trait{1,1},' FDR'])
            colormap(ctx2)
            SurfStatColLim([-0.5 0.5])
            exportfigbo(f,[RPATH, 'pheno_gencorr', trait{1,1}{1}, 'wholebrainctxFDR.png'],'png', 10)
            close(f)
            
        end
    
end

for load_coher_area = 1
    area_CH = readtable([dir 'affect_cognition/solar_SA/saXcog_aff_gencorr.csv']);
    
    %rearrange order of traits, so it matches with order above
    area_CH = area_CH([1001:1200 401:600 201:400 801:1000 601:800 1:200],:);
    
    
    for t = 1:6
        if t == 1
            trait = area_CH(1:t*200,:);
        else
            trait = area_CH(((t-1)*200)+1:t*200,:);
        end
        
        %rearrange order of parcels, so it matches with order above
        idx = zeros(size(names200));
        for n = 1:length(names200)
            idx(n,1) = find(startsWith(trait.name, [names200{n} '_']));            
        end
        trait = trait(idx,:);
        
        % write table of overlap between phenotypic and genetic corr.
        rn = names200(var_area_FDR(t,:)==1);
        T = table(trait{find(var_area_FDR(t,:)==1),4:7},...
            'RowNames',rn);
        writetable(T, [RPATH, 'gencorr_', trait{1,1}{1}, '_area.csv'], 'WriteRowNames', true)
        
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(trait{:,7},0.05);
        max(h);
        
        heri_ct = zeros(1,20484);
        for x = 1:100
            heri_ct(:,find(parcels200==x+1)) = trait{x,6}*h(x);
        end
        for x = 1:100
            heri_ct(:,find(parcels200==x+1001)) = trait{x+100,6}.*h(x+100);
        end
        
        
        
        ctx2 = flipud(cbrewer('div','RdBu',11));
        f = figure
        BoSurfStatViewData(heri_ct,SN,[trait{1,1},' FDR'])
        colormap((ctx2))
        SurfStatColLim([-0.5 0.5])
        exportfigbo(f,[RPATH, 'gencorr', trait{1,1}{1}, 'wholebrainareaFDR.png'],'png', 10)
        close(f)
        
        heri_ct = zeros(1,20484);
        for x = 1:100
            heri_ct(:,find(parcels200==x+1)) = trait{x,6};
        end
        for x = 1:100
            heri_ct(:,find(parcels200==x+1001)) = trait{x+100,6};
        end
        
        ctx2 = flipud(cbrewer('div','RdBu',11));
        f = figure,
        BoSurfStatViewData(heri_ct,SN,trait{1,1})
        colormap((ctx2))
        SurfStatColLim([-0.5 0.5])
        exportfigbo(f,[RPATH, 'gencorr', trait{1,1}{1}, 'wholebrainarea.png'],'png', 10)
        close(f)
        
            % plot only with significant phenotypic and genetic correlations
            idx_rE = var_area_FDR(t,:)'.*h;            
            heri_ct = zeros(1,20484);
            for x = 1:100
                heri_ct(:,find(parcels200==x+1)) = trait{x,6}*idx_rE(x);
            end
            for x = 1:100
                heri_ct(:,find(parcels200==x+1001)) = trait{x+100,6}.*idx_rE(x+100);
            end
            
            ctx2 = flipud(cbrewer('div','RdBu',11));
            f = figure,
            BoSurfStatViewData(heri_ct,SN,[trait{1,1},' FDR'])
            colormap(ctx2)
            SurfStatColLim([-0.5 0.5])
            exportfigbo(f,[RPATH, 'pheno_gencorr', trait{1,1}{1}, 'wholebrainareaFDR.png'],'png', 10)
            close(f)        
        
    end
    
end

for load_coher_sub = 1
    subc = readtable([dir 'affect_cognition/solar_subc/subcXcog_aff_gencorr.csv']);
    
    %rearrange order of traits, so it matches with order above
    subc = subc([81:96 33:48 17:32 65:80 49:64 1:16],:);
       
    
        for t = 1:6
            if t == 1
                trait = subc(1:t*16,:);
            else
                trait = subc(((t-1)*16)+1:t*16,:);
            end
                        
            % write table of overlap between phenotypic and genetic corr.
            rn = subcort_lina.Properties.VariableNames(var_sub_FDR(t,:)==1);
            T = table(trait{find(var_sub_FDR(t,:)==1),4:7},...
             'RowNames',rn);
            %T.Properties.VariableNames([1 2 3 4]) = {'rE' 'rp' 'rG' 'p'};
            writetable(T, [RPATH, 'gencorr_', trait{1,1}{1}, '_subc.csv'], 'WriteRowNames', true)            
            
            [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(trait{:,7},0.05);
            max(h);       
            
        end
    
end


%% Supplementary Figures: Correlation matrix of all behavioural measures
for fig_s1 = 1
    % BEHAVIORAL MEASURES    
    
    list_of_vars = [HCP.Flanker_Unadj, HCP.CardSort_Unadj, HCP.PicSeq_Unadj, HCP.ListSort_Unadj, HCP.PMAT24_A_RTCR,...
        HCP.PicVocab_Unadj, HCP.ReadEng_Unadj,...
        HCP.LifeSatisf_Unadj,HCP.MeanPurp_Unadj,HCP.PosAffect_Unadj,...
        HCP.AngAffect_Unadj, HCP.AngHostil_Unadj, HCP.FearAffect_Unadj, HCP.Sadness_Unadj, HCP.PercStress_Unadj];
    
    comptitlevar = {'Flanker', 'Card sort.','Picture seq.','List sort.','PMAT','Picture vocab.','Reading',...
        'Life satisfaction','Mean purpose','Positive affect',...
        'Anger-affect','Anger-hostil.','Fear-affect','Sadness','Perc. stress'};
    

   % Run linear model for each behavioral variable to create correlation matrix
   keep = studykeep;
   [idx, cidx] = find(isnan(list_of_vars(keep,:)));
   keep(idx) = [];
   beta_var2 = zeros(15);
   
    for i = 1:size(list_of_vars,2)        
        vark   = zscore(list_of_vars(keep,i));
        agek   = zscore(HCP_r.Age_in_Yrs(keep));
        sexk   = cellstr(HCP.Gender(keep));
        sex_num = grp2idx(sexk);
        
        M    = 1 + term(agek) + term(sexk) + (term(agek) * term(sexk)) + (term(agek) * term(agek))+ term(vark);
        slm  = SurfStatLinMod(zscore(list_of_vars(keep,:)),M);
        slm  = SurfStatT(slm, vark);
        t_var2(i,:) = slm.t;
        beta_var2(i,:) = slm.coef(end,:);
        pp_var2(i,:) = 1- tcdf(slm.t,slm.df);
        pn_var2(i,:) = 1- tcdf(-slm.t,slm.df);     
    end
    
    t_var2(eye(15)==1) = 0;
    beta_var2(eye(15)==1) = 0;
    pp_var2(eye(15)==1) = 0;
    pn_var2(eye(15)==1) = 0;
    
    pp_u = pp_var2(triu(ones(15),1)==true);
    pn_u = pn_var2(triu(ones(15),1)==true);
    fdr_bh(pp_u, 0.025)
    ppn = squareform(ans);
    fdr_bh(pn_u, 0.025)
    pnn = squareform(ans);   
    ppp = ppn+pnn;
    
    ppp(eye(15)==1) =1;
     
        
   %make correlation matrix plot 
   for explo = 1 
       psignr = ppp == 1;
        %rleft = t_var2.*psignr;
        %rleft = beta_var2.*psignr;
        rleft = beta_var2;
        mat = rleft;
        f = figure;
        set(gcf, 'Position',  [100, 100, 1200, 800])
        imagesc(mat)
        colorbar;
        caxis([-1,1]);
        textStrings = num2str(mat(:), '%0.2f');       % Create strings from the matrix values
        textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
        textStrings = replace(textStrings,'0.00','-');
        [x, y] = meshgrid(1:15);  % Create x and y coordinates for the strings
        hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
            'HorizontalAlignment', 'center');
        midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
        textColors = repmat(abs(mat(:)) > 0.5, 1, 3);  % Choose white or black for the text color of the strings so they can be easily seen over the background color
        
        set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors
        set(gca, 'XTick', 1:(size(list_of_vars,2)));
        set(gca, 'XTickLabel', comptitlevar)
        xtickangle(-45)
        colormap(flipud(cbrewer('div','RdBu',99)));
        colorbar;
        set(gca, 'YTick', 1:(size(list_of_vars,2)));
        set(gca, 'YTickLabel', comptitlevar)
        title('Phenotypic correlation: standardized beta values')
        
        print(f,[RPATH 'Corr_all_values'],'-dpng')
        close(f)
   end
end

for figs2 = 1
    % plot brain data    
    buckner = ([102 51 51; 51 51 102; 102 102 153; 153 153 204; 153 255 51; ...
        102 204 51; 255 255 51; 255 153 51; 255 102 51; 204 51 0])./255;
    
    % mean thickness on surface
    heri_ct = zeros(1,20484);
    for i = 1:100
        heri_ct(:,find(parcels200==i+1)) = mean(HCP200_CT(studykeep,i));
    end
    for i = 1:100
        heri_ct(:,find(parcels200==i+1001)) = mean(HCP200_CT(studykeep,i+100));
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(buckner)
    SurfStatColLim([1 4])
    exportfigbo(f,[RPATH, 'F0.HCP.ctx.png'],'png', 10)
    close(f)
    
    f = figure,
    hist( mean(HCP200_CT(studykeep,:)))
    exportfigbo(f,[RPATH, 'F0.HCP.total.ctx.png'],'png', 10)
    close(f)
    
    % mean surface area on surface
    heri_ct = zeros(1,20484);
    for i = 1:100
        heri_ct(:,find(parcels200==i+1)) = mean(area200(i,studykeep));
    end
    for i = 1:100
        heri_ct(:,find(parcels200==i+1001)) = mean(area200(i+100,studykeep));
    end
    
    f = figure,
    BoSurfStatViewData(heri_ct,SN,'')
    colormap(cbrewer('seq','YlGnBu',99))
    SurfStatColLim([500 1700])
    exportfigbo(f,[RPATH, 'F0.area.HCP.png'],'png', 10)
    close(f)
    
    f = figure,
    hist( sum(area200(:,studykeep)))
    exportfigbo(f,[RPATH, 'F0.area.HCP.total.png'],'png', 10)
    close(f)
end
