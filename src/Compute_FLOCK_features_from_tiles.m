%store_patient_info_tiles = x.tiles_ZOOM2000
%cluster_names2 = {'Tumor' ; 'Stroma' ; 'Immune cells'};


function [ feat_tile_ALL , feat_tile_TUMOR_STROMA  ,  feat_tile_TUMOR_IMMUNE     ,  feats_names] = ...
    Compute_FLOCK_features_from_tiles(store_patient_info_tiles , cluster_names2  , radius_vector, names_var)

    %%% SOS apo haralick_no_img
    feats_names = load('fedeg_names.mat');
    feats_names = feats_names.fedeg_gman_names;

    feat_tile_ALL=struct;
    feat_tile_TUMOR_STROMA=struct;
    feat_tile_TUMOR_IMMUNE=struct;
    
    for j=1:length(feats_names)
        feat_tile_ALL.(feats_names{j}) = NaN;
        feat_tile_TUMOR_STROMA.(feats_names{j}) = NaN;
        feat_tile_TUMOR_IMMUNE.(feats_names{j}) = NaN;
    end
    
    
    %%%%% c =y , r =x
    para = [];
    para.feature_space='Centroid-Area';%%%%para.feature_space='Centroid-Area-MeanIntensity'; %para.feature_space='Centroid-MeanIntensity';
    para.bandWidth_space=80;% bandwidth in the spacial space, higher of the value, bigger of the FeDeG
    para.bandWidth_features=[100];% bandwidth in the corresponding feature space
    %para.debug=1; % turn the debug mode on
    para.shownuclei_thiness=1; % the visual effect of nuclei 
    para.shownucleiconnection_thiness=1; % the visual effect of FeDeG
    para.debug=0;
    para.num_fixed_types=0; 
    
    for KK = 1 : size(store_patient_info_tiles,1)
        %disp(KK)
        para.I = store_patient_info_tiles{KK,1}.image;
        
        para.data_other_attribute=[];
        para.clust2types=[];
        para.typeCent=[];
        
        %------------------------  Tumor Stroma ------------------------------------
        [mask_for_rad ,bounds ,para.nuclei , index] = create_bounds_per_celltype(KK , store_patient_info_tiles , 'TUMOR_STROMA');
        
        if(max(max(mask_for_rad))>0 && size(bounds,2)>1)
            
            properties  =    store_patient_info_tiles{KK,1}.props(index);
            [~, idx] = ismember([properties.cellClass], cluster_names2);
            for zzz = 1 : size(properties,1)
                properties(zzz,1).Centroid = [ properties(zzz,1).cellCentroidX   properties(zzz,1).cellCentroidY  ];
                properties(zzz,1).Classified = idx(zzz)*100; 
            end
            para.properties=properties; % you can check out the nuclear properites 
            
            
            [clustCent,~,cluster2dataCell,para.data_other_attribute,para.clust2types,para.typeCent]=Lconstruct_FeDeG_v3(para.I,para);
            [~,~  ,  fedeg_gman_features ,  ~]= L_get_FeDeG_features_v3_no_pheno(clustCent,cluster2dataCell,para,[para.properties.cellClass], radius_vector, names_var);
            
            for iii = 1:length(feats_names)
                feat_tile_TUMOR_STROMA.(feats_names{iii}) = [feat_tile_TUMOR_STROMA.(feats_names{iii}) fedeg_gman_features{iii} ] ;
            end
        end
        %-----------------------------------------------------------------------------------------------
        
        para.data_other_attribute=[];
        para.clust2types=[];
        para.typeCent=[];
        
        %------------------------  TUMOR_IMMUNE ------------------------------------
        [mask_for_rad ,bounds ,para.nuclei , index] = create_bounds_per_celltype(KK , store_patient_info_tiles , 'TUMOR_IMMUNE');
        
        if(max(max(mask_for_rad))>0 && size(bounds,2)>1)
            
            properties  =    store_patient_info_tiles{KK,1}.props(index);
            [~, idx] = ismember([properties.cellClass], cluster_names2);
            for zzz = 1 : size(properties,1)
                properties(zzz,1).Centroid = [ properties(zzz,1).cellCentroidX   properties(zzz,1).cellCentroidY  ];
                properties(zzz,1).Classified = idx(zzz)*100; 
            end
            para.properties=properties; % you can check out the nuclear properites 
            
            
            [clustCent,~,cluster2dataCell,para.data_other_attribute,para.clust2types,para.typeCent]=Lconstruct_FeDeG_v3(para.I,para);
            [~,~  ,  fedeg_gman_features ,  ~]= L_get_FeDeG_features_v3_no_pheno(clustCent,cluster2dataCell,para,[para.properties.cellClass], radius_vector, names_var);
            
            for iii = 1:length(feats_names)
                feat_tile_TUMOR_IMMUNE.(feats_names{iii}) = [feat_tile_TUMOR_IMMUNE.(feats_names{iii}) fedeg_gman_features{iii} ] ;
            end
        end
        %-----------------------------------------------------------------------------------------------
        
        para.data_other_attribute=[];
        para.clust2types=[];
        para.typeCent=[];
        
        %------------------------  ALL ------------------------------------
        [mask_for_rad ,bounds ,para.nuclei , index] = create_bounds_per_celltype(KK , store_patient_info_tiles , 'ALL');
        
        if(max(max(mask_for_rad))>0 && size(bounds,2)>1)
            
            properties  =    store_patient_info_tiles{KK,1}.props(index);
            [~, idx] = ismember([properties.cellClass], cluster_names2);
            for zzz = 1 : size(properties,1)
                properties(zzz,1).Centroid = [ properties(zzz,1).cellCentroidX   properties(zzz,1).cellCentroidY  ];
                properties(zzz,1).Classified = idx(zzz)*100; 
            end
            para.properties=properties; % you can check out the nuclear properites 
            
            
            [clustCent,~,cluster2dataCell,para.data_other_attribute,para.clust2types,para.typeCent]=Lconstruct_FeDeG_v3(para.I,para);
            [~,~  ,  fedeg_gman_features ,  ~]= L_get_FeDeG_features_v3_no_pheno(clustCent,cluster2dataCell,para,[para.properties.cellClass], radius_vector, names_var);
            
            for iii = 1:length(feats_names)
                feat_tile_ALL.(feats_names{iii}) = [feat_tile_ALL.(feats_names{iii}) fedeg_gman_features{iii} ] ;
            end
        end
        %-----------------------------------------------------------------------------------------------
        
    end     
        %
        %[mask_for_rad ,bounds ,para.nuclei , index] = create_bounds_per_celltype(KK , store_patient_info_tiles , 'ALL');
        %store_patient_info_tiles{KK,1}.cell_class(find(cell_cluster_number == 1))
        %store_patient_info_tiles{KK,1}.cell_class(find(cell_cluster_number == 2))
        %store_patient_info_tiles{KK,1}.cell_class(find(cell_cluster_number == 3))

end
