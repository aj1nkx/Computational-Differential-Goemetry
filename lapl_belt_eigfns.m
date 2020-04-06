function lapl_belt_eigfns()
    % Add the library paths
    % Parts of this code are borrowed from Anand Joshi @ USC
    addpath(genpath('~/Downloads/toolbox_graph'));
    addpath('./surface_laplacian_sample_code/')
    n_evecs = 10;
    n_clust = 4;
    n_surf = 1;
    
    % Load Laplacians from the surfaces
    all_surfs = cell(1,n_surf);
    temp_load = load('s11.mat');
    all_surfs{1} = temp_load.s11;
    
    % Iterating over the different surfaces
    for i = 1:length(all_surfs)
        fprintf('Surface %d\n',i);
        % Find and plot the eigenvectors
        Si = all_surfs{i};
        [Vj_temp, Dj_temp] = eigs(surface_laplacian(all_surfs{i}), n_evecs, 'sm');
        [~, ind_sort] = sort(diag(Dj_temp),'ascend');
        V = Vj_temp(:,ind_sort);
        D = Dj_temp(:,ind_sort);
        grd_width = floor(sqrt(n_evecs));
        figure;
        for jj = 1:n_evecs
            subplot(grd_width+1,grd_width,jj); patch('faces',Si.faces,'vertices',Si.vertices,'facevertexcdata',V(:,jj),'edgecolor','none','facecolor','interp');
            axis equal; material dull;view(90,0);camlight;material dull;axis off;axis tight;
        end
        % title(['First ',num2str(n_evecs),' eigenvectors of the LB operator.']);
        saveas(gcf,['./data/Evecs_',num2str(i),'.png']);
        fprintf('Finished computing and plotting eigenvectors\n');
        
        % Compute the GPS co-ordinates
        gps_coords = V(:,2:n_evecs)*pinv(sqrt(D(2:n_evecs,2:n_evecs)));

        % Segment the region into 4 parts and visualize it
        cluster_id = kmeans(gps_coords, n_clust);
        figure;
        patch('faces',Si.faces,'vertices',Si.vertices,'facevertexcdata',cluster_id,'edgecolor','none','facecolor','interp');
        axis equal; material dull; view(90,0); camlight; material dull; axis off; axis tight;
        % title('Segmented surface');
        saveas(gcf,['./data/Segmentation_',num2str(i),'.png']);
        fprintf('Finished clustering\n');

        % Project the 3d coordinates into some percentage of basis vectors
        app_perc = [0.02, 0.1, 0.3, 0.7];
        figure;
        for j = 1:length(app_perc)
            fprintf('%f fraction of eigenvectors\n',app_perc(j));
            % Do the eigendecomposition
            [Vj_temp, Dj_temp] = eigs(surface_laplacian(all_surfs{i}), round(size(V,1)*app_perc(j)), 'sm');
            [~, ind_sort] = sort(diag(Dj_temp),'ascend');
            Vj = Vj_temp(:,ind_sort);
            
            % Project and reconstruct X, Y and Z coordinates
            Xc_new = Vj*((Si.vertices(:,1)).'*Vj).';
            Yc_new = Vj*((Si.vertices(:,2)).'*Vj).';
            Zc_new = Vj*((Si.vertices(:,3)).'*Vj).';
            
            % Plot the reconstruction            
            subplot(1,length(app_perc),j);patch('faces',Si.faces,'vertices',[Xc_new, Yc_new, Zc_new],'edgecolor','none');
            axis equal; material dull;view(90,0);camlight;material dull;axis off;axis tight;
        end
        % title(['Surface approximations with ',num2str(app_perc(1)*100),',',num2str(app_perc(2)*100),',',num2str(app_perc(3)*100),',',num2str(app_perc(4)*100),' percent basis functions']);
        saveas(gcf,['./data/Surf_approx_',num2str(i),'.png']);
    end
    
    % Delete all the custom paths
    restoredefaultpath;
end