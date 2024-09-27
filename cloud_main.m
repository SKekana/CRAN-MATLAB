classdef cloud_main<handle
    
    
    properties( GetAccess = 'public', SetAccess = 'public' )
      RB_sectors;
      BICM_capacity_tables_path='data_files/BICM_capacity_tables_20000_realizations.mat';
      BICM_capacity_tables;
      CQI_SINR_mapping;
      users_allocated_RBs_SINR_values;
      users_allocated_RBs_capacity_values;
    end
    
    methods
        function obj=cloud_main()
        end
        function [cellsArray]=creating_cells(obj,numOfTiers,cellRad)
            [cellCenters,cells_vertices,numOfCells]=cells_centers_array(numOfTiers, cellRad);
            cellsArray=sectors(cellCenters,cells_vertices);
            [minimumX,maximumX]=obj.vertix(cells_vertices,1);
            [minimumY,maximumY]=obj.vertix(cells_vertices,2);
            ROI_x=[minimumX maximumX];
            ROI_y=[minimumY maximumY];
            [num_RRH]=obj.numberOfRRH(numOfCells);
            [macro_RRH_positions]=macro_RRH(num_RRH/3,cells_vertices,cellRad);
            [micro_RRH_positions]=micro_RRH(minimumX,maximumX,cells_vertices,cellRad,macro_RRH_positions);
            pico_RRH(cellCenters);
            [r1 c1]=size(macro_RRH_positions);
            [r2 c2]=size(micro_RRH_positions);
            [r3 c3]=size(cellCenters);
            total_RRH_positions(1:r1,1:2)=macro_RRH_positions;
            total_RRH_positions(r1+1:r1+r2,1:2)=micro_RRH_positions;
            total_RRH_positions(r1+r2+1:r1+r2+r3,1:2)=cellCenters;
            channel_model=channelModel();
             RRHs_array=RRH_cloud(total_RRH_positions,r1,r2,r3,channel_model);
             RRHs_array(1:r1,1).RRH_sectors(cellCenters,cellsArray,cellRad);
             channel_model.calculate_pathloss_maps(RRHs_array(1:r1,1),ROI_x,ROI_y);
%              figure(2)
%              imagesc(ROI_x,ROI_y,channel_model.pathloss(:,:,1));
%              set(gca,'YDir','normal');
%              colorbar;
%              RRHs_array.plotting_RRH_ID();
             channel_model.shadowFadingMapClaussen( channel_model.shadow_fading_map_resolution,ROI_x,ROI_y,channel_model.shadow_fading_n_neighbors,RRHs_array(1:r1,1),channel_model.shadow_fading_mean,channel_model.shadow_fading_sd,channel_model.inter_RRH_correlation);
             RX_power_all_sectors=RRHs_array(1,1).antenna_tx_power./10.^(channel_model.pathloss/10)./10.^(channel_model.shadowFading.pathlossNew_times_3/10);
             thermal_noise_w=channel_model.calculate_thermal_noise_w();
             [SNR_dB SINR_dB]=obj.SNR_SINR(RX_power_all_sectors,thermal_noise_w);
             RRHs_array(13).antenna
             array_users=user_cloud(100,minimumX,maximumX,minimumY,maximumY);
             array_users.RRH_sector_ID(1000,RRHs_array(1:r1),cellCenters,cellsArray);
             cellsArray.sectorFirstTierNeighbour(cellRad);
             obj.sectorsResourcesInfo(cellsArray,1);
             obj.roundRobin(cellsArray);
             obj.load_BICM_tables();
             users_array_positions=[array_users(:,1).xPosition;array_users(:,1).yPosition]';
             users_positions_in_pixels=obj.pos_to_pixel(users_array_positions,[ROI_x(1) ROI_y(1)],channel_model.data_resolution);
%              circle(array_users(1).xPosition,array_users(1).yPosition,1000);
            array_users(1).user_near_RRHs_Info(1,1)
             
                           
        end 
        function [num_RRH]=numberOfRRH(obj,numOfCells)
            factor=1;
            incrementer=3;
            currentNum=factor*3;
            while currentNum<numOfCells
                factor=factor+incrementer;
                incrementer=incrementer+2;
                currentNum=factor*3;
            end
            num_RRH=currentNum;
          
        
        end
       
        function [minimum,maximum]=vertix(obj,cells_vertices,xORy)
            if(xORy==1)
                start=1;
            else
                start=2;
            end
            begin =1;
            [rows col]=size(cells_vertices);
            while(start<=col)
                vector(1,begin:begin+rows-1)=transpose(cells_vertices(1:rows,start));
                begin=begin+rows;
                start=start+2;
            end
            minimum=min(vector);
            maximum=max(vector);
          
        end
        function [SNR_dB ,SINR_dB]=SNR_SINR(obj,RX_power_all_sectors,thermal_noise_w)
            for ii=1:1:size(RX_power_all_sectors,3)
             SNR_dB(:,:,ii)=10*log10(RX_power_all_sectors(:,:,ii)./thermal_noise_w);
             SINR_dB(:,:,ii)=10*log10(RX_power_all_sectors(:,:,ii)./(sum(RX_power_all_sectors,3)-RX_power_all_sectors(:,:,ii)+thermal_noise_w));
             end
        end
        function sectorsResourcesInfo(obj,sectors,reuse)
            for ii=1:1:length(sectors)
                switch reuse
                    case 1
                        sectors(ii,1).Resources_Info.RB_number=100;
                        sectors(ii,1).Resources_Info.bandwidth=20;
                        sectors(ii,1).Resources_Info.total_RB_allocation='false';
                        sectors(ii,1).Resources_Info.RB_array_allocation_indictor=[];
                        sectors(ii,1).Resources_Info.frequency='f1';
                        sectors(ii,1).Resources_Info.area_covered='total area';
                        sectors(ii,1).Resources_Info.numberOfAllocatedRB=0;
            end
            end
        end
        function roundRobin(obj,sectors)
            obj.RB_sectors=cell(1,sectors(1,1).Resources_Info.RB_number);
            for ii=1:1:length(sectors)
                sectors(ii,1).Resources_Info.RB_allocation_type='round robin';
                sectors(ii,1).numOfUsersInSector=length(sectors(ii).users_array);
                total_RB_allocated_per_sector=min(sum([sectors(ii).users_array(1,:).user_Required_RB_number]),sectors(ii).Resources_Info.RB_number);
                sectors(ii,1).Resources_Info.numberOfAllocatedRB=total_RB_allocated_per_sector;
                RB_allocated=1;
                while RB_allocated<=total_RB_allocated_per_sector
                    for jj=1:1:sectors(ii).numOfUsersInSector
                    if(sectors(ii).users_array(jj).userCurrentAllocatedResources<sectors(ii).users_array(jj).user_Required_RB_number &&RB_allocated<=total_RB_allocated_per_sector)
                    sectors(ii).users_array(jj).userCurrentAllocatedResources=sectors(ii).users_array(jj).userCurrentAllocatedResources+1;
                    sectors(ii).users_array(jj).userAllocatedResources=[sectors(ii).users_array(jj).userAllocatedResources RB_allocated];
                    obj.RB_sectors{1,RB_allocated}=[obj.RB_sectors{1,RB_allocated} ii];
                    RB_allocated=RB_allocated+1;
                    end
                end
            end
        end
    
    
    
        end
    function [ pos_pixel, pos_pixel_exact] = pos_to_pixel( obj,pos,roi_min, data_res)
                   pos_pixel(:,1) = floor((pos(:,1)-roi_min(1))/data_res)+1;
                   pos_pixel(:,2) = floor((pos(:,2)-roi_min(2))/data_res)+1;
                   pos_pixel_exact(:,1) = ((pos(:,1)-roi_min(1))/data_res)+1;
                   pos_pixel_exact(:,2) = ((pos(:,2)-roi_min(2))/data_res)+1;
    end
    function CQI_SINR_mapping= CQI_equivalent_SINR_vector(obj)
        CQI_SINR_mapping(1,1)=-7;
        CQI_SINR_mapping(1,2)=-5;
        CQI_SINR_mapping(1,3)=-3;
        CQI_SINR_mapping(1,4)=-1;
        CQI_SINR_mapping(1,5)=1;
        CQI_SINR_mapping(1,6)=3;
        CQI_SINR_mapping(1,7)=5;
        CQI_SINR_mapping(1,8)=7;
        CQI_SINR_mapping(1,9)=9;
        CQI_SINR_mapping(1,10)=11;
        CQI_SINR_mapping(1,11)=13;
        CQI_SINR_mapping(1,12)=15;
        CQI_SINR_mapping(1,13)=17;
        CQI_SINR_mapping(1,14)=19;
        CQI_SINR_mapping(1,15)=20;
    end
    function CQI_BICM_table_mapping=CQI_to_BICM_table_mapping(obj)
        CQI_BICM_table_mapping(1,1:6)=1;
        CQI_BICM_table_mapping(1,7:9)=2;
        CQI_BICM_table_mapping(1,10:15)=3;
    end
    function CQI_to_SINR_mapping(obj,users_array)
        for ii=1:1:length(users_array)
            obj.users_allocated_RBs_SINR(ii,:)=obj.CQI_SINR_mapping(users_array(ii,1).userAllocated_RB_CQI);
            obj.users_allocated_RBs_modulation_type(ii,:)=obj.CQI_BICM_table_mapping(users_array(ii,1).userAllocated_RB_CQI);
        end
           
    end
    function load_BICM_tables(obj)
        load(obj.BICM_capacity_tables_path);
        BICM_capacity_tables(1).I(1)=0;
        BICM_capacity_tables(1).I(end)=ceil(BICM_capacity_tables(1).I(end));
        BICM_capacity_tables(2).I(1)=0;
        BICM_capacity_tables(2).I(end)=ceil(BICM_capacity_tables(1).I(end));
        BICM_capacity_tables(3).I(1)=0;
        BICM_capacity_tables(3).I(end)=ceil(BICM_capacity_tables(1).I(end));
        obj.BICM_capacity_tables=BICM_capacity_tables;
    end
    function SINR_to_capacity_mapping(obj)
       for ii=1:1:size(obj.users_allocated_RBs_SINR,1)
           for jj=1:1:size(obj.users_allocated_RBs_SINR(ii,:),2)
               table_index=obj.user_allocated_RBs_modulation_type(ii,jj);
               snr_index=find(obj.BICM_capacity_tables(table_index).SNR==obj.users_allocated_RBs_SINR(ii,jj));
               obj.users_allocated_RBs_capacities(ii,jj)=obj.BICM_capacity__tables(table_index).I(snr_index);
           end
           obj.users_average_capacity(ii,1)=sum(users_allocated_RBs_capacities(ii,:))/length(users_allocated_RBs_capacities(ii,:));
           table1_capacity_index=find(obj.BICM_capacity_tables(1).I==obj.users_average_capacity(ii,1));
           if(isempty(table1_capacity_inex))
               capacity_difference=obj.BICM_capacity_tables(1).I-obj.users_average_capacity(ii,1);
               table1_capacity_index=find(capacity_difference==(max(capacity_difference(find(capacity_difference<0)))));
              
           end
            SNR_effective_vector=[SNR_effective_vector obj.BICM_capacity_tables(1).SNR(capacity_index)];
           table2_capacity_index=find(obj.BICM_capacity_tables(2).I==obj.users_average_capacity(ii,1));
           if(isempty(table2_capacity_inex))
               capacity_difference=obj.BICM_capacity_tables(2).I-obj.users_average_capacity(ii,1);
               table2_capacity_index=find(capacity_difference==(max(capacity_difference(find(capacity_difference<0)))));
               
           end
           SNR_effective_vector=[SNR_effective_vector obj.BICM_capacity_tables(2).SNR(capacity_index)];
           table3_capacity_index=find(obj.BICM_capacity_tables(3).I==obj.users_average_capacity(ii,1));
           if(isempty(table3_capacity_inex))
               capacity_difference=obj.BICM_capacity_tables(2).I-obj.users_average_capacity(ii,1);
               table3_capacity_index=find(capacity_difference==(max(capacity_difference(find(capacity_difference<0)))));
               
           end
           SNR_effective_vector=[SNR_effective_vector obj.BICM_capacity_tables(2).SNR(capacity_index)];
           obj.users_effective_SINR_vector(ii,1)=min(SNR_effective_vector);
       end
    end
    end
end