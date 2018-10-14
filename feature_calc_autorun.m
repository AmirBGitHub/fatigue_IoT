%Feature Script
%Seamus Lombardo (6-9-17)
%Amir Baghdadi (6-11-17)
%Seamus Lombardo (6-15-17)
%Amir Baghdadi (6-16-17)
%Seamus Lombardo (6-16-17)
%Seamus Lombardo (6-17-17)
%seamus lombardo (6-18-17)
%Amir Baghdadi (6-23-17)

% MMH
% Dorisha sampling rate: 46 Hz
% Alycia sampling rate: 49.61 Hz
% LarryR sampling rate: 50.92 Hz
% Greg sampling rate: 48.92 Hz

% WLK
% Dave sampling rate: 49.18 Hz
% Gary sampling rate: 49.61 Hz
% Shivam sampling rate: 46.93 Hz
% Lisa sampling rate: 49.80 Hz

clear all
clc
%% Create matrix of individual sheet names
% for each participant put the name of foot, back, and hip data files here 
% to be analyzed 

dir1 = 'D:/GitHub/Fatigue Changepoint GitHub/fatigue-changepoint/Data/Raw';
dir2 = 'D:/GitHub/Fatigue IoT GitHub/fatigue_IoT/data';

Subject = {'Jim' , 'Anastasi', 'Tarun', 'Angela', 'Darya', 'Gary', 'Holly', 'Larry R', 'Lisa', 'Nicholas', 'Po', 'Alycia', 'Dave', 'Dorisha', 'Sean'};
SubjectNum=[ 1   ,     2     ,    3   ,    4    ,    5   ,    6  ,    7   ,     8    ,   9   ,     10    ,  11 ,    12   ,   13  ,     14   ,   15   ];                     

% subject ordering in FPCA data set
% Po        1 #
% Marcus	2 
% Nicholas	3 #
% Alycia	4 #
% Larry R	5 #
% Greg      6 
% Gary      7 #
% Dorisha	8 #
% Dave      9 #
% Anastasia	10 #
% Jim       11 #
% Sean      12 #
% Shivam	13 
% Hamid     14 
% Lisa      15 #
% Thomas	16 
% Sevack	17 
% Brian     18 
% Kayla     19 
% Angela	20 #
% Darya     21 #

for sbjct = 6:15

    sheet_name_foot_mat = fullfile([dir2, '/MMH-', char(Subject(sbjct)), '_Shimmer_CEA2_Calibrated_SD.csv']);
    sheet_name_back_mat = fullfile([dir2, '/MMH-', char(Subject(sbjct)), '_Shimmer_CE9B_Calibrated_SD.csv']);

    samplingrate = 51.2/2;

    loadData = load([dir1, '/Subject', num2str(sbjct),'.mat']);
    segSteps = loadData.M_i_k;
    
    % spreading the segments into 2 min intervals
    twoMinSegCnt=1;
    k=1;
    for i = 1:size(segSteps,1)
        tS=segSteps{i,17}(1);
        tE=segSteps{i,17}(end);
        S_s_all(i) = ceil(samplingrate*tS);
        S_e_all(i) = ceil(samplingrate*tE);
        if tE>twoMinSegCnt*120
            segSteps_2min{twoMinSegCnt,:} = segSteps(k:i,:);
            S_s_2min{twoMinSegCnt,:} = S_s_all(k:i);
            S_e_2min{twoMinSegCnt,:} = S_e_all(k:i);
            k=i+1;
            twoMinSegCnt=twoMinSegCnt+1;
        end
    end
    
    % making segments available for all 2 min intervals
    if size(segSteps_2min,1) < 86
       for i = size(segSteps_2min,1)+1:86
           segSteps_2min{i,1} = NaN;
           S_s_2min{i,1} = NaN;
           S_e_2min{i,1} = NaN;
       end
    end
    

    % Read data

    % iterate through all sheets listed above
    % select sheet name
    sheet_name_foot = char(sheet_name_foot_mat);
    sheet_name_back = char(sheet_name_back_mat);

    % determine the sheet with the smallest size and use the end cell of that sheet as the end value for all trials   
    whole_foot_array = xlsread(sheet_name_foot);
    whole_back_array = xlsread(sheet_name_back);

    %%
    % Load Foot Data Non-Fatigue

    footData = downsample(whole_foot_array(600*samplingrate:end,2:7),2);
    backData = downsample(whole_back_array(600*samplingrate:end,2:7),2);

    size_foot_array = size(footData);
    size_back_array = size(backData);
    end_dimension_array_foot = size_foot_array(1);
    end_dimension_array_back = size_back_array(1);
    size_array = [end_dimension_array_foot end_dimension_array_back];
    min_size = min(size_array);

    foot_accel_x = footData(1:min_size,1);
    foot_accel_y = footData(1:min_size,2);
    foot_accel_z = footData(1:min_size,3);
    foot_vel_x = footData(1:min_size,4);
    foot_vel_y = footData(1:min_size,5);
    foot_vel_z = footData(1:min_size,6);

    back_accel_x = backData(1:min_size,1);
    back_accel_y = backData(1:min_size,2);
    back_accel_z = backData(1:min_size,3);
    back_vel_x = backData(1:min_size,4);
    back_vel_y = backData(1:min_size,5);
    back_vel_z = backData(1:min_size,6);

    % concatenate to fit Amir's filter (IMU_Kalman.m)
    imu_data_foot = [foot_accel_x foot_accel_y foot_accel_z foot_vel_x foot_vel_y foot_vel_z];
    imu_data_back = [back_accel_x back_accel_y back_accel_z back_vel_x back_vel_y back_vel_z];

    %% make 2-min portion of data to be analyze
    
    twoMinIntCnts = min([86, floor(min_size/floor(samplingrate*120))]);
    for n = 1:twoMinIntCnts
        fatigue_state = 'whole';
        twoMinint = n
        data_in = [imu_data_foot(floor(samplingrate*120)*(n-1)+1:floor(samplingrate*120)*(n),:), imu_data_back(floor(samplingrate*120)*(n-1)+1:floor(samplingrate*120)*(n),:)];
        dataCnt = range(floor(samplingrate*120)*(n-1)+1:floor(samplingrate*120)*(n))+1;
        [number_of_steps(n),step_time_average(n),step_time_SD(n),step_distance_average(n),... 
            step_distance_SD(n),step_mean_vel_average(n),step_peak_vel_average(n),step_mean_acc_average(n),...
            step_peak_acc_average(n),step_mean_jrk_average(n),step_peak_jrk_average(n),overall_average_angle_back_bent(n),...
            overall_SD_angle_back_bent(n),leg_average_rot_velocity(n),leg_SD_rot_velocity(n),vertical_impact_average(n),...
            vertical_impact_SD(n),back_rotational_pos_range_in_sag_plane(n),back_rotational_pos_mean_in_sag_plane(n),...
            back_rotational_pos_SD_in_sag_plane(n)] = ...
        Kinematics_feature(data_in,segSteps_2min{n,1},S_s_2min{n,1}-(n-1)*dataCnt+1,S_e_2min{n,1}-(n-1)*dataCnt+1,samplingrate,sheet_name_foot,n);
    end

    % write to excel sheet
    xlswrite('testsheet_whole.xlsx',number_of_steps(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'A2:A86');
    xlswrite('testsheet_whole.xlsx',step_time_average(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'B2:B86');
    xlswrite('testsheet_whole.xlsx',step_time_SD(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'C2:C86');
    xlswrite('testsheet_whole.xlsx',step_distance_average(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'D2:D86');
    xlswrite('testsheet_whole.xlsx',step_distance_SD(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'E2:E86');
    xlswrite('testsheet_whole.xlsx',step_mean_vel_average(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'F2:F86');
    xlswrite('testsheet_whole.xlsx',step_peak_vel_average(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'G2:G86');
    xlswrite('testsheet_whole.xlsx',step_mean_acc_average(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'H2:H86');
    xlswrite('testsheet_whole.xlsx',step_peak_acc_average(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'I2:I86');
    xlswrite('testsheet_whole.xlsx',step_mean_jrk_average(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'J2:J86');
    xlswrite('testsheet_whole.xlsx',step_peak_jrk_average(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'K2:K86');
    xlswrite('testsheet_whole.xlsx',overall_average_angle_back_bent(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'L2:L86');
    xlswrite('testsheet_whole.xlsx',overall_SD_angle_back_bent(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'M2:M86');
    xlswrite('testsheet_whole.xlsx',leg_average_rot_velocity(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'N2:N86');
    xlswrite('testsheet_whole.xlsx',leg_SD_rot_velocity(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'O2:O86');
    xlswrite('testsheet_whole.xlsx',vertical_impact_average(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'P2:P86');
    xlswrite('testsheet_whole.xlsx',vertical_impact_SD(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'Q2:Q86');
    xlswrite('testsheet_whole.xlsx',back_rotational_pos_range_in_sag_plane(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'R2:R86');
    xlswrite('testsheet_whole.xlsx',back_rotational_pos_mean_in_sag_plane(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'S2:S86');
    xlswrite('testsheet_whole.xlsx',back_rotational_pos_SD_in_sag_plane(1:twoMinIntCnts)',strcat('subject_',num2str(sbjct),'_whole'),'T2:T86');

end



 function [number_of_steps,step_time_average,step_time_SD,step_distance_average,... 
           step_distance_SD,step_mean_vel_average,step_peak_vel_average,step_mean_acc_average,...
           step_peak_acc_average,step_mean_jrk_average,step_peak_jrk_average,overall_average_angle_back_bent,...
           overall_SD_angle_back_bent,leg_average_rot_velocity,leg_SD_rot_velocity,vertical_impact_average,...
           vertical_impact_SD,back_rotational_pos_range_in_sag_plane,back_rotational_pos_mean_in_sag_plane,...
           back_rotational_pos_SD_in_sag_plane] = ...
               Kinematics_feature(imu_data,M_i_k_in,S_s_in,S_e_in,samplingrate,sheet_name_foot,index_value) 
    
        S_e = S_e_in(S_e_in<=size(imu_data,1));
        S_s = S_s_in(1:length(S_e));
        M_i_k = M_i_k_in(1:length(S_e),:);

        % get time vector
        data_portion = 1:floor(samplingrate*120);
        stepCnt = size(M_i_k,1);
        time_length = 1:1:size(imu_data(data_portion),1);
        time = time_length.*(1/samplingrate);
        
        imu_data_foot = imu_data(:,1:6);
        imu_data_back = imu_data(:,7:12);

        %% Filter raw data 
        acc_g = 9.55;
        [state_foot, R_foot] = IMU_Kalman(imu_data_foot(data_portion,:),samplingrate);
        [state_back, R_back] = IMU_Kalman(imu_data_back(data_portion,:),samplingrate);
        %reorganize state for ease of use

        accel_angle_foot = state_foot(:,1:3);
        gyro_angle_foot = state_foot(:,4:6);
        filtered_angle_foot = state_foot(:,7:9);
        unbiased_gyro_foot = state_foot(:,10:12);

        accel_angle_back = state_back(:,1:3);
        gyro_angle_back = state_back(:,4:6);
        filtered_angle_back = state_back(:,7:9);
        unbiased_gyro_back = state_back(:,10:12);
        
        %%
        
        acc_glob_foot = Quaternion(imu_data_foot(data_portion,1:3),R_foot);
        acc_lin_foot = acc_glob_foot + acc_g; 

        % acceleration filter parameters
        n=4;
        fc=4;
        Fs=35;
        Wn = (2/Fs)*fc;
        [b,a]=butter(n,Wn,'low');
        
        Acceleration_filt_foot = filter(b,a,acc_lin_foot);
        Acceleration_filt_magn_foot = filter(b,a,sqrt(sum(Acceleration_filt_foot.^2,2)));
        %Acceleration_filt_magn_foot = Acceleration_filt_magn_foot + min(Acceleration_filt_magn_foot);
              
        %%
        %**************************************************************
        % Features
        %**************************************************************
        %% Number of steps taken per 2 min
        
        number_of_steps = size(M_i_k,1);
        %% Step Features
        
        if size(M_i_k,1)>1 || (size(M_i_k,1)==1 && size(M_i_k,2)>1)
            for i = 1:size(M_i_k,1)
                StepLength_Time(i) = range(M_i_k{i,17});
                StepLength_Dist(i) = range(M_i_k{i,1});
                StepVelocity_avg(i) = mean(M_i_k{i,4});
                StepVelocity_peak(i) = max(M_i_k{i,4});
                StepAcceleration_avg(i) = mean(M_i_k{i,6});
                StepAcceleration_peak(i) = max(M_i_k{i,6});
                StepJerk_avg(i) = mean(M_i_k{i,8});
                StepJerk_peak(i) = max(M_i_k{i,8});
            end
            step_time_average = mean(StepLength_Time);
            step_time_SD = std(StepLength_Time);
            step_distance_average = mean(StepLength_Dist);
            step_distance_SD = std(StepLength_Dist);
            step_mean_vel_average = mean(StepVelocity_avg);
            step_peak_vel_average = mean(StepVelocity_peak);
            step_mean_acc_average = mean(StepAcceleration_avg);
            step_peak_acc_average = mean(StepAcceleration_peak);
            step_mean_jrk_average = mean(StepJerk_avg);
            step_peak_jrk_average = mean(StepJerk_peak);
        else
            step_time_average = 0;
            step_time_SD = 0;
            step_distance_average = 0;
            step_distance_SD = 0;
            step_mean_vel_average = 0;
            step_peak_vel_average = 0;
            step_mean_acc_average = 0;
            step_peak_acc_average = 0;
            step_mean_jrk_average = 0;
            step_peak_jrk_average = 0;
        end
        
        %% Angle of back during bent position
        % set threshold
        bent_threshold = 50;
        % Check if there are ANY values below this threshold
        check_threshold = filtered_angle_back(:,2) < bent_threshold;
        if 1 && all(check_threshold == 0)  
            TF_back = 0;
            disp(strcat(sheet_name_foot,'segment',num2str(index_value),'_no bent data'))
        else
            TF_back = 1;
        end
        
        if TF_back == 0 
            overall_average_angle_back_bent = 0;
            overall_SD_angle_back_bent = 0;
        else
        % get logical value corresponding to if acceleration magnitude is below
        % threshold
        below_threshold = filtered_angle_back(:,2) < bent_threshold;
        % convert from boolean to double
        below_threshold_back = +below_threshold;    
        % find indices
        below_threshold_indicies = find(below_threshold_back);
        % get overall average angle without excluding periods of less than 1
        % second
        for i = 1:length(below_threshold_indicies)
            index = below_threshold_indicies(i);
            angle_y_at_indicies(i) = filtered_angle_back(index,2);
        end

        mean_back_angle_y_unfiltered = mean(angle_y_at_indicies);

        %get average angles excluding time periods of less than 1 second

        %splits into subarrays of consecutive indices
        outputs_back_angle_indicies = SplitVec(below_threshold_indicies,'consecutive');

        k = 1;

        for i = 1:length(outputs_back_angle_indicies)
            matrix_back = cell2mat(outputs_back_angle_indicies(i));
            index = matrix_back;
            % if there are around the sampling rate consecutive indices, the bent period is longer than 1
            % second
            if length(matrix_back) > round(samplingrate)
                for j = 1:length(index)
                    % get the angles during that period
                    angle_y_at_index(j) = filtered_angle_back(index(j),2);
                end
                % average the angles during that period
                average_back_y_filtered(k) = mean(angle_y_at_index);
                k = k+1;
            else
                average_back_y_filtered(k) = 0;
            end
        end
        
            %average the mean for each period to get one value for each axis using the
            %above 1 second periods
            overall_average_angle_back_bent = mean(average_back_y_filtered);
            overall_SD_angle_back_bent = std(average_back_y_filtered);
        end
      
        %% Leg rotational velocity in Sagittal Plane
        %use data that was corrected by Kalman filter
        %segment data by each step
        if size(M_i_k,1)>1 || (size(M_i_k,1)==1 && size(M_i_k,2)>1)
            for i = 1:size(M_i_k,1)
                foot_rot_vel_segmented(i,1) = range(M_i_k{i,16});
            end

            %get mean values
            leg_average_rot_velocity = mean(foot_rot_vel_segmented);
            leg_SD_rot_velocity = std(foot_rot_vel_segmented);
        else
            leg_average_rot_velocity = 0;
            leg_SD_rot_velocity = 0;            
        end

      
        %% Vertical impact on trunk
        
        descending_accel_mag = sort(Acceleration_filt_magn_foot,'descend');
        % to rule out large random peaks
        max_accel_mag = descending_accel_mag(1);
        % define peaks relative to maximum value
        [Peak_points,Peak_mag] = peakfinder(Acceleration_filt_magn_foot,.7,max_accel_mag*.8);
        impact_acc = Peak_points(Peak_mag>0.7*(mean(Peak_mag)));
        
        vertical_impact_average = mean(Acceleration_filt_magn_foot(impact_acc));
        vertical_impact_SD = std(Acceleration_filt_magn_foot(impact_acc));

        %% Trunk Forward/Backward range of motion
        % assume that data needs to be recentered on zero
        % get average
        average_value = mean(filtered_angle_back(:,3));
        filtered_angle_back_z = filtered_angle_back(:,3) +(-average_value);

        % segment data by each step
         if size(S_s)>=1 
            for i = 1:length(S_s)
                if S_s(i)>=1
                    back_rom(i) = range(filtered_angle_back_z(S_s(i):S_e(i)));
                    back_avg(i) = mean(filtered_angle_back_z(S_s(i):S_e(i)));
                end
            end
            %get mean value
            back_rotational_pos_range_in_sag_plane = mean(back_rom);
            back_rotational_pos_mean_in_sag_plane = mean(back_avg);
            back_rotational_pos_SD_in_sag_plane = std(back_avg);
        else
            back_rotational_pos_range_in_sag_plane = 0;
            back_rotational_pos_mean_in_sag_plane = 0;
            back_rotational_pos_SD_in_sag_plane = 0;  
        end
     
 end
 
 
function [State, R] = IMU_Kalman(imu_data, SamplingRate)

duration = (length(imu_data)-1)/SamplingRate;
dt = 1/SamplingRate;
t = (0:dt:duration)';

acc_data = imu_data(:,1:3);
gyro_data = imu_data(:,4:6);
gyro = gyro_data;

v_gyro = sqrt(0.3)*randn(length(t),3); % measurement noise for gyro,  variance = 0.3
v_acc = sqrt(0.3)*randn(length(t),3); % measurement noise for accellerometer, variance = 0.3

% Compute the angles computed by using only accelerometers of gyroscope
x = [-atan(acc_data(:,1)./sqrt(sum([acc_data(:,2).^2 acc_data(:,3).^2],2)))*(180/pi), ...
     atan(acc_data(:,2)./sqrt(sum([acc_data(:,1).^2 acc_data(:,3).^2],2)))*(180/pi), ...
     atan(acc_data(:,3)./sqrt(sum([acc_data(:,1).^2 acc_data(:,2).^2],2)))*(180/pi)];
x_acc = x + v_acc; % Angle computed from accelerometer measurement

% figure;
% plot(x(:,1));
% hold on
% plot(x(:,2));
% hold on
% plot(x(:,3));
% hold on

x_gyro = cumsum(gyro*dt,1); % Angle computed by integration of gyro measurement

% figure;
% plot(x_gyro(:,1));
% hold on
% plot(x_gyro(:,2));
% hold on
% plot(x_gyro(:,3));
% hold on


P = [0.1 0; 0 0.1];  
R_angle = 0.5;   

Q_angle = 0.9;
Q_gyro = 0.1;
Q = [Q_angle 0; 0 Q_gyro];

A = [0 1; 0 1];
q_bias = [0 0 0]; % Initialize gyro bias
angle = [0 90 0]; % Initialize gyro angle
q_m = 0;
X = [0 90 0; 0 0 0];

    for i=1:length(t)

         % Gyro update 

         q_m = gyro(i,:);

         q = q_m - q_bias; % gyro bias removal

         Pdot = A*P + P*A' + Q;

         rate = q;

         angle = angle + q*dt;

         P = P + Pdot*dt;

         % Kalman (Accelerometer) update 

         H = [0 1]; 
         angle_err = x_acc(i,:)-angle;
         E = H*P*H' + R_angle;

         K = P*H'*inv(E);

         P = P - K*H*P;
         X = X + K * angle_err;
         
         
         x1(i,:) = X(1,:);
         R(:,:,i) = ang2orth(x1(i,:));
         x2(i,:) = X(2,:);
         angle = x1(i,:);
         q_bias = x2(i,:);

         x3(i,:) = q;  % unbiased gyro rate
    end

State = [x_acc x_gyro x1 x3];

x_acc = State(:,1:3);   % angle calculated from accelerometer
x_gyro = State(:,4:6);  % angle calculated from gyro 
filtered_angle = State(:,7:9);      % filtered angle
unbiased_gyro = State(:,10:12);     % unbiased gyro
rotation_mat = R;

% figure;
% plot(filtered_angle(:,1));
% hold on
% plot(filtered_angle(:,2));
% hold on
% plot(filtered_angle(:,3));
% hold on

    function orthm = ang2orth(ang) 
        sa = sind(ang(1)); ca = cosd(ang(1)); 
        sb = sind(ang(2)-80); cb = cosd(ang(2)-80); 
        sc = sind(ang(3)); cc = cosd(ang(3)); 

        ra = [  ca,  sa,  0; ... 
               -sa,  ca,  0; ... 
                 0,   0,  1]; 
        rb = [  cb,  0, -sb; ... 
                 0,  1,  0; ... 
                sb,  0,  cb]; 
        rc = [  1,   0,   0; ... 
                0,   cc, sc;... 
                0,  -sc, cc]; 
        orthm = (ra*rb*rc)';
    end

end


        

function X_glob = Quaternion(X_body,R)

Q = dcm2q(R);
X_glob = qvqc(Q,X_body);

%% Functions


    function q=dcm2q(R)
        % DCM2Q(R) converts direction cosine matrices into quaternions.
        %
        %     The resultant quaternion(s) will perform the equivalent vector
        %     transformation as the input DCM(s), i.e.:
        %
        %       qconj(Q)*V*Q = R*V
        %
        %     where R is the DCM, V is a vector, and Q is the quaternion.  Note that
        %     for purposes of quaternion-vector multiplication, a vector is treated
        %     as a quaterion with a scalar element of zero.
        %
        %     If the input is a 3x3xN array, the output will be a vector of
        %     quaternions where input direction cosine matrix R(:,:,k) corresponds
        %     to the output quaternion Q(k,:).
        %
        %     Note that this function is meaningless for non-orthonormal matrices!
        %
        % See also Q2DCM.

        % Release: $Name: quaternions-1_3 $
        % $Revision: 1.11 $
        % $Date: 2009-07-25 04:28:18 $

        % Copyright (c) 2000-2009, Jay A. St. Pierre.  All rights reserved.

        % Thanks to Tatsuki Kashitani for pointing out the numerical instability in
        % the original implementation.  His suggested fix also included a check for
        % the "sr" values below being zero.  But I think I've convinced myself that
        % this isn't necessary if the input matrices are orthonormal (or at least
        % very close to orthonormal).

        if nargin~=1
          error('One input argument required');
        else
          size_R=size(R);
          if ( size_R(1)~=3 || size_R(2)~=3 || length(size_R)>3 )
            error('Invalid input: must be a 3x3xN array')
          end
        end

        q = zeros( 4, size( R, 3 ) );

        for id_dcm = 1 : size( R, 3 )
          dcm = R( :, :, id_dcm );
          if trace( dcm ) > 0
            % Positve Trace Algorithm
            sr  = sqrt( 1 + trace( dcm ));
            sr2 = 2*sr;
            q(1,id_dcm) = ( dcm(2,3) - dcm(3,2) ) / sr2;
            q(2,id_dcm) = ( dcm(3,1) - dcm(1,3) ) / sr2;
            q(3,id_dcm) = ( dcm(1,2) - dcm(2,1) ) / sr2;
            q(4,id_dcm) = 0.5 * sr;
          else
            % Negative Trace Algorithm
            if ( dcm(1,1) > dcm(2,2) ) && ( dcm(1,1) > dcm(3,3) )
              % Maximum Value at DCM(1,1)
              sr  = sqrt( 1 + (dcm(1,1) - ( dcm(2,2) + dcm(3,3) )) );
              sr2 = 2*sr;
              q(1,id_dcm) = 0.5 * sr;
              q(2,id_dcm) = ( dcm(2,1) + dcm(1,2) ) / sr2;
              q(3,id_dcm) = ( dcm(3,1) + dcm(1,3) ) / sr2;
              q(4,id_dcm) = ( dcm(2,3) - dcm(3,2) ) / sr2;
            elseif dcm(2,2) > dcm(3,3)
              % Maximum Value at DCM(2,2)
              sr  = sqrt( 1 + (dcm(2,2) - ( dcm(3,3) + dcm(1,1) )) );
              sr2 = 2*sr;
              q(1,id_dcm) = ( dcm(2,1) + dcm(1,2) ) / sr2;
              q(2,id_dcm) = 0.5 * sr;
              q(3,id_dcm) = ( dcm(2,3) + dcm(3,2) ) / sr2;
              q(4,id_dcm) = ( dcm(3,1) - dcm(1,3) ) / sr2;
            else
              % Maximum Value at DCM(3,3)
              sr  = sqrt( 1 + (dcm(3,3) - ( dcm(1,1) + dcm(2,2) )) );
              sr2 = 2*sr;
              q(1,id_dcm) = ( dcm(3,1) + dcm(1,3) ) / sr2;
              q(2,id_dcm) = ( dcm(2,3) + dcm(3,2) ) / sr2;
              q(3,id_dcm) = 0.5 * sr;
              q(4,id_dcm) = ( dcm(1,2) - dcm(2,1) ) / sr2;
            end
          end
        end

        % Make quaternion vector a column of quaternions
        q=q.';

        q=real(q);
    end


    function qtype=isq(q)
        % ISQ(Q) checks to see if Q is a quaternion or set of quaternions.
        %     ISQ returns a value accordingly:
        %
        %        0 if Q is not a quaternion or vector of quaternions:
        %          has more than 2 dimensions or neither dimension is of length 4
        %       
        %        1 if the component quaternions of Q are column vectors:
        %          Q is 4xN, where N~=4, or
        %          Q is 4x4 and only the columns are normalized 
        %
        %        2 if the component quaternions of Q are row vectors:
        %          Q is Nx4, where N~=4, or
        %          Q is 4x4 and only the rows are normalized 
        %
        %        3 if the shape of the component quaternions is indeterminant:
        %          Q is 4x4, and either both the columns and rows are normalized
        %          or neither the columns nor rows are normalized.
        %
        %     In other words, if Q is 4x4, ISQ attempts to discern the shape of
        %     component quaternions by determining whether the rows or the columns
        %     are normalized (i.e., it assumes that normalized quaternions are
        %     the more typical use of quaternions).
        %
        %     The test for normalization uses 2*EPS as a tolerance.
        %
        % See also ISNORMQ, EPS.

        % Release: $Name: quaternions-1_3 $
        % $Revision: 1.7 $
        % $Date: 2009-07-26 20:05:12 $

        % Copyright (c) 2001-2009, Jay A. St. Pierre.  All rights reserved.

        if nargin~=1

          error('isq() requires one input argument');

        else

          tol=2*eps;

          size_q=size(q);

          if ( length(size_q)~=2 || max(size_q==4)~=1 )
            qtype=0; % Not a quaternion or quaternion vector

          elseif ( size_q(1)==4 && ...
                   ( size_q(2)~=4 || ( ~sum((sum(q.^2,1)-ones(1,4))>tol) &&   ...
                                        sum((sum(q.^2,2)-ones(4,1))>tol)    ) ...
                     ) ...
                   )
            qtype=1; % Component q's are column vectors

          elseif ( size_q(2)==4 && ...
                   ( size_q(1)~=4 || ( ~sum((sum(q.^2,2)-ones(4,1))>tol) &&   ...
                                        sum((sum(q.^2,1)-ones(1,4))>tol)    ) ...
                     ) ...
                   )
            qtype=2; % Component q's are row vectors

          else
            qtype=3; % Component q's are either columns or rows (indeterminate)

          end

        end
    end


    function qout=qconj(qin)
        % QCONJ(Q) calculates the conjugate of the quaternion Q.
        %     Works on "vectors" of quaterions as well.  Will return the same shape
        %     vector as input.  If input is a vector of four quaternions, QCONJ will
        %     determine whether the quaternions are row or column vectors according
        %     to ISQ.
        %
        % See also ISQ.

        % Release: $Name: quaternions-1_3 $
        % $Revision: 1.16 $
        % $Date: 2009-07-26 20:05:12 $

        % Copyright (c) 2001-2009, Jay A. St. Pierre.  All rights reserved.


        if nargin~=1
          error('qconj() requires one input argument');
        else
          qtype = isq(qin);
          if ( qtype==0 )
            error(...
              'Invalid input: must be a quaternion or a vector of quaternions')
          elseif ( qtype==3 )
            warning(...
              'qconj:indeterminateShape', ...
              'Component quaternion shape indeterminate, assuming row vectors')
          end
        end

        % Make sure component quaternions are row vectors
        if( qtype == 1 )
          qin=qin.';
        end

        qout(:,1)=-qin(:,1);
        qout(:,2)=-qin(:,2);
        qout(:,3)=-qin(:,3);
        qout(:,4)= qin(:,4);

        % Make sure output is same shape as input
        if( qtype == 1 )
          qout=qout.';
        end
    end

    function v_out=qcvq(q,v)
        % QcVQ(Q,V) performs the operation qconj(Q)*V*Q
        %     where the vector is treated as a quaternion with a scalar element of
        %     zero.
        %
        %     Q and V can be vectors of quaternions and vectors, but they must
        %     either be the same length or one of them must have a length of one.
        %     The output will have the same shape as V.  Q will be passed through
        %     QNORM to ensure it is normalized.
        %
        % See also QVQc, QNORM, QMULT.

        % Note that QNORM is invoked by QMULT, therefore QcQV does not invoke
        % it directly.

        % Release: $Name: quaternions-1_3 $
        % $Revision: 1.2 $
        % $Date: 2009-07-26 20:05:12 $

        % Copyright (c) 2000-2009, Jay A. St. Pierre.  All rights reserved.


        if nargin~=2
          error('Two input arguments required.');
        else

          qtype=isq(q);
          if ( qtype == 0 )
            error('Input Q must be a quaternion or a vector of quaternions')
          end

          size_v=size(v);
          if ( length(size_v)~=2 || max(size_v==3)~=1 )
            error(['Invalid input: second input must be a 3-element vector', 10, ...
                   'or a vector of 3-element vectors'])
          end

        end

        % Make sure q is a column of quaternions
        if ( qtype==1 )
          q=q.';
        end

        % Make sure v is a column of vectors
        row_of_vectors = (size_v(2) ~= 3);
        if ( row_of_vectors )
          v=v.';
          size_v=size(v);
        end

        size_q=size(q);

        if (  size_q(1)~=size_v(1) && size_q(1)~=1 && size_v(1)~=1 )
          error(['Inputs do not have the same number of elements:', 10, ...
                 '   number of quaternions in q = ', num2str(size_q(1)), 10,...
                 '   number of vectors in v     = ', num2str(size_v(1)), 10,...
                 'Inputs must have the same number of elements, or', 10, ...
                 'one of the inputs must have a single element.']) 
        elseif ( size_q(1)==1 && size_v(1)==3 )
          if ( qtype==1 )
            warning(...
              'qcvq:assumingVcols', ...
              'Q is 4x1 and V is 3x3: assuming vectors are column vectors')
            row_of_vectors = 1;
            v=v.';
          else
            warning(...
              'qcvq:assumingVrows', ...
              'Q is 1x4 and V is 3x3: assuming vectors are row vectors')
          end
        elseif ( qtype==3 && size_v(1)==1 )
          if ( row_of_vectors )
            warning(...
              'qcvq:assumingQcols', ...
              'Q is 4x4 and V is 3x1: assuming quaternions are column vectors')
            q=q.';
          else
            warning(...
              'qcvq:assumingQrows', ...
              'Q is 4x4 and V is 1x3: assuming quaternions are row vectors')
          end  
        end

        % Build up full vectors if one input is a singleton
        if ( size_q(1) ~= size_v(1) )
          ones_length = ones(max(size_q(1),size_v(1)),1);
          if ( size_q(1) == 1 )
            q = [q(1)*ones_length ...
                 q(2)*ones_length ...
                 q(3)*ones_length ...
                 q(4)*ones_length ];
          else % size_v(1) == 1
            v = [v(1)*ones_length ...
                 v(2)*ones_length ...
                 v(3)*ones_length ];    
          end
        end

        % Add an element to V
        v(:,4)=zeros(size_v(1),1);

        % Turn off warnings before calling qconj (it has simillar warnings as
        % qvxform, so all warnings would just be duplicated).  Save current state of
        % warnings, though.
        warning_state = warning; warning('off', 'qconj:indeterminateShape');
        local_warning = lastwarn;

        % Perform transform
        vt=qmult(qconj(q),qmult(v,q));

        % Restore warning state to original state
        warning(warning_state);
        lastwarn(local_warning);

        % Eliminate last element of vt for output
        v_out = vt(:,1:3);

        % Make sure output vectors are the same shape as input vectors
        if ( row_of_vectors )
          v_out = v_out.';
        end
    end
    
    
    function q_out=qmult(q1,q2)
        % QMULT(Q1,Q2) calculates the product of two quaternions Q1 and Q2.
        %    Inputs can be vectors of quaternions, but they must either have the
        %    same number of component quaternions, or one input must be a single
        %    quaternion.  QMULT will determine whether the component quaternions of
        %    the inputs are row or column vectors according to ISQ.
        %  
        %    The output will have the same shape as Q1.  If the component
        %    quaternions of either Q1 or Q2 (but not both) are of indeterminate
        %    shape (see ISQ), then the shapes will be assumed to be the same for
        %    both inputs.  If both Q1 and Q2 are of indeterminate shape, then both
        %    are assumed to be composed of row vector quaternions.
        %
        % See also ISQ.

        % Release: $Name: quaternions-1_3 $
        % $Revision: 1.14 $
        % $Date: 2009-07-26 20:05:12 $

        % Copyright (c) 2001-2009, Jay A. St. Pierre.  All rights reserved.


        if nargin~=2
          error('qmult() requires two input arguments');
        else
          q1type = isq(q1);
          if ( q1type == 0 )
            error(['Invalid input: q1 must be a quaternion or a vector of' ...
                  ' quaternions'])
          end
          q2type = isq(q2);
          if ( q2type == 0 )
            error(['Invalid input: q2 must be a quaternion or a vector of' ...
                  ' quaternions'])
          end
        end

        % Make sure q1 is a column of quaternions (components are rows)
        if ( q1type==1 || (q1type==3 && q2type==1) )
          q1=q1.';
        end

        % Make sure q2 is a column of quaternions (components are rows)
        if ( q2type==1 || (q2type==3 && q1type==1) )
          q2=q2.';
        end

        num_q1=size(q1,1);
        num_q2=size(q2,1);

        if (  num_q1~=num_q2 && num_q1~=1 && num_q2~=1 )
          error(['Inputs do not have the same number of elements:', 10, ...
                 '   number of quaternions in q1 = ', num2str(num_q1), 10,...
                 '   number of quaternions in q2 = ', num2str(num_q2), 10,...
                 'Inputs must have the same number of elements, or', 10, ...
                 'one of the inputs must be a single quaternion (not a', 10, ...
                 'vector of quaternions).']) 
        end

        % Build up full quaternion vector if one input is a single quaternion
        if ( num_q1 ~= num_q2 )
          ones_length = ones(max(num_q1,num_q2),1);
          if ( num_q1 == 1 )
            q1 = [q1(1)*ones_length ...
                  q1(2)*ones_length ...
                  q1(3)*ones_length ...
                  q1(4)*ones_length ];
          else % num_q2 == 1
            q2 = [q2(1)*ones_length ...
                  q2(2)*ones_length ...
                  q2(3)*ones_length ...
                  q2(4)*ones_length ];    
          end
        end

        % Products

        % If q1 and q2 are not vectors of quaternions, then:
        %
        %   q1*q2 = q1*[ q2(4) -q2(3)  q2(2) -q2(1)
        %                q2(3)  q2(4) -q2(1) -q2(2)
        %               -q2(2)  q2(1)  q2(4) -q2(3)
        %                q2(1)  q2(2)  q2(3)  q2(4) ]
        %
        % But to deal with vectorized quaternions, we have to use the ugly
        % commands below.

        prod1 = ...
            [ q1(:,1).*q2(:,4) -q1(:,1).*q2(:,3)  q1(:,1).*q2(:,2) -q1(:,1).*q2(:,1)];
        prod2 = ...
            [ q1(:,2).*q2(:,3)  q1(:,2).*q2(:,4) -q1(:,2).*q2(:,1) -q1(:,2).*q2(:,2)];
        prod3 = ...
            [-q1(:,3).*q2(:,2)  q1(:,3).*q2(:,1)  q1(:,3).*q2(:,4) -q1(:,3).*q2(:,3)];
        prod4 = ...
            [ q1(:,4).*q2(:,1)  q1(:,4).*q2(:,2)  q1(:,4).*q2(:,3)  q1(:,4).*q2(:,4)];

        q_out = prod1 + prod2 + prod3 + prod4;

        % Make sure output is same format as q1
        if ( q1type==1 || (q1type==3 && q2type==1) )
          q_out=q_out.';
        end

        % NOTE that the following algorithm proved to be slower than the one used
        % above:
        %
        % q_out = zeros(size(q1));
        % 
        % q_out(:,1:3) = ...
        %     [q1(:,4) q1(:,4) q1(:,4)].*q2(:,1:3) + ...
        %     [q2(:,4) q2(:,4) q2(:,4)].*q1(:,1:3) + ...
        %     cross(q1(:,1:3), q2(:,1:3));
        % 
        % q_out(:,4) = q1(:,4).*q2(:,4) - dot(q1(:,1:3), q2(:,1:3), 2);
    end

    
    function v_out=qvqc(q,v)
        % QVQc(Q,V) performs the operation Q*V*qconj(Q)
        %     where the vector is treated as a quaternion with a scalar element of
        %     zero. 
        %
        %     Q and V can be vectors of quaternions and vectors, but they must
        %     either be the same length or one of them must have a length of one.
        %     The output will have the same shape as V.  Q will be passed through
        %     QNORM to ensure it is normalized.
        %
        % See also QcQV, QNORM

        % Release: $Name: quaternions-1_3 $
        % $Revision: 1.1 $
        % $Date: 2009-07-24 19:14:44 $

        % Copyright (c) 2000-2009, Jay A. St. Pierre.  All rights reserved.


        if nargin~=2
          error('Two input arguments required');
        else
          q     = qconj(q);
          v_out = qcvq(q, v);
        end
    end
    
end


function varargout = peakfinder(x0, sel, thresh, extrema, includeEndpoints, interpolate)
    %PEAKFINDER Noise tolerant fast peak finding algorithm
    %   INPUTS:
    %       x0 - A real vector from the maxima will be found (required)
    %       sel - The amount above surrounding data for a peak to be,
    %           identified (default = (max(x0)-min(x0))/4). Larger values mean
    %           the algorithm is more selective in finding peaks.
    %       thresh - A threshold value which peaks must be larger than to be
    %           maxima or smaller than to be minima.
    %       extrema - 1 if maxima are desired, -1 if minima are desired
    %           (default = maxima, 1)
    %       includeEndpoints - If true the endpoints will be included as
    %           possible extrema otherwise they will not be included
    %           (default = true)
    %       interpolate - If true quadratic interpolation will be performed
    %           around each extrema to estimate the magnitude and the
    %           position of the peak in terms of fractional indicies. Note that
    %           unlike the rest of this function interpolation assumes the
    %           input is equally spaced. To recover the x_values of the input
    %           rather than the fractional indicies you can do:
    %           peakX = x0 + (peakLoc - 1) * dx
    %           where x0 is the first x value and dx is the spacing of the
    %           vector. Output peakMag to recover interpolated magnitudes.
    %           See example 2 for more information.
    %           (default = false)
    %
    %   OUTPUTS:
    %       peakLoc - The indicies of the identified peaks in x0
    %       peakMag - The magnitude of the identified peaks
    %
    %   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
    %       are at least 1/4 the range of the data above surrounding data.
    %
    %   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
    %       that are at least sel above surrounding data.
    %
    %   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local
    %       maxima that are at least sel above surrounding data and larger
    %       (smaller) than thresh if you are finding maxima (minima).
    %
    %   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
    %       data if extrema > 0 and the minima of the data if extrema < 0
    %
    %   [peakLoc] = peakfinder(x0,sel,thresh,extrema, includeEndpoints)
    %       returns the endpoints as possible extrema if includeEndpoints is
    %       considered true in a boolean sense
    %
    %   [peakLoc, peakMag] = peakfinder(x0,sel,thresh,extrema,interpolate)
    %       returns the results of results of quadratic interpolate around each
    %       extrema if interpolate is considered to be true in a boolean sense
    %
    %   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
    %       local maxima as well as the magnitudes of those maxima
    %
    %   If called with no output the identified maxima will be plotted along
    %       with the input data.
    %
    %   Note: If repeated values are found the first is identified as the peak
    %
    % Example 1:
    % t = 0:.0001:10;
    % x = 12*sin(10*2*pi*t)-3*sin(.1*2*pi*t)+randn(1,numel(t));
    % x(1250:1255) = max(x);
    % peakfinder(x)
    %
    % Example 2:
    % ds = 100;  % Downsample factor
    % dt = .001; % Time step
    % ds_dt = ds*dt; % Time delta after downsampling
    % t0 = 1;
    % t = t0:dt:5 + t0;
    % x = 0.2-sin(0.01*2*pi*t)+3*cos(7/13*2*pi*t+.1)-2*cos((1+pi/10)*2*pi*t+0.2)-0.2*t;
    % x(end) = min(x);
    % x_ds = x(1:ds:end); % Downsample to test interpolation
    % [minLoc, minMag] = peakfinder(x_ds, .8, 0, -1, false, true);
    % minT = t0 + (minLoc - 1) * ds_dt; % Take into account 1 based indexing
    % p = plot(t,x,'-',t(1:ds:end),x_ds,'o',minT,minMag,'rv');
    % set(p(2:end), 'linewidth', 2); % Show the markers more clearly
    % legend('Actual Data', 'Input Data', 'Estimated Peaks');
    % Copyright Nathanael C. Yoder 2015 (nyoder@gmail.com)

    % Perform error checking and set defaults if not passed in
    narginchk(1, 6);
    nargoutchk(0, 2);

    s = size(x0);
    flipData =  s(1) < s(2);
    len0 = numel(x0);
    if len0 ~= s(1) && len0 ~= s(2)
        error('PEAKFINDER:Input','The input data must be a vector')
    elseif isempty(x0)
        varargout = {[],[]};
        return;
    end
    if ~isreal(x0)
        warning('PEAKFINDER:NotReal','Absolute value of data will be used')
        x0 = abs(x0);
    end

    if nargin < 2 || isempty(sel)
        sel = (max(x0)-min(x0))/4;
    elseif ~isnumeric(sel) || ~isreal(sel)
        sel = (max(x0)-min(x0))/4;
        warning('PEAKFINDER:InvalidSel',...
            'The selectivity must be a real scalar.  A selectivity of %.4g will be used',sel)
    elseif numel(sel) > 1
        warning('PEAKFINDER:InvalidSel',...
            'The selectivity must be a scalar.  The first selectivity value in the vector will be used.')
        sel = sel(1);
    end

    if nargin < 3 || isempty(thresh)
        thresh = [];
    elseif ~isnumeric(thresh) || ~isreal(thresh)
        thresh = [];
        warning('PEAKFINDER:InvalidThreshold',...
            'The threshold must be a real scalar. No threshold will be used.')
    elseif numel(thresh) > 1
        thresh = thresh(1);
        warning('PEAKFINDER:InvalidThreshold',...
            'The threshold must be a scalar.  The first threshold value in the vector will be used.')
    end

    if nargin < 4 || isempty(extrema)
        extrema = 1;
    else
        extrema = sign(extrema(1)); % Should only be 1 or -1 but make sure
        if extrema == 0
            error('PEAKFINDER:ZeroMaxima','Either 1 (for maxima) or -1 (for minima) must be input for extrema');
        end
    end

    if nargin < 5 || isempty(includeEndpoints)
        includeEndpoints = true;
    end

    if nargin < 6 || isempty(interpolate)
        interpolate = false;
    end

    x0 = extrema*x0(:); % Make it so we are finding maxima regardless
    thresh = thresh*extrema; % Adjust threshold according to extrema.
    dx0 = diff(x0); % Find derivative
    dx0(dx0 == 0) = -eps; % This is so we find the first of repeated values
    ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1; % Find where the derivative changes sign

    % Include endpoints in potential peaks and valleys as desired
    if includeEndpoints
        x = [x0(1);x0(ind);x0(end)];
        ind = [1;ind;len0];
        minMag = min(x);
        leftMin = minMag;
    else
        x = x0(ind);
        minMag = min(x);
        leftMin = min(x(1), x0(1));
    end

    % x only has the peaks, valleys, and possibly endpoints
    len = numel(x);

    if len > 2 % Function with peaks and valleys
        % Set initial parameters for loop
        tempMag = minMag;
        foundPeak = false;

        if includeEndpoints
            % Deal with first point a little differently since tacked it on
            % Calculate the sign of the derivative since we tacked the first
            %  point on it does not neccessarily alternate like the rest.
            signDx = sign(diff(x(1:3)));
            if signDx(1) <= 0 % The first point is larger or equal to the second
                if signDx(1) == signDx(2) % Want alternating signs
                    x(2) = [];
                    ind(2) = [];
                    len = len-1;
                end
            else % First point is smaller than the second
                if signDx(1) == signDx(2) % Want alternating signs
                    x(1) = [];
                    ind(1) = [];
                    len = len-1;
                end
            end
        end

        % Skip the first point if it is smaller so we always start on a
        %   maxima
        if x(1) >= x(2)
            ii = 0;
        else
            ii = 1;
        end

        % Preallocate max number of maxima
        maxPeaks = ceil(len/2);
        peakLoc = zeros(maxPeaks,1);
        peakMag = zeros(maxPeaks,1);
        cInd = 1;
        % Loop through extrema which should be peaks and then valleys
        while ii < len
            ii = ii+1; % This is a peak
            % Reset peak finding if we had a peak and the next peak is bigger
            %   than the last or the left min was small enough to reset.
            if foundPeak
                tempMag = minMag;
                foundPeak = false;
            end

            % Found new peak that was lager than temp mag and selectivity larger
            %   than the minimum to its left.
            if x(ii) > tempMag && x(ii) > leftMin + sel
                tempLoc = ii;
                tempMag = x(ii);
            end

            % Make sure we don't iterate past the length of our vector
            if ii == len
                break; % We assign the last point differently out of the loop
            end

            ii = ii+1; % Move onto the valley
            % Come down at least sel from peak
            if ~foundPeak && tempMag > sel + x(ii)
                foundPeak = true; % We have found a peak
                leftMin = x(ii);
                peakLoc(cInd) = tempLoc; % Add peak to index
                peakMag(cInd) = tempMag;
                cInd = cInd+1;
            elseif x(ii) < leftMin % New left minima
                leftMin = x(ii);
            end
        end

        % Check end point
        if includeEndpoints
            if x(end) > tempMag && x(end) > leftMin + sel
                peakLoc(cInd) = len;
                peakMag(cInd) = x(end);
                cInd = cInd + 1;
            elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
                peakLoc(cInd) = tempLoc;
                peakMag(cInd) = tempMag;
                cInd = cInd + 1;
            end
        elseif ~foundPeak
            if x(end) > tempMag && x(end) > leftMin + sel
                peakLoc(cInd) = len;
                peakMag(cInd) = x(end);
                cInd = cInd + 1;
            elseif tempMag > min(x0(end), x(end)) + sel
                peakLoc(cInd) = tempLoc;
                peakMag(cInd) = tempMag;
                cInd = cInd + 1;
            end
        end

        % Create output
        if cInd > 1
            peakInds = ind(peakLoc(1:cInd-1));
            peakMags = peakMag(1:cInd-1);
        else
            peakInds = [];
            peakMags = [];
        end
    else % This is a monotone function where an endpoint is the only peak
        [peakMags,xInd] = max(x);
        if includeEndpoints && peakMags > minMag + sel
            peakInds = ind(xInd);
        else
            peakMags = [];
            peakInds = [];
        end
    end

    % Apply threshold value.  Since always finding maxima it will always be
    %   larger than the thresh.
    if ~isempty(thresh)
        m = peakMags>thresh;
        peakInds = peakInds(m);
        peakMags = peakMags(m);
    end

    if interpolate && ~isempty(peakMags)
        middleMask = (peakInds > 1) & (peakInds < len0);
        noEnds = peakInds(middleMask);

        magDiff = x0(noEnds + 1) - x0(noEnds - 1);
        magSum = x0(noEnds - 1) + x0(noEnds + 1)  - 2 * x0(noEnds);
        magRatio = magDiff ./ magSum;

        peakInds(middleMask) = peakInds(middleMask) - magRatio/2;
        peakMags(middleMask) = peakMags(middleMask) - magRatio .* magDiff/8;
    end

    % Rotate data if needed
    if flipData
        peakMags = peakMags.';
        peakInds = peakInds.';
    end

    % Change sign of data if was finding minima
    if extrema < 0
        peakMags = -peakMags;
        x0 = -x0;
    end

    % Plot if no output desired
    if nargout == 0
        if isempty(peakInds)
            disp('No significant peaks found')
        else
            figure;
            plot(1:len0,x0,'.-',peakInds,peakMags,'ro','linewidth',2);
        end
    else
        varargout = {peakInds,peakMags};
    end

end
