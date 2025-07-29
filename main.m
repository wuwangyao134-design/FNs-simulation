% =========================================================================
% NSGA-II、MyNSGA-II、MOPSO、Baseline 主执行脚本 (main_compare.m) - 最终完善版
% 描述: 对比原始 NSGA-II、带有记忆机制及自适应参数的 MyNSGA-II、
%       纯血 MOPSO 和几何距离贪婪基线。
%       新增：多场景测试与多轮统计运行，并构建近似参考前沿。
%       重要更新：每次统计运行使用固定的随机阴影衰落，以提高优化算法收敛性。
%       更新：自动保存最终统计结果到指定文件夹。
% =========================================================================
clc;
clear;
close all;

%% 0. 全局参数设定
nSlots = 5; % 仿真时隙数量
num_stat_runs = 30; % 统计运行次数：每个场景、每个算法都独立运行30次
N_RUNS_BASELINE = 500; 

% --- 定义阴影衰落的标准差 (新增，用于在 run_idx 循环中生成固定值) ---
shadow_std_dev.LoS = 3.0;  % LoS场景阴影衰落标准差 (dB)
shadow_std_dev.NLoS = 8.29; % NLoS场景阴影衰落标准差 (dB)

% --- 定义实验场景参数 ---
% 每个 cell 包含 {I, M, R2_x, R2_y}
experimental_scenarios = {
    % {20, 2, 100, 100};  % S1
    % {20, 2, 200, 200};  % S2
    % {45, 3, 100, 100};  % S3
    % {45, 3, 200, 200};  % S4
    % {80, 4, 100, 100};  % S5
    {80, 4, 200, 200};  % S6
};
num_scenarios = length(experimental_scenarios); % 在手动单场景模式下，num_scenarios将是1

% --- 初始化全局结果存储结构 (在手动单场景模式下，它将只存储当前场景的数据) ---
all_scenario_results = struct();
alg_names_for_results = {'MyNSGA_II', 'NSGA_II', 'MOPSO', 'Baseline'};
metric_names = {'IGD', 'HV', 'Spacing', 'Spread', 'Runtime', 'NumFeasibleSolutions'}; 

for a_name = alg_names_for_results
    for m_name = metric_names
        % all_scenario_results.AlgName.MetricName{scenario_idx}(run_idx, t_slot)
        all_scenario_results.(a_name{1}).(m_name{1}) = cell(num_scenarios, 1);
        for s = 1:num_scenarios
            all_scenario_results.(a_name{1}).(m_name{1}){s} = zeros(num_stat_runs, nSlots);
        end
    end
end

% --- 新增：定义结果文件和图表的输出文件夹 ---
% 你可以在这里指定你想要的路径。例如：
output_base_folder = 'E:\guthub-matlab\第二篇论文代码\leachkmeans-master'; % 指定项目根目录

% 在项目根目录下，创建一个专门存放实验结果的子目录名
results_subfolder_name = 'ExperimentResults_Output'; % 你可以改成任何你喜欢的名字

% 组合成完整的输出目录路径
output_main_folder = fullfile(output_base_folder, results_subfolder_name);

% 创建一个带有时间戳的子文件夹，以便每次运行都有独立的输出目录
output_timestamp_folder = fullfile(output_main_folder, string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));

% 检查并创建目录
if ~exist(output_timestamp_folder, 'dir')
    mkdir(output_timestamp_folder);
    fprintf('已创建输出目录: %s\n', output_timestamp_folder);
else
    fprintf('输出目录已存在: %s\n', output_timestamp_folder);
end

%% 5. 最外层循环：遍历实验场景 (Scenario Loop)
% 定义 problem_base，包含所有场景通用的参数
problem_base.objFunc = @EvaluateParticle;
problem_base.Tslot = 2.0;%根据规模可调
problem_base.systemTotalBandwidth = 225e6;
problem_base.nObj = 2; % 两个目标 (G1 和 G2)

for s_idx = 1:num_scenarios % 场景循环 (在手动单场景模式下，此循环只运行一次)
    % 每次场景循环开始时，将 problem 初始化为 problem_base，然后覆盖场景特定参数
    problem = problem_base; 
    current_scenario_config = experimental_scenarios{s_idx};
    
    % --- 根据当前场景配置 Problem 结构体中的变化参数 ---
    problem.nTerminals = current_scenario_config{1};
    problem.nFogNodes = current_scenario_config{2};
    problem.area = [0 current_scenario_config{3}; 0 current_scenario_config{4}]; % R2_x = R2_y

    fprintf('\n========== 开始处理场景 %d (I=%d, M=%d, R2=%dx%d) ==========\n', ...
            s_idx, problem.nTerminals, problem.nFogNodes, current_scenario_config{3}, current_scenario_config{4});
            
    % --- 参数配置 (为所有算法分别配置参数) ---
    % 这些参数在每个场景中是固定的，所以放在场景循环内部是合适的
    params_nsga = struct('N', 200, 'T_max', 500, 'pc', 0.9, 'pm', 0.05, 'mu', 20, 'mum', 20);
    params_my_nsga = params_nsga;
    params_my_nsga.adaptive_enabled = true; params_my_nsga.pm_min = 0.01; params_my_nsga.pm_max = 0.1;
    params_my_nsga.mum_min = 5; params_my_nsga.mum_max = 10; params_my_nsga.memory_ratio = 0.1;

    params_mopso = struct('N', 200, 'T_max', 500, 'w_max', 0.9, 'w_min', 0.4, 'c1', 1.5, 'c2', 1.5, ...
        'nGrid', 10, 'ArchiveSize', 150, 'pm_mopso', 0.1, 'mum_mopso', 10, 'eta', 2);
    % 基准方法无特殊参数 (Baseline_Greedy 函数内部逻辑决定)

    % --- 存储当前场景下所有算法在所有时隙的所有 num_stat_runs 次运行的最终前沿数据 ---
    % (在手动单场景模式下，这个变量只存储当前场景的数据)
    current_scenario_all_run_final_fronts = cell(length(alg_names_for_results), nSlots, num_stat_runs);

    % --- 新增：用于构建全局 PF* 的临时存储 ---
    % (在手动单场景模式下，这个变量只存储当前场景的数据，并在场景结束后被清理)
    all_combined_solutions_for_global_pf = cell(num_stat_runs, nSlots);

%% 6. 中层循环：统计运行 (Statistical Run Loop)
    for run_idx = 1:num_stat_runs % 统计运行循环 (num_stat_runs 次)
        % --- 随机种子管理 (关键！) ---
        % 为每一次独立的统计运行设置一个确定且独特的随机种子
        % 这确保了每次运行的随机序列是不同的，但整个实验是可复现的。
        seed_value = (s_idx - 1) * 1000 + run_idx; % 更简洁且不易重复的种子生成方式
        rng(seed_value); % 设置随机种子
        
        fprintf('\n------- 场景 %d/%d, 统计运行 %d/%d -------\n', ...
                s_idx, num_scenarios, run_idx, num_stat_runs);

        % --- 生成每次运行独特的场景内随机元素 (终端位置，任务量，初始雾节点位置等) ---
        % 这些元素在每次 run_idx 循环中都会重新生成，确保每次运行都是新的随机实例
        problem.terminalProperties.positions = problem.area(2,1) + (problem.area(2,2) - problem.area(2,1)) .* rand(problem.nTerminals, 2);
        problem.terminalProperties.Pt_dbm = linspace(10, 15, problem.nTerminals); 
        problem.terminalProperties.fc = linspace(2.4e9, 5.8e9, problem.nTerminals); 
        problem.fogNodeProperties.cpu_cycle_rate = linspace(2e9, 5e9, problem.nFogNodes); 
        
        problem.initial_fog_deployment_flat = problem.area(2,1) + (problem.area(2,2) - problem.area(2,1)) .* rand(1, problem.nFogNodes*2);
        problem.initial_fog_positions_matrix = reshape(problem.initial_fog_deployment_flat, [2, problem.nFogNodes])';
        
        min_task_size = 0.1e6; max_task_size = 1e6;
        problem.terminalProperties.task_sizes = min_task_size + (max_task_size - min_task_size) .* rand(1, problem.nTerminals);
        
        n_half = floor(problem.nTerminals / 2);
        Bmin_vec = [ones(1, n_half) * 0.2e6, ones(1, problem.nTerminals - n_half) * 0.2e6];
        Bmax_vec = [ones(1, n_half) * 10e6, ones(1, problem.nTerminals - n_half) * 10e6];
        problem.bounds.bandwidth = [Bmin_vec; Bmax_vec];
        
        % --- 新增：为本次统计运行生成固定阴影衰落值，并存储到 problem 中 ---
        % EvaluateParticle 将使用这些固定值。确保 randn() 的参数是 nTerminals
        problem.fixed_shadow_LoS_val = shadow_std_dev.LoS * randn(1, problem.nTerminals);
        problem.fixed_shadow_NLoS_val = shadow_std_dev.NLoS * randn(1, problem.nTerminals);
        % -----------------------------------------------------------------
        
        % --- 初始化时隙间的记忆机制 (每轮统计运行开始时重置) ---
        LastSlotArchive_MyNSGA = []; 

%% 7. 最内层循环：时隙循环 (Time Slot Loop)
        for t_slot = 1:nSlots
            fprintf('  开始执行时隙: %d / %d\n', t_slot, nSlots);
            % --- 运行 MyNSGA-II ---
            tic;
            FinalPopulation_MyNSGA_Slot = MyNSGA_II(problem, params_my_nsga, LastSlotArchive_MyNSGA);
            time_my_nsga = toc;
            % 立即存储运行时间
            all_scenario_results.MyNSGA_II.Runtime{s_idx}(run_idx, t_slot) = time_my_nsga;
            Fronts_MyNSGA_Slot = FindAllFronts(FinalPopulation_MyNSGA_Slot);
            FinalArchive_MyNSGA_Slot = getFirstFront(Fronts_MyNSGA_Slot);
            if ~isempty(FinalArchive_MyNSGA_Slot), [FinalArchive_MyNSGA_Slot.RunTime] = deal(time_my_nsga); end
            LastSlotArchive_MyNSGA = FinalArchive_MyNSGA_Slot;

            % --- 运行 原始 NSGA-II ---
            tic;
            FinalPopulation_NSGA_Slot = NSGA_II(problem, params_nsga, []);
            time_nsga = toc;
            % 立即存储运行时间
            all_scenario_results.NSGA_II.Runtime{s_idx}(run_idx, t_slot) = time_nsga;
            Fronts_NSGA_Slot = FindAllFronts(FinalPopulation_NSGA_Slot);
            FinalArchive_NSGA_Slot = getFirstFront(Fronts_NSGA_Slot);
            if ~isempty(FinalArchive_NSGA_Slot), [FinalArchive_NSGA_Slot.RunTime] = deal(time_nsga); end

            % --- 运行 纯血 MOPSO ---
            tic;
            FinalArchive_MOPSO_Slot = MOPSO(problem, params_mopso);
            time_mopso = toc;
            % 立即存储运行时间
            all_scenario_results.MOPSO.Runtime{s_idx}(run_idx, t_slot) = time_mopso;
            if ~isempty(FinalArchive_MOPSO_Slot), [FinalArchive_MOPSO_Slot.RunTime] = deal(time_mopso); end

            % --- 运行 基准方法 (Baseline) ---
            tic;
            % 修正：将预分配大小设置为 N_RUNS_BASELINE
            All_Baseline_Solutions_for_this_slot = repmat(struct('Position',[],'Objectives',[],'Tmax',[]), N_RUNS_BASELINE, 1);
            for b_run_greedy = 1:N_RUNS_BASELINE
                All_Baseline_Solutions_for_this_slot(b_run_greedy) = Baseline_Greedy(problem);
            end
            time_baseline = toc;
            % 立即存储运行时间
            all_scenario_results.Baseline.Runtime{s_idx}(run_idx, t_slot) = time_baseline;
            Fronts_Baseline = FindAllFronts(All_Baseline_Solutions_for_this_slot);
            FinalArchive_Baseline_Slot = getFirstFront(Fronts_Baseline);
            if ~isempty(FinalArchive_Baseline_Slot), [FinalArchive_Baseline_Slot.RunTime] = deal(time_baseline); end

            % --- 收集当前运行的每个算法的最终前沿数据 (用于计算该算法的指标) ---
            % 这部分保持不变，current_scenario_all_run_final_fronts 仍然存储各算法的独立前沿
            current_scenario_all_run_final_fronts{1, t_slot, run_idx} = getObjectivesMatrix(FinalArchive_MyNSGA_Slot);
            current_scenario_all_run_final_fronts{2, t_slot, run_idx} = getObjectivesMatrix(FinalArchive_NSGA_Slot);
            current_scenario_all_run_final_fronts{3, t_slot, run_idx} = getObjectivesMatrix(FinalArchive_MOPSO_Slot);
            current_scenario_all_run_final_fronts{4, t_slot, run_idx} = getObjectivesMatrix(FinalArchive_Baseline_Slot);

            % --- 新增：为构建全局 PF*，收集当前 run_idx 和 t_slot 下所有算法的非支配解 ---
            % 这些解在后续会被合并，然后进行一次全局非支配排序
            
            % 使用 cell 数组收集各个算法的 ObjectivesMatrix，避免动态 vertcat
            temp_objs_collector_for_current_slot = cell(1, length(alg_names_for_results)); 
            
            % 逐一获取各算法的 ObjectivesMatrix 并存入 cell 数组
            if ~isempty(FinalArchive_MyNSGA_Slot)
                temp_objs_collector_for_current_slot{1} = getObjectivesMatrix(FinalArchive_MyNSGA_Slot);
            end
            if ~isempty(FinalArchive_NSGA_Slot)
                temp_objs_collector_for_current_slot{2} = getObjectivesMatrix(FinalArchive_NSGA_Slot);
            end
            if ~isempty(FinalArchive_MOPSO_Slot)
                temp_objs_collector_for_current_slot{3} = getObjectivesMatrix(FinalArchive_MOPSO_Slot);
            end
            if ~isempty(FinalArchive_Baseline_Slot)
                temp_objs_collector_for_current_slot{4} = getObjectivesMatrix(FinalArchive_Baseline_Slot);
            end
            
            % 移除空的 cell 元素，并一次性 vertcat 合并所有有效矩阵
            temp_objs_collector_for_current_slot = temp_objs_collector_for_current_slot(~cellfun('isempty', temp_objs_collector_for_current_slot));
            
            if ~isempty(temp_objs_collector_for_current_slot)
                combined_objs_for_this_slot_run = vertcat(temp_objs_collector_for_current_slot{:});
            else
                combined_objs_for_this_slot_run = []; % 如果所有算法都没找到解，则为空矩阵
            end
            
            % 将这些合并后的解 (不一定是全局非支配的，但已经去除了算法内部的被支配解) 存储起来
            all_combined_solutions_for_global_pf{run_idx, t_slot} = combined_objs_for_this_slot_run;
        end % 时隙循环结束
    end % 统计运行循环结束

%% 8. 后处理：构建 PF* 并计算指标 (在所有 num_stat_runs 次运行结束后，针对当前场景和时隙)
    fprintf('\n===== 场景 %d 运行完毕，正在计算指标... =====\n', s_idx);
    
    for t_slot = 1:nSlots % 遍历每个时隙
        fprintf('  处理时隙 %d 的指标计算...\n', t_slot);
    
        % --- 1. 为当前场景和时隙构建近似参考前沿 (PF*) 和 HV 参考点 ---
    
        % 现在从 all_combined_solutions_for_global_pf 中收集数据
        % 这个集合只包含各算法在各次运行各时隙中找到的非支配解 (相对于自身算法而言)
        % 而不是所有原始种群的解，大大减少了数据量。
        all_objs_for_pf_star = [];
        % 遍历所有统计运行，收集当前时隙的所有算法的非支配解
        for r_idx = 1:num_stat_runs
            current_run_combined_objs = all_combined_solutions_for_global_pf{r_idx, t_slot};
            if ~isempty(current_run_combined_objs)
                all_objs_for_pf_star = [all_objs_for_pf_star; current_run_combined_objs]; %#ok<AGROW> % 合并
            end
        end
        
        % --- 关键修正：过滤惩罚值 (通用化处理) ---
        % 确保 PF* 只包含有效的非惩罚解
        if ~isempty(all_objs_for_pf_star)
            % 假设惩罚值为 1e9 或更高，且没有NaN/Inf。
            valid_objs_mask = ~any(all_objs_for_pf_star >= 1e9 | isnan(all_objs_for_pf_star) | isinf(all_objs_for_pf_star), 2);
            all_objs_for_pf_star = all_objs_for_pf_star(valid_objs_mask, :);
        end
        
        % 构建 IGD 参考前沿 (igd_reference_front_obj)
        igd_reference_front_obj = [];
        if size(all_objs_for_pf_star, 1) > 1
            % 对筛选后的集合进行最终的全局非支配排序
            igd_ref_front_idx = FindNonDominatedSolutions(all_objs_for_pf_star); 
            igd_reference_front_obj = all_objs_for_pf_star(igd_ref_front_idx, :);
            if size(igd_reference_front_obj, 1) < 2
                warning('场景 %d, 时隙 %d 的 IGD 参考前沿点数少于2个，无法计算 IGD。', s_idx, t_slot);
                igd_reference_front_obj = [];
            else
                [~, sort_idx_igd] = sort(igd_reference_front_obj(:,1)); % 排序以便IGD计算
                igd_reference_front_obj = igd_reference_front_obj(sort_idx_igd, :);
            end
        else
            warning('场景 %d, 时隙 %d 的所有运行未能生成足够解来构建 IGD 参考前沿。', s_idx, t_slot);
            igd_reference_front_obj = []; % 如果总解数不足，则无法构建参考前沿
        end
        
        % 确定 HV 参考点
        hv_reference_point = [];
        if ~isempty(all_objs_for_pf_star)
             % 取所有目标的最大值加一个小的余量，作为 HV 的最差参考点（适用于最小化问题）
            hv_reference_point = max(all_objs_for_pf_star, [], 1) * 1.1; 
        end
        
        if isempty(hv_reference_point) || any(isnan(hv_reference_point)) || any(isinf(hv_reference_point))
            warning('场景 %d, 时隙 %d 的 HV 参考点计算出错或为空，使用默认值。', s_idx, t_slot);
            hv_reference_point = [0.8, 0.35]; % 预设一个保守的默认值，请根据问题实际范围调整
        end
        % --- 2. 遍历所有运行，计算指标并存储 ---
        for alg_idx = 1:length(alg_names_for_results)
            current_alg_name = alg_names_for_results{alg_idx};
            for r_idx = 1:num_stat_runs
                % current_run_archive_obj 是一个目标矩阵 (来源于 current_scenario_all_run_final_fronts)
                current_run_archive_obj_matrix = current_scenario_all_run_final_fronts{alg_idx, t_slot, r_idx}; 
                
                % 为了兼容 CalculateMetricsOnly，需要把目标矩阵和运行时间包在一个结构体里
                % RunTime 值之前已存储在 all_scenario_results.(current_alg_name).Runtime{s_idx}(r_idx, t_slot) 中
                runtime_for_this_run = all_scenario_results.(current_alg_name).Runtime{s_idx}(r_idx, t_slot);
                temp_archive_for_metrics = createArchiveFromObjectives(current_run_archive_obj_matrix, runtime_for_this_run);
                % 调用 CalculateMetricsOnly 来获取指标值
                metrics_current_run = CalculateMetricsOnly(temp_archive_for_metrics, ...
                    igd_reference_front_obj, hv_reference_point);
                
                % 存储指标到最终结果结构体中
                all_scenario_results.(current_alg_name).IGD{s_idx}(r_idx, t_slot) = metrics_current_run.IGD;
                all_scenario_results.(current_alg_name).HV{s_idx}(r_idx, t_slot) = metrics_current_run.HV;
                all_scenario_results.(current_alg_name).Spacing{s_idx}(r_idx, t_slot) = metrics_current_run.Spacing;
                all_scenario_results.(current_alg_name).Spread{s_idx}(r_idx, t_slot) = metrics_current_run.Spread;
                all_scenario_results.(current_alg_name).NumFeasibleSolutions{s_idx}(r_idx, t_slot) = metrics_current_run.NumFeasibleSolutions;
                % Runtime 已经在之前的循环中直接存储了
            end
        end
    end % 时隙循环 (指标计算) 结束
    fprintf('\n===== 场景 %d 指标计算完成 =====\n', s_idx);
    
    % 清除本场景的临时全局PF*数据，释放内存
    clear all_combined_solutions_for_global_pf; 
end % 场景循环结束
%% 9. 最终结果展示和保存 (所有场景运行和指标计算完毕后)
fprintf('\n============== 所有场景执行完毕，正在输出最终结果 ===============\n');
% --- 新增：准备保存的结构体 ---
% 这个结构体将存储所有场景、所有算法的平均值和标准差
statistical_results_to_save = struct(); 
% --- 选择要展示和绘图的场景和时隙 ---
s_idx_to_display_plot = num_scenarios; % 只选择最后一个场景
t_slot_to_display_plot = nSlots;       % 只选择最后一个时隙
% 获取指定场景的配置信息，用于后续保存和识别
current_scenario_display_info = experimental_scenarios{s_idx_to_display_plot};
current_scenario_display_name = sprintf('S%d_I%d_M%d_R2%dx%d', ...
                                        s_idx_to_display_plot, current_scenario_display_info{1}, ...
                                        current_scenario_display_info{2}, ...
                                        current_scenario_display_info{3}, ...
                                        current_scenario_display_info{4});
% 在 statistical_results_to_save 中为指定场景创建一个字段
statistical_results_to_save.(current_scenario_display_name) = struct();
fprintf('\n--- 场景 %d (I=%d, M=%d, R2=%dx%d) 的平均性能 (时隙 %d) ---\n', ...
        s_idx_to_display_plot, current_scenario_display_info{1}, ...
        current_scenario_display_info{2}, current_scenario_display_info{3}, ...
        current_scenario_display_info{4}, t_slot_to_display_plot);
fprintf('------------------------------------------------------------------\n');
% 在当前场景结构体中为指定时隙创建一个字段
slot_field_name_for_save = sprintf('Slot_%d', t_slot_to_display_plot);
statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save) = struct();
for alg_name = alg_names_for_results
    current_alg_name = alg_name{1};
    
    % 计算均值和标准差，并忽略 NaN 值
    avg_igd = mean(all_scenario_results.(current_alg_name).IGD{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    std_igd = std(all_scenario_results.(current_alg_name).IGD{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    
    avg_hv = mean(all_scenario_results.(current_alg_name).HV{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    std_hv = std(all_scenario_results.(current_alg_name).HV{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    
    avg_spacing = mean(all_scenario_results.(current_alg_name).Spacing{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    std_spacing = std(all_scenario_results.(current_alg_name).Spacing{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    
    avg_spread = mean(all_scenario_results.(current_alg_name).Spread{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    std_spread = std(all_scenario_results.(current_alg_name).Spread{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    avg_runtime = mean(all_scenario_results.(current_alg_name).Runtime{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    std_runtime = std(all_scenario_results.(current_alg_name).Runtime{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    avg_num_feasible = mean(all_scenario_results.(current_alg_name).NumFeasibleSolutions{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    
    fprintf('    %s:\n', current_alg_name);
    fprintf('      IGD: %.4f (%.4f)\n', avg_igd, std_igd);
    fprintf('      HV: %.4f (%.4f)\n', avg_hv, std_hv);
    fprintf('      Spacing: %.4f (%.4f)\n', avg_spacing, std_spacing);
    fprintf('      Spread: %.4f (%.4f)\n', avg_spread, std_spread);
    fprintf('      Runtime: %.4f (%.4f) s\n', avg_runtime, std_runtime);
    fprintf('      Avg Feasible Solutions: %.1f\n', avg_num_feasible);
    % --- 将计算出的平均值和标准差保存到结构体中 ---
    % 使用算法名称作为字段名
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).IGD_avg = avg_igd;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).IGD_std = std_igd;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).HV_avg = avg_hv;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).HV_std = std_hv;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Spacing_avg = avg_spacing;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Spacing_std = std_spacing;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Spread_avg = avg_spread;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Spread_std = std_spread;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Runtime_avg = avg_runtime;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Runtime_std = std_runtime;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).NumFeasibleSolutions_avg = avg_num_feasible;
end
fprintf('  ----------------------------------------\n');
% --- 保存统计结果到 .mat 文件 ---
output_filename = fullfile(output_timestamp_folder, sprintf('Experiment_Statistical_Results_%s_Slot%d_%s.mat', current_scenario_display_name, t_slot_to_display_plot, string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
save(output_filename, 'statistical_results_to_save', 'alg_names_for_results', 'metric_names', 'experimental_scenarios', 'nSlots', 'num_stat_runs');
fprintf('\n所有统计结果已保存到文件: %s\n', output_filename);
fprintf('\n============== 所有场景执行完毕，正在输出最终结果 ===============\n');
% --- 辅助函数：提取第一非支配前沿 (通常 FindAllFronts 已经返回 cell 数组) ---
function first_front_archive = getFirstFront(fronts_cell_array)
    if ~isempty(fronts_cell_array) && ~isempty(fronts_cell_array{1})
        first_front_archive = fronts_cell_array{1};
    else
        first_front_archive = [];
    end
end

% --- 辅助函数：从存档结构体中提取目标值矩阵 ---
function obj_matrix = getObjectivesMatrix(archive_struct)
    if ~isempty(archive_struct)
        obj_matrix = vertcat(archive_struct.Objectives);
        % 确保只包含有效解，这里的过滤方式应与计算指标时保持一致
        obj_matrix = obj_matrix(obj_matrix(:,1) < 1e9 & obj_matrix(:,2) < 1e9, :); 
    else
        obj_matrix = [];
    end
end

% --- 辅助函数：创建Archive结构体以兼容CalculateMetricsOnly的输入 ---
% 这个函数用于将纯粹的目标矩阵转换为 CalculateMetricsOnly 所需的 archive 结构体格式
function archive_struct_out = createArchiveFromObjectives(obj_matrix, runtime_val)
    if isempty(obj_matrix)
        archive_struct_out = [];
        return;
    end
    % 根据 obj_matrix 的行数创建结构体数组
    archive_struct_out = repmat(struct('Objectives', [], 'RunTime', NaN), size(obj_matrix, 1), 1);
    
    for i = 1:size(obj_matrix, 1)
        archive_struct_out(i).Objectives = obj_matrix(i,:);
        archive_struct_out(i).RunTime = runtime_val; % 所有解共享本次运行的运行时间
    end
end

% --- 辅助函数：合并两个结构体 (MATLAB R2016b+ 支持) ---
function s1 = mergeStructs(s1, s2)
    names = fieldnames(s2);
    for i = 1:length(names)
        s1.(names{i}) = s2.(names{i});
    end
end