%% LMM Model Selection Helper
% This code helps you choose the best Linear Mixed-Effects Model
% for your pressure experiment with nested Subject/IC structure

function best_model = select_best_lmm_model(data)
    % INPUT: data table with columns: Response, Pressure, Subject, IC
    % OUTPUT: best_model structure with model and selection info
    
    fprintf('=== LMM Model Selection Process ===\n\n');
    
    %% Step 1: Data Preparation
    fprintf('Step 1: Preparing data...\n');
    
    % Ensure categorical variables
    data.Pressure = categorical(data.Pressure);
    data.Subject = categorical(data.Subject);
    data.IC = categorical(data.IC);
    
    % Create nested IC identifier (unique across subjects)
    data.IC_nested = categorical(strcat(string(data.Subject), '_', string(data.IC)));
    
    % Display data summary
    fprintf('Data summary:\n');
    fprintf('  - Subjects: %d\n', length(unique(data.Subject)));
    fprintf('  - ICs per subject: %.1f (average)\n', length(unique(data.IC_nested))/length(unique(data.Subject)));
    fprintf('  - Pressure levels: %d\n', length(unique(data.Pressure)));
    fprintf('  - Total observations: %d\n\n', height(data));
    
    %% Step 2: Define Models to Test
    fprintf('Step 2: Testing different model structures...\n');
    
    models = {};
    model_names = {};
    model_descriptions = {};
    
    % Model 1: No random effects (baseline)
    models{1} = 'Response ~ Pressure';
    model_names{1} = 'No Random Effects';
    model_descriptions{1} = 'Basic ANOVA-like model (for comparison)';
    
    % Model 2: Random intercepts only - Subject
    models{2} = 'Response ~ Pressure + (1|Subject)';
    model_names{2} = 'Subject Intercept';
    model_descriptions{2} = 'Random intercepts for subjects only';
    
    % Model 3: Random intercepts only - IC nested
    models{3} = 'Response ~ Pressure + (1|IC_nested)';
    model_names{3} = 'IC Intercept';
    model_descriptions{3} = 'Random intercepts for ICs only';
    
    % Model 4: Random intercepts - Both Subject and IC
    models{4} = 'Response ~ Pressure + (1|Subject) + (1|IC_nested)';
    model_names{4} = 'Subject + IC Intercepts';
    model_descriptions{4} = 'Random intercepts for both subjects and ICs';
    
    % Model 5: Random slopes - Subject only
    models{5} = 'Response ~ Pressure + (Pressure|Subject)';
    model_names{5} = 'Subject Intercept + Slopes';
    model_descriptions{5} = 'Random intercepts and slopes for subjects';
    
    % Model 6: Random slopes - IC only  
    models{6} = 'Response ~ Pressure + (Pressure|IC_nested)';
    model_names{6} = 'IC Intercept + Slopes';
    model_descriptions{6} = 'Random intercepts and slopes for ICs';
    
    % Model 7: Random intercepts + Subject slopes
    models{7} = 'Response ~ Pressure + (1|Subject) + (1|IC_nested) + (Pressure-1|Subject)';
    model_names{7} = 'Intercepts + Subject Slopes';
    model_descriptions{7} = 'Random intercepts (both) + subject slopes only';
    
    % Model 8: Random intercepts + IC slopes
    models{8} = 'Response ~ Pressure + (1|Subject) + (1|IC_nested) + (Pressure-1|IC_nested)';
    model_names{8} = 'Intercepts + IC Slopes';
    model_descriptions{8} = 'Random intercepts (both) + IC slopes only';
    
    % Model 9: Full model - both intercepts and slopes
    models{9} = 'Response ~ Pressure + (Pressure|Subject) + (Pressure|IC_nested)';
    model_names{9} = 'Full Random Effects';
    model_descriptions{9} = 'Random intercepts and slopes for both subjects and ICs';
    
    %% Step 3: Fit Models and Collect Results
    fprintf('Fitting models...\n');
    
    fitted_models = cell(length(models), 1);
    model_stats = table();
    
    for i = 1:length(models)
        fprintf('  Fitting Model %d: %s...', i, model_names{i});
        
        try
            if i == 1
                % Use fitlm for the no random effects model
                fitted_models{i} = fitlm(data, models{i});
                aic_val = fitted_models{i}.ModelCriterion.AIC;
                bic_val = fitted_models{i}.ModelCriterion.BIC;
                loglik_val = fitted_models{i}.LogLikelihood;
                converged = true;
                n_params = fitted_models{i}.NumCoefficients + 1; % +1 for residual variance
            else
                % Use fitlme for mixed effects models
                fitted_models{i} = fitlme(data, models{i}, 'FitMethod', 'REML');
                aic_val = fitted_models{i}.ModelCriterion.AIC;
                bic_val = fitted_models{i}.ModelCriterion.BIC;
                loglik_val = fitted_models{i}.LogLikelihood;
                
                % Check convergence - different ways depending on MATLAB version
                try
                    if isprop(fitted_models{i}, 'Optimizer')
                        converged = fitted_models{i}.Optimizer.Converged;
                    elseif isfield(fitted_models{i}, 'Optimizer')
                        converged = fitted_models{i}.Optimizer.Converged;
                    else
                        % Assume converged if model was successfully created
                        converged = true;
                    end
                catch
                    % If we can't check convergence, assume it converged
                    converged = true;
                end
                
                n_params = fitted_models{i}.NumParameters;
            end
            
            % Store results
            model_stats = [model_stats; table(i, model_names(i), aic_val, bic_val, loglik_val, converged, n_params, ...
                'VariableNames', {'Model_Num', 'Model_Name', 'AIC', 'BIC', 'LogLikelihood', 'Converged', 'NumParams'})];
            
            fprintf(' ✓\n');
            
        catch ME
            fprintf(' ✗ (Failed: %s)\n', ME.message);
            fitted_models{i} = [];
            model_stats = [model_stats; table(i, model_names(i), NaN, NaN, NaN, false, NaN, ...
                'VariableNames', {'Model_Num', 'Model_Name', 'AIC', 'BIC', 'LogLikelihood', 'Converged', 'NumParams'})];
        end
    end
    
    %% Step 4: Model Comparison and Selection
    fprintf('\n=== MODEL COMPARISON RESULTS ===\n\n');
    
    % Filter only converged models
    converged_models = model_stats(model_stats.Converged == true, :);
    
    if isempty(converged_models)
        error('No models converged successfully!');
    end
    
    % Sort by AIC (lower is better)
    converged_models = sortrows(converged_models, 'AIC');
    
    fprintf('Converged Models (ranked by AIC):\n');
    fprintf('%-4s %-25s %-10s %-10s %-12s %-8s\n', 'Rank', 'Model', 'AIC', 'BIC', 'LogLik', 'Params');
    fprintf('%-4s %-25s %-10s %-10s %-12s %-8s\n', '----', '-------------------------', '------', '------', '--------', '------');
    
    for i = 1:height(converged_models)
        fprintf('%-4d %-25s %-10.2f %-10.2f %-12.2f %-8d\n', ...
            i, converged_models.Model_Name{i}, converged_models.AIC(i), ...
            converged_models.BIC(i), converged_models.LogLikelihood(i), converged_models.NumParams(i));
    end
    
    %% Step 5: Statistical Tests for Model Selection
    fprintf('\n=== MODEL COMPARISON TESTS ===\n\n');
    
    % Compare nested models using likelihood ratio tests
    best_model_idx = converged_models.Model_Num(1);
    best_model = fitted_models{best_model_idx};
    
    fprintf('Best model by AIC: %s\n', converged_models.Model_Name{1});
    fprintf('Best model by BIC: %s\n\n', converged_models.Model_Name{converged_models.BIC == min(converged_models.BIC)});
    
    % Test if random effects are needed
    if any(converged_models.Model_Num == 4) && any(converged_models.Model_Num == 1)
        fprintf('Testing if random effects are necessary:\n');
        model_basic = fitted_models{1};
        model_random = fitted_models{4};
        
        % Note: Direct comparison between lm and lme is not straightforward
        % This is a conceptual comparison
        fprintf('  AIC difference (Random vs Basic): %.2f\n', ...
            converged_models.AIC(converged_models.Model_Num == 4) - ...
            converged_models.AIC(converged_models.Model_Num == 1));
    end
    
    %% Step 6: Final Recommendations
    fprintf('\n=== RECOMMENDATIONS ===\n\n');
    
    fprintf('1. BEST MODEL: %s\n', converged_models.Model_Name{1});
    fprintf('   Formula: %s\n', models{best_model_idx});
    fprintf('   Reason: Lowest AIC among converged models\n\n');
    
    % Check for overfitting warning
    if converged_models.NumParams(1) > height(data)/10
        fprintf('⚠️  WARNING: Best model has many parameters relative to data size.\n');
        fprintf('   Consider simpler model: %s\n\n', converged_models.Model_Name{2});
    end
    
    % BIC recommendation (more conservative)
    [~, bic_best_idx] = min(converged_models.BIC);
    if bic_best_idx ~= 1
        fprintf('2. CONSERVATIVE CHOICE: %s\n', converged_models.Model_Name{bic_best_idx});
        fprintf('   Reason: Lowest BIC (penalizes complexity more)\n\n');
    end
    
    %% Step 7: Diagnostic Plots for Best Model
    fprintf('3. DIAGNOSTIC PLOTS for best model:\n');
    
    if best_model_idx > 1  % Mixed effects model
        figure('Name', 'Model Diagnostics');
        subplot(2,2,1);
        plotResiduals(best_model, 'fitted');
        title('Residuals vs Fitted');
        
        subplot(2,2,2);
        plotResiduals(best_model, 'probability');
        title('Normal Probability Plot');
        
        subplot(2,2,3);
        plotResiduals(best_model, 'histogram');
        title('Residual Histogram');
        
        subplot(2,2,4);
        plotResiduals(best_model, 'lagged');
        title('Residual Autocorrelation');
        
        fprintf('   Check diagnostic plots for model assumptions.\n\n');
    end
    
    %% Step 8: Test Your Research Question
    fprintf('4. TESTING PRESSURE EFFECT:\n');
    
    try
        if best_model_idx == 1
            % Linear model
            anova_result = anova(best_model);
            % For fitlm, find pressure row
            pressure_idx = contains(string(anova_result.Properties.RowNames), 'Pressure', 'IgnoreCase', true);
            if any(pressure_idx)
                p_value = anova_result.pValue(pressure_idx);
            else
                p_value = NaN;
            end
        else
            % Mixed effects model
            anova_result = anova(best_model, 'DFMethod', 'satterthwaite');
            % For fitlme, check if Term column exists, otherwise use row names
            if ismember('Term', anova_result.Properties.VariableNames)
                pressure_idx = strcmp(anova_result.Term, 'Pressure');
                p_value = anova_result.pValue(pressure_idx);
            else
                % Try using row names
                pressure_idx = contains(string(anova_result.Properties.RowNames), 'Pressure', 'IgnoreCase', true);
                if any(pressure_idx)
                    p_value = anova_result.pValue(pressure_idx);
                else
                    p_value = NaN;
                end
            end
        end
        
        if ~isnan(p_value) && ~isempty(p_value)
            fprintf('   Pressure effect p-value: %.6f\n', p_value);
            if p_value < 0.05
                fprintf('   ✓ Significant pressure effect found!\n');
            else
                fprintf('   ✗ No significant pressure effect detected.\n');
            end
        else
            fprintf('   Could not extract pressure p-value automatically.\n');
            fprintf('   Run: anova(best_model) to see full results.\n');
            p_value = NaN;
        end
        
    catch ME
        fprintf('   Error testing pressure effect: %s\n', ME.message);
        fprintf('   Run: anova(best_model) manually to see results.\n');
        p_value = NaN;
    end
    
    %% Return Results
    best_model_info.model = best_model;
    best_model_info.formula = models{best_model_idx};
    best_model_info.name = converged_models.Model_Name{1};
    best_model_info.stats = converged_models;
    best_model_info.pressure_pvalue = p_value;
    
    if nargout > 0
        best_model = best_model_info;
    end
    
    fprintf('\nModel selection complete!\n');
end

%% Example Usage:
% Assuming your data is in a table called 'mydata' with columns:
% Response, Pressure, Subject, IC
%
% best_model = select_best_lmm_model(mydata);
%
% % Access the best model:
% final_model = best_model.model;
% 
% % Get detailed results:
% anova(final_model, 'DFMethod', 'satterthwaite')
% fixedEffects(final_model)