classdef CRMTool < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        CRMToolUIFigure                 matlab.ui.Figure
        CollectiveRiskModelLabel        matlab.ui.control.Label
        SeverityDistributionPanel       matlab.ui.container.Panel
        SeverityLabel                   matlab.ui.control.Label
        SeverityDropDown                matlab.ui.control.DropDown
        FrequencyDistributionPanel      matlab.ui.container.Panel
        FrequencyLabel                  matlab.ui.control.Label
        FrequencyDropDown               matlab.ui.control.DropDown
        NumericalInversionAlgorithmPanel  matlab.ui.container.Panel
        AlgorithmLabel                  matlab.ui.control.Label
        AlgorithmDropDown               matlab.ui.control.DropDown
        OptionsPanel                    matlab.ui.container.Panel
        NLabel                          matlab.ui.control.Label
        NEditField                      matlab.ui.control.EditField
        SixSigmaRuleLabel               matlab.ui.control.Label
        SixSigmaRuleEditField           matlab.ui.control.EditField
        SeverityDistributionParametersPanel  matlab.ui.container.Panel
        Parameter2EditFieldLabel        matlab.ui.control.Label
        Parameter2EditField             matlab.ui.control.EditField
        Parameter3EditFieldLabel        matlab.ui.control.Label
        Parameter3EditField             matlab.ui.control.EditField
        Label                           matlab.ui.control.Label
        Label_2                         matlab.ui.control.Label
        Label_3                         matlab.ui.control.Label
        Parameter1EditFieldLabel        matlab.ui.control.Label
        Parameter1EditField             matlab.ui.control.EditField
        FrequencyDistributionParametersPanel  matlab.ui.container.Panel
        Parameter1Label                 matlab.ui.control.Label
        Parameter1EditField_2           matlab.ui.control.EditField
        Parameter2Label                 matlab.ui.control.Label
        Parameter2EditField_2           matlab.ui.control.EditField
        Parameter3Label                 matlab.ui.control.Label
        Parameter3EditField_2           matlab.ui.control.EditField
        Label_4                         matlab.ui.control.Label
        Label_5                         matlab.ui.control.Label
        Label_6                         matlab.ui.control.Label
        InputArgumentsPanel             matlab.ui.container.Panel
        VaRprobabilitiesEditFieldLabel  matlab.ui.control.Label
        VaRprobabilitiesEditField       matlab.ui.control.EditField
        XvaluesLabel                    matlab.ui.control.Label
        XvaluesEditField                matlab.ui.control.EditField
        ResultsPanel                    matlab.ui.container.Panel
        UIAxes                          matlab.ui.control.UIAxes
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxes3                         matlab.ui.control.UIAxes
        CALCULATEButton                 matlab.ui.control.Button
        EXITButton                      matlab.ui.control.Button
        SaveSwitchLabel                 matlab.ui.control.Label
        SaveSwitch                      matlab.ui.control.Switch
        VaRquantilesLabel               matlab.ui.control.Label
        VaRquantilesEditField           matlab.ui.control.EditField
        TableSwitchLabel                matlab.ui.control.Label
        TableSwitch                     matlab.ui.control.Switch
        AggregateLossDistributionCalculatorLabel  matlab.ui.control.Label
        AboutButton                     matlab.ui.control.Button
    end

    
    properties (Access = private)
        cfXParameter1     % Parameter 1 of the Severity Distribution
        cfXParameter2     % Parameter 2 of the Severity Distribution
        cfXParameter3     % Parameter 3 of the Severity Distribution
        cfNParameter1     % Parameter 1 of the Frequency Distribution
        cfNParameter2     % Parameter 2 of the Frequency Distribution
        cfNParameter3     % Parameter 3 of the Frequency Distribution
        options           % Options parameters for the Algorithm
        x                 % x values
        prob              % prob values for computing the VaRs
        cfX               % CF of the severity distribution
        cf                % CF of the compound (frequency & severity) distribution
        Result            % Result of the Inversion Algorithm
        Save              % Save result status
        Table             % Table generate status
        SeverityID        % Description
        FrequencyID       % Description
    end
    

    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
            % Ver.: 30-Jun-2018 11:40:35
            warning('off','MATLAB:subscripting:noSubscriptsSpecified')
            warning('off','MATLAB:Figure:FigureSavedToMATFile')
            warning('off','MATLAB:ui:uifigure:UnsupportedAppDesignerFunctionality')
            warning('off','MATLAB:lang:cannotClearExecutingFunction')
            warning('off','VW:CharFunTool:cf2DistGP')
        end

        % Value changed function: AlgorithmDropDown, 
        % FrequencyDropDown, SeverityDropDown
        function ValueChanged(app, event)
            % (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
            % Ver.: 30-Jun-2018 11:40:35
            value = app.SeverityDropDown.Value;
            if ~strcmp(value,app.SeverityID)
                app.Parameter1EditField.Value = '';
                app.Parameter2EditField.Value = '';
                app.Parameter3EditField.Value = '';
            end
            app.SeverityID = value;
            app.Parameter1EditField.Enable = 'off';
            app.Parameter2EditField.Enable = 'off';
            app.Parameter3EditField.Enable = 'off';
            app.Parameter1EditFieldLabel.Enable = 'off';
            app.Parameter2EditFieldLabel.Enable = 'off';
            app.Parameter3EditFieldLabel.Enable = 'off';
            app.Parameter1EditFieldLabel.Text = '';
            app.Parameter2EditFieldLabel.Text = '';
            app.Parameter3EditFieldLabel.Text = '';
            switch value
                case app.SeverityDropDown.Items{1}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = '5';
                    end
                    app.Parameter2EditField.Value = '';
                    app.Parameter3EditField.Value = '';
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'off';
                    app.Parameter3EditField.Enable = 'off';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'off';
                    app.Parameter3EditFieldLabel.Enable = 'off';
                    app.Parameter1EditFieldLabel.Text = 'lambda';
                    app.Parameter2EditFieldLabel.Text = '';
                    app.Parameter3EditFieldLabel.Text = '';
                    app.Label.Text = 'rate parameter, lambda > 0';
                    app.Label_2.Text = '';
                    app.Label_3.Text = '';
                    app.cfXParameter1 = str2double(app.Parameter1EditField.Value);
                    % lambda = app.cfXParameter1
                    app.cfX = @(t) cf_Exponential(t,app.cfXParameter1);
                case app.SeverityDropDown.Items{2}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = '3';
                    end
                    if isempty(app.Parameter2EditField.Value)
                        app.Parameter2EditField.Value = '5';
                    end
                    app.Parameter3EditField.Value = '';
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'on';
                    app.Parameter3EditField.Enable = 'off';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'on';
                    app.Parameter3EditFieldLabel.Enable = 'off';
                    app.Parameter1EditFieldLabel.Text = 'DF1';
                    app.Parameter2EditFieldLabel.Text = 'DF2';
                    app.Parameter3EditFieldLabel.Text = '';
                    app.Label.Text = 'degrees of freedom, DF1 > 0';
                    app.Label_2.Text = 'degrees of freedom, DF2 > 0';
                    app.Label_3.Text = '';
                    app.cfXParameter1 = str2double(app.Parameter1EditField.Value);
                    app.cfXParameter2 = str2double(app.Parameter2EditField.Value);
                    % DF1 = app.cfXParameter1
                    % DF2 = app.cfXParameter2
                    app.cfX = @(t) cf_FisherSnedecor(t,app.cfXParameter1,app.cfXParameter2);
                case app.SeverityDropDown.Items{3}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = '2';
                    end
                    if isempty(app.Parameter2EditField.Value)
                        app.Parameter2EditField.Value = '5';
                    end
                    app.Parameter3EditField.Value = '';
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'on';
                    app.Parameter3EditField.Enable = 'off';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'on';
                    app.Parameter3EditFieldLabel.Enable = 'off';
                    app.Parameter1EditFieldLabel.Text = 'alpha';
                    app.Parameter2EditFieldLabel.Text = 'beta';
                    app.Parameter3EditFieldLabel.Text = '';
                    app.Label.Text = 'shape, alpha > 0';
                    app.Label_2.Text = 'rate, beta > 0';
                    app.Label_3.Text = '';
                    app.cfXParameter1 = str2double(app.Parameter1EditField.Value);
                    app.cfXParameter2 = str2double(app.Parameter2EditField.Value);
                    % alpha = app.cfXParameter1
                    % beta = app.cfXParameter2
                    app.cfX = @(t) cf_Gamma(t,app.cfXParameter1,app.cfXParameter2);
                case app.SeverityDropDown.Items{4}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = '1';
                    end
                    if isempty(app.Parameter2EditField.Value)
                        app.Parameter2EditField.Value = '1';
                    end
                    if isempty(app.Parameter3EditField.Value)
                        app.Parameter3EditField.Value = '0';
                    end
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'on';
                    app.Parameter3EditField.Enable = 'on';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'on';
                    app.Parameter3EditFieldLabel.Enable = 'on';
                    app.Parameter1EditFieldLabel.Text = 'xi';
                    app.Parameter2EditFieldLabel.Text = 'sigma';
                    app.Parameter3EditFieldLabel.Text = 'theta';
                    app.Label.Text = 'shape, here xi >= 0';
                    app.Label_2.Text = 'scale, sigma > 0';
                    app.Label_3.Text = 'threshold, theta >= 0';
                    app.cfXParameter1 = str2double(app.Parameter1EditField.Value);
                    app.cfXParameter2 = str2double(app.Parameter2EditField.Value);
                    app.cfXParameter3 = str2double(app.Parameter3EditField.Value);
                    % xi = app.cfXParameter1
                    % sigma = app.cfXParameter2
                    % theta = app.cfXParameter3
                    app.cfX = @(t) cfX_GeneralizedPareto(t,app.cfXParameter1,app.cfXParameter2,app.cfXParameter3);
                case app.SeverityDropDown.Items{5}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = '3';
                    end
                    app.Parameter2EditField.Value = '';
                    app.Parameter3EditField.Value = '';
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'off';
                    app.Parameter3EditField.Enable = 'off';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'off';
                    app.Parameter3EditFieldLabel.Enable = 'off';
                    app.Parameter1EditFieldLabel.Text = 'DF';
                    app.Parameter2EditFieldLabel.Text = '';
                    app.Parameter3EditFieldLabel.Text = '';
                    app.Label.Text = 'degrees of freedom, DF > 0';
                    app.Label_2.Text = '';
                    app.Label_3.Text = '';
                    app.cfXParameter1 = str2double(app.Parameter1EditField.Value);
                    % DF = app.cfXParameter1
                    app.cfX = @(t) cf_ChiSquare(t,app.cfXParameter1);
                case app.SeverityDropDown.Items{6}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = '2';
                    end
                    if isempty(app.Parameter2EditField.Value)
                        app.Parameter2EditField.Value = '5';
                    end
                    app.Parameter3EditField.Value = '';
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'on';
                    app.Parameter3EditField.Enable = 'off';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'on';
                    app.Parameter3EditFieldLabel.Enable = 'off';
                    app.Parameter1EditFieldLabel.Text = 'alpha';
                    app.Parameter2EditFieldLabel.Text = 'beta';
                    app.Parameter3EditFieldLabel.Text = '';
                    app.Label.Text = 'shape, alpha > 0';
                    app.Label_2.Text = 'rate, beta > 0';
                    app.Label_3.Text = '';
                    app.cfXParameter1 = str2double(app.Parameter1EditField.Value);
                    app.cfXParameter2 = str2double(app.Parameter2EditField.Value);
                    % alpha = app.cfXParameter1
                    % beta = app.cfXParameter2
                    app.cfX = @(t) cf_InverseGamma(t,app.cfXParameter1,app.cfXParameter2);
                case app.SeverityDropDown.Items{7}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = '1';
                    end
                    if isempty(app.Parameter2EditField.Value)
                        app.Parameter2EditField.Value = '1';
                    end
                    app.Parameter3EditField.Value = '';
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'on';
                    app.Parameter3EditField.Enable = 'off';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'on';
                    app.Parameter3EditFieldLabel.Enable = 'off';
                    app.Parameter1EditFieldLabel.Text = 'alpha';
                    app.Parameter2EditFieldLabel.Text = 'beta';
                    app.Parameter3EditFieldLabel.Text = '';
                    app.Label.Text = 'scale, alpha > 0';
                    app.Label_2.Text = 'shape, beta > 0';
                    app.Label_3.Text = '';
                    app.cfXParameter1 = str2double(app.Parameter1EditField.Value);
                    app.cfXParameter2 = str2double(app.Parameter2EditField.Value);
                    % alpha = app.cfXParameter1
                    % beta = app.cfXParameter2
                    app.cfX = @(t) cfX_LogLogistic(t,app.cfXParameter1,app.cfXParameter2);
                case app.SeverityDropDown.Items{8}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = '0';
                    end
                    if isempty(app.Parameter2EditField.Value)
                        app.Parameter2EditField.Value = '1';
                    end
                    app.Parameter3EditField.Value = '';
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'on';
                    app.Parameter3EditField.Enable = 'off';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'on';
                    app.Parameter3EditFieldLabel.Enable = 'off';
                    app.Parameter1EditFieldLabel.Text = 'mu';
                    app.Parameter2EditFieldLabel.Text = 'sigma';
                    app.Parameter3EditFieldLabel.Text = '';
                    app.Label.Text = 'parameters mu (real)';
                    app.Label_2.Text = 'parameter sigma > 0';
                    app.Label_3.Text = '';
                    app.cfXParameter1 = str2double(app.Parameter1EditField.Value);
                    app.cfXParameter2 = str2double(app.Parameter2EditField.Value);
                    % mu = app.cfXParameter1
                    % sigma = app.cfXParameter2
                    app.cfX = @(t) cfX_LogNormal(t,app.cfXParameter1,app.cfXParameter2);
                case app.SeverityDropDown.Items{9}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = '1';
                    end
                    if isempty(app.Parameter2EditField.Value)
                        app.Parameter2EditField.Value = '1';
                    end
                    if isempty(app.Parameter3EditField.Value)
                        app.Parameter3EditField.Value = 'TypeI';
                    end
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'on';
                    app.Parameter3EditField.Enable = 'on';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'on';
                    app.Parameter3EditFieldLabel.Enable = 'on';
                    app.Parameter1EditFieldLabel.Text = 'alpha';
                    app.Parameter2EditFieldLabel.Text = 'sigma';
                    app.Parameter3EditFieldLabel.Text = 'type';
                    app.Label.Text = 'shape, alpha > 0';
                    app.Label_2.Text = 'scale, sigma > 0';
                    app.Label_3.Text = 'type parameter (TypeI or TypeII)';
                    app.cfXParameter1 = str2double(app.Parameter1EditField.Value);
                    app.cfXParameter2 = str2double(app.Parameter2EditField.Value);
                    app.cfXParameter3 = app.Parameter3EditField.Value;
                    % alpha = app.cfXParameter1
                    % sigma = app.cfXParameter2
                    % type = app.cfXParameter3
                    app.cfX = @(t) cfX_Pareto(t,app.cfXParameter1,app.cfXParameter2,app.cfXParameter3);
                case app.SeverityDropDown.Items{10}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = '1';
                    end
                    if isempty(app.Parameter2EditField.Value)
                        app.Parameter2EditField.Value = '1';
                    end
                    app.Parameter3EditField.Value = '';
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'on';
                    app.Parameter3EditField.Enable = 'off';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'on';
                    app.Parameter3EditFieldLabel.Enable = 'off';
                    app.Parameter1EditFieldLabel.Text = 'alpha';
                    app.Parameter2EditFieldLabel.Text = 'beta';
                    app.Parameter3EditFieldLabel.Text = '';
                    app.Label.Text = 'shape, alpha > 0';
                    app.Label_2.Text = 'scale, beta > 0';
                    app.Label_3.Text = '';
                    app.cfXParameter1 = str2double(app.Parameter1EditField.Value);
                    app.cfXParameter2 = str2double(app.Parameter2EditField.Value);
                    % alpha = app.cfXParameter1
                    % beta = app.cfXParameter2
                    app.cfX = @(t) cfX_PearsonV(t,app.cfXParameter1,app.cfXParameter2);
                case app.SeverityDropDown.Items{11}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = '1';
                    end
                    if isempty(app.Parameter2EditField.Value)
                        app.Parameter2EditField.Value = '1';
                    end
                    app.Parameter3EditField.Value = '';
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'on';
                    app.Parameter3EditField.Enable = 'off';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'on';
                    app.Parameter3EditFieldLabel.Enable = 'off';
                    app.Parameter1EditFieldLabel.Text = 'alpha';
                    app.Parameter2EditFieldLabel.Text = 'beta';
                    app.Parameter3EditFieldLabel.Text = '';
                    app.Label.Text = 'shape, alpha > 0';
                    app.Label_2.Text = 'scale, beta > 0';
                    app.Label_3.Text = '';
                    app.cfXParameter1 = str2double(app.Parameter1EditField.Value);
                    app.cfXParameter2 = str2double(app.Parameter2EditField.Value);
                    % alpha = app.cfXParameter1
                    % beta = app.cfXParameter2
                    app.cfX = @(t) cfX_PearsonVI(t,app.cfXParameter1,app.cfXParameter2);
                case app.SeverityDropDown.Items{12}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = '1';
                    end
                    if isempty(app.Parameter2EditField.Value)
                        app.Parameter2EditField.Value = '1';
                    end
                    app.Parameter3EditField.Value = '';
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'on';
                    app.Parameter3EditField.Enable = 'off';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'on';
                    app.Parameter3EditFieldLabel.Enable = 'off';
                    app.Parameter1EditFieldLabel.Text = 'alpha';
                    app.Parameter2EditFieldLabel.Text = 'beta';
                    app.Parameter3EditFieldLabel.Text = '';
                    app.Label.Text = 'scale, alpha > 0';
                    app.Label_2.Text = 'shape, beta > 0';
                    app.Label_3.Text = '';
                    app.cfXParameter1 = str2double(app.Parameter1EditField.Value);
                    app.cfXParameter2 = str2double(app.Parameter2EditField.Value);
                    % alpha = app.cfXParameter1
                    % beta = app.cfXParameter2
                    app.cfX = @(t) cfX_Weibull(t,app.cfXParameter1,app.cfXParameter2);
                case app.SeverityDropDown.Items{13}
                    if isempty(app.Parameter1EditField.Value)
                        app.Parameter1EditField.Value = 'dlmread(''SeverityData.txt'')';
                    end
                    app.Parameter2EditField.Value = '';
                    app.Parameter3EditField.Value = '';
                    app.Parameter1EditField.Enable = 'on';
                    app.Parameter2EditField.Enable = 'off';
                    app.Parameter3EditField.Enable = 'off';
                    app.Parameter1EditFieldLabel.Enable = 'on';
                    app.Parameter2EditFieldLabel.Enable = 'off';
                    app.Parameter3EditFieldLabel.Enable = 'off';
                    app.Parameter1EditFieldLabel.Text = 'X data';
                    app.Parameter2EditFieldLabel.Text = '';
                    app.Parameter3EditFieldLabel.Text = '';
                    app.Label.Text = 'set / load the data';
                    app.Label_2.Text = '';
                    app.Label_3.Text = '';
                    app.cfXParameter1 = eval(app.Parameter1EditField.Value);
                    app.cfX = @(t) cfE_Empirical(t,app.cfXParameter1);
                case app.SeverityDropDown.Items{14}
                    app.Parameter1EditField.Value = '';
                    app.Parameter2EditField.Value = '';
                    app.Parameter3EditField.Value = '';
                    app.Label.Text = '';
                    app.Label_2.Text = '';
                    app.Label_3.Text = '';
                    app.cfX = [];
            end
            %
            value2 = app.FrequencyDropDown.Value;
            if ~strcmp(value2,app.FrequencyID)
                app.Parameter1EditField_2.Value = '';
                app.Parameter2EditField_2.Value = '';
                app.Parameter3EditField_2.Value = '';
            end
            app.FrequencyID = value2;
            app.Parameter1EditField_2.Enable = 'off';
            app.Parameter2EditField_2.Enable = 'off';
            app.Parameter3EditField_2.Enable = 'off';
            app.Parameter1Label.Enable = 'off';
            app.Parameter2Label.Enable = 'off';
            app.Parameter3Label.Enable = 'off';
            app.Parameter1Label.Text = '';
            app.Parameter2Label.Text = '';
            app.Parameter3Label.Text = '';
            switch value2
                case app.FrequencyDropDown.Items{1}
                    if isempty(app.Parameter1EditField_2.Value)
                        app.Parameter1EditField_2.Value = '10';
                    end
                    if isempty(app.Parameter2EditField_2.Value)
                        app.Parameter2EditField_2.Value = '0.5';
                    end
                    app.Parameter3EditField_2.Value = '';
                    app.Parameter1EditField_2.Enable = 'on';
                    app.Parameter2EditField_2.Enable = 'on';
                    app.Parameter3EditField_2.Enable = 'off';
                    app.Parameter1Label.Enable = 'on';
                    app.Parameter2Label.Enable = 'on';
                    app.Parameter3Label.Enable = 'off';
                    app.Parameter1Label.Text = 'n';
                    app.Parameter2Label.Text = 'p';
                    app.Parameter3Label.Text = '';
                    app.Label_4.Text = 'number of trials, n = 1,2,...';
                    app.Label_5.Text = 'success probability, p in [0,1]';
                    app.Label_6.Text = '';
                    app.cfNParameter1 = str2double(app.Parameter1EditField_2.Value);
                    app.cfNParameter2 = str2double(app.Parameter2EditField_2.Value);
                    % n = app.cfNParameter1
                    % p = app.cfNParameter2
                    app.cf = @(t) cfN_Binomial(t,app.cfNParameter1,app.cfNParameter2,app.cfX);
                case app.FrequencyDropDown.Items{2}
                    if isempty(app.Parameter1EditField_2.Value)
                        app.Parameter1EditField_2.Value = '1';
                    end
                    if isempty(app.Parameter2EditField_2.Value)
                        app.Parameter2EditField_2.Value = '1';
                    end
                    if isempty(app.Parameter3EditField_2.Value)
                        app.Parameter3EditField_2.Value = '1';
                    end
                    app.Parameter1EditField_2.Enable = 'on';
                    app.Parameter2EditField_2.Enable = 'on';
                    app.Parameter3EditField_2.Enable = 'on';
                    app.Parameter1Label.Enable = 'on';
                    app.Parameter2Label.Enable = 'on';
                    app.Parameter3Label.Enable = 'on';
                    app.Parameter1Label.Text = 'a';
                    app.Parameter2Label.Text = 'b';
                    app.Parameter3Label.Text = 'c';
                    app.Label_4.Text = 'variable mean parameter, a > 0';
                    app.Label_5.Text = 'variable mean parameter, b > 0';
                    app.Label_6.Text = 'fixed mean parameter, c > 0';
                    app.cfNParameter1 = str2double(app.Parameter1EditField_2.Value);
                    app.cfNParameter2 = str2double(app.Parameter2EditField_2.Value);
                    app.cfNParameter3 = str2double(app.Parameter3EditField_2.Value);
                    % a = app.cfNParameter1
                    % b = app.cfNParameter2
                    % c = app.cfNParameter2
                    app.cf = @(t) cfN_Delaporte(t,app.cfNParameter1,app.cfNParameter2,app.cfNParameter3,app.cfX);
                case app.FrequencyDropDown.Items{3}
                    if isempty(app.Parameter1EditField_2.Value)
                        app.Parameter1EditField_2.Value = '3';
                    end
                    if isempty(app.Parameter2EditField_2.Value)
                        app.Parameter2EditField_2.Value = '0.5';
                    end
                    app.Parameter3EditField_2.Value = '';
                    app.Parameter1EditField_2.Enable = 'on';
                    app.Parameter2EditField_2.Enable = 'on';
                    app.Parameter3EditField_2.Enable = 'off';
                    app.Parameter1Label.Enable = 'on';
                    app.Parameter2Label.Enable = 'on';
                    app.Parameter3Label.Enable = 'off';
                    app.Parameter1Label.Text = 'a';
                    app.Parameter2Label.Text = 'p';
                    app.Parameter3Label.Text = '';
                    app.Label_4.Text = 'variable mean parameter, a > 0';
                    app.Label_5.Text = 'success probability, p in [0,1]';
                    app.Label_6.Text = '';
                    app.cfNParameter1 = str2double(app.Parameter1EditField_2.Value);
                    app.cfNParameter2 = str2double(app.Parameter2EditField_2.Value);
                    % a = app.cfNParameter1
                    % p = app.cfNParameter2
                    app.cf = @(t) cfN_GeneralizedPoisson(t,app.cfNParameter1,app.cfNParameter2,app.cfX);
                case app.FrequencyDropDown.Items{4}
                    if isempty(app.Parameter1EditField_2.Value)
                        app.Parameter1EditField_2.Value = '0.5';
                    end
                    if isempty(app.Parameter2EditField_2.Value)
                        app.Parameter2EditField_2.Value = 'standard';
                    end                    
                    app.Parameter3EditField_2.Value = '';
                    app.Parameter1EditField_2.Enable = 'on';
                    app.Parameter2EditField_2.Enable = 'on';
                    app.Parameter3EditField_2.Enable = 'off';
                    app.Parameter1Label.Enable = 'on';
                    app.Parameter2Label.Enable = 'on';
                    app.Parameter3Label.Enable = 'off';
                    app.Parameter1Label.Text = 'p';
                    app.Parameter2Label.Text = 'type';
                    app.Parameter3Label.Text = '';
                    app.Label_4.Text = 'probability, p in [0,1]';
                    app.Label_5.Text = 'type (''standard'' or ''shifted'')';
                    app.Label_6.Text = '';
                    app.cfNParameter1 = str2double(app.Parameter1EditField_2.Value);
                    app.cfNParameter2 = app.Parameter2EditField_2.Value;
                    % p = app.cfNParameter1
                    app.cf = @(t) cfN_Geometric(t,app.cfNParameter1,app.cfNParameter2,app.cfX);
                case app.FrequencyDropDown.Items{5}
                    if isempty(app.Parameter1EditField_2.Value)
                        app.Parameter1EditField_2.Value = '0.5';
                    end
                    app.Parameter2EditField_2.Value = '';
                    app.Parameter3EditField_2.Value = '';
                    app.Parameter1EditField_2.Enable = 'on';
                    app.Parameter2EditField_2.Enable = 'off';
                    app.Parameter3EditField_2.Enable = 'off';
                    app.Parameter1Label.Enable = 'on';
                    app.Parameter2Label.Enable = 'off';
                    app.Parameter3Label.Enable = 'off';
                    app.Parameter1Label.Text = 'p';
                    app.Parameter2Label.Text = '';
                    app.Parameter3Label.Text = '';
                    app.Label_4.Text = 'probability, p in [0,1]';
                    app.Label_5.Text = '';
                    app.Label_6.Text = '';
                    app.cfNParameter1 = str2double(app.Parameter1EditField_2.Value);
                    % p = app.cfNParameter1
                    app.cf = @(t) cfN_Logarithmic(t,app.cfNParameter1,app.cfX);
                case app.FrequencyDropDown.Items{6}
                    if isempty(app.Parameter1EditField_2.Value)
                        app.Parameter1EditField_2.Value = '10';
                    end
                    if isempty(app.Parameter2EditField_2.Value)
                        app.Parameter2EditField_2.Value = '0.5';
                    end
                    app.Parameter3EditField_2.Value = '';
                    app.Parameter1EditField_2.Enable = 'on';
                    app.Parameter2EditField_2.Enable = 'on';
                    app.Parameter3EditField_2.Enable = 'off';
                    app.Parameter1Label.Enable = 'on';
                    app.Parameter2Label.Enable = 'on';
                    app.Parameter3Label.Enable = 'off';
                    app.Parameter1Label.Text = 'r';
                    app.Parameter2Label.Text = 'p';
                    app.Parameter3Label.Text = '';
                    app.Label_4.Text = 'number of failures, r = 1,2,...';
                    app.Label_5.Text = 'success probability, p in [0,1]';
                    app.Label_6.Text = '';
                    app.cfNParameter1 = str2double(app.Parameter1EditField_2.Value);
                    app.cfNParameter2 = str2double(app.Parameter2EditField_2.Value);
                    % r = app.cfNParameter1
                    % p = app.cfNParameter2
                    app.cf = @(t) cfN_NegativeBinomial(t,app.cfNParameter1,app.cfNParameter2,app.cfX);
                case app.FrequencyDropDown.Items{7}
                    if isempty(app.Parameter1EditField_2.Value)
                        app.Parameter1EditField_2.Value = '10';
                    end
                    app.Parameter2EditField_2.Value = '';
                    app.Parameter3EditField_2.Value = '';
                    app.Parameter1EditField_2.Enable = 'on';
                    app.Parameter2EditField_2.Enable = 'off';
                    app.Parameter3EditField_2.Enable = 'off';
                    app.Parameter1Label.Enable = 'on';
                    app.Parameter2Label.Enable = 'off';
                    app.Parameter3Label.Enable = 'off';
                    app.Parameter1Label.Text = 'lambda';
                    app.Parameter2Label.Text = '';
                    app.Parameter3Label.Text = '';
                    app.Label_4.Text = 'rate parameter, lambda > 0';
                    app.Label_5.Text = '';
                    app.Label_6.Text = '';
                    app.cfNParameter1 = str2double(app.Parameter1EditField_2.Value);
                    % lambda = app.cfNParameter1
                    app.cf = @(t) cfN_Poisson(t,app.cfNParameter1,app.cfX);
                case app.FrequencyDropDown.Items{8}
                    if isempty(app.Parameter1EditField_2.Value)
                        app.Parameter1EditField_2.Value = '1';
                    end
                    if isempty(app.Parameter2EditField_2.Value)
                        app.Parameter2EditField_2.Value = '1';
                    end
                    if isempty(app.Parameter3EditField_2.Value)
                        app.Parameter3EditField_2.Value = '5';
                    end
                    app.Parameter1EditField_2.Enable = 'on';
                    app.Parameter2EditField_2.Enable = 'on';
                    app.Parameter3EditField_2.Enable = 'on';
                    app.Parameter1Label.Enable = 'on';
                    app.Parameter2Label.Enable = 'on';
                    app.Parameter3Label.Enable = 'on';
                    app.Parameter1Label.Text = 'a';
                    app.Parameter2Label.Text = 'b';
                    app.Parameter3Label.Text = 'm';
                    app.Label_4.Text = 'a parameter, (real)';
                    app.Label_5.Text = 'b parameter, (real)';
                    app.Label_6.Text = 'm parameter, m = 1,2,...';
                    app.cfNParameter1 = str2double(app.Parameter1EditField_2.Value);
                    app.cfNParameter2 = str2double(app.Parameter2EditField_2.Value);
                    app.cfNParameter3 = str2double(app.Parameter3EditField_2.Value);
                    % a = app.cfNParameter1
                    % b = app.cfNParameter2
                    % m = app.cfNParameter2
                    app.cf = @(t) cfN_PolyaEggenberger(t,app.cfNParameter1,app.cfNParameter2,app.cfNParameter3,app.cfX);
                case app.FrequencyDropDown.Items{9}
                    if isempty(app.Parameter1EditField_2.Value)
                        app.Parameter1EditField_2.Value = '5';
                    end
                    if isempty(app.Parameter2EditField_2.Value)
                        app.Parameter2EditField_2.Value = '0.5';
                    end
                    app.Parameter3EditField_2.Value = '';
                    app.Parameter1EditField_2.Enable = 'on';
                    app.Parameter2EditField_2.Enable = 'on';
                    app.Parameter3EditField_2.Enable = 'off';
                    app.Parameter1Label.Enable = 'on';
                    app.Parameter2Label.Enable = 'on';
                    app.Parameter3Label.Enable = 'off';
                    app.Parameter1Label.Text = 'a';
                    app.Parameter2Label.Text = 'b';
                    app.Parameter3Label.Text = '';
                    app.Label_4.Text = 'a parameter, a > 0';
                    app.Label_5.Text = 'b parameter, b > 0';
                    app.Label_6.Text = '';
                    app.cfNParameter1 = str2double(app.Parameter1EditField_2.Value);
                    app.cfNParameter2 = str2double(app.Parameter2EditField_2.Value);
                    % a = app.cfNParameter1
                    % b = app.cfNParameter2
                    app.cf = @(t) cfN_Quinkert(t,app.cfNParameter1,app.cfNParameter2,app.cfX);
                case app.FrequencyDropDown.Items{10}
                    if isempty(app.Parameter1EditField_2.Value)
                        app.Parameter1EditField_2.Value = '2';
                    end
                    if isempty(app.Parameter2EditField_2.Value)
                        app.Parameter2EditField_2.Value = '2';
                    end
                    if isempty(app.Parameter3EditField_2.Value)
                        app.Parameter3EditField_2.Value = '5';
                    end
                    app.Parameter1EditField_2.Enable = 'on';
                    app.Parameter2EditField_2.Enable = 'on';
                    app.Parameter3EditField_2.Enable = 'on';
                    app.Parameter1Label.Enable = 'on';
                    app.Parameter2Label.Enable = 'on';
                    app.Parameter3Label.Enable = 'on';
                    app.Parameter1Label.Text = 'a';
                    app.Parameter2Label.Text = 'b';
                    app.Parameter3Label.Text = 'r';
                    app.Label_4.Text = 'a parameter, a > 0';
                    app.Label_5.Text = 'b parameter, b > 0';
                    app.Label_6.Text = 'r parameter, r > 0';
                    app.cfNParameter1 = str2double(app.Parameter1EditField_2.Value);
                    app.cfNParameter2 = str2double(app.Parameter2EditField_2.Value);
                    app.cfNParameter3 = str2double(app.Parameter3EditField_2.Value);
                    % a = app.cfNParameter1
                    % b = app.cfNParameter2
                    % r = app.cfNParameter2
                    if isempty(app.cfX)
                        warning('VW:CRMTool:CfDidNotConverge','CF of the Waring distribution DO NOT converge!');
                    end
                    app.cf = @(t) cfN_Waring(t,app.cfNParameter1,app.cfNParameter2,app.cfNParameter3,app.cfX);
                case app.FrequencyDropDown.Items{11}
                    if isempty(app.Parameter1EditField_2.Value)
                        app.Parameter1EditField_2.Value = 'dlmread(''FrequencyData.txt'')';
                    end
                    app.Parameter2EditField_2.Value = '';
                    app.Parameter3EditField_2.Value = '';
                    app.Parameter1EditField_2.Enable = 'on';
                    app.Parameter2EditField_2.Enable = 'off';
                    app.Parameter3EditField_2.Enable = 'off';
                    app.Parameter1Label.Enable = 'on';
                    app.Parameter2Label.Enable = 'off';
                    app.Parameter3Label.Enable = 'off';
                    app.Parameter1Label.Text = 'N data';
                    app.Parameter2Label.Text = '';
                    app.Parameter3Label.Text = '';
                    app.Label_4.Text = 'set / load the data';
                    app.Label_5.Text = '';
                    app.Label_6.Text = '';
                    app.cfNParameter1 = eval(app.Parameter1EditField_2.Value);
                    app.cf = @(t) cfE_Empirical(t,app.cfNParameter1,app.cfX);
                case app.FrequencyDropDown.Items{12}
                    app.Parameter1EditField_2.Value = '';
                    app.Parameter2EditField_2.Value = '';
                    app.Parameter3EditField_2.Value = '';
                    app.Label_4.Text = '';
                    app.Label_5.Text = '';
                    app.Label_6.Text = '';
                    if isempty(app.cfX)
                        app.cf = [];
                    else
                        app.cf = app.cfX;
                    end
            end
            %
            value3 = app.AlgorithmDropDown.Value;
            switch value3
                case app.AlgorithmDropDown.Items{2}
                    app.VaRprobabilitiesEditField.Value = '[0.9, 0.95, 0.99]';
                    app.VaRprobabilitiesEditField.Enable = 'off';
                    app.VaRprobabilitiesEditFieldLabel.Enable = 'off';
                    app.XvaluesEditField.Value = 'auto';
                    app.XvaluesEditField.Enable = 'off';
                    app.XvaluesLabel.Enable = 'off';
                    app.NEditField.Value = '2^14';
                    app.NEditField.Enable = 'off';
                    app.NLabel.Enable = 'off';
                    app.SixSigmaRuleEditField.Value = '10';
                    app.SixSigmaRuleEditField.Enable = 'off';
                    app.SixSigmaRuleLabel.Enable = 'off';
                    app.Result = [];
                otherwise
                    app.VaRprobabilitiesEditField.Enable = 'on';
                    app.VaRprobabilitiesEditFieldLabel.Enable = 'on';
                    app.XvaluesLabel.Enable = 'on';
                    app.XvaluesEditField.Enable = 'on';
                    app.SixSigmaRuleEditField.Enable = 'on';
                    app.SixSigmaRuleLabel.Enable = 'on';
                    app.NEditField.Enable = 'on';
                    app.NLabel.Enable = 'on';
                    app.options.isCompound = 1;
                    app.options.isPlot = 0;
                    app.options.xN = 501;
                    app.options.N = str2num(app.NEditField.Value);
                    app.options.SixSigmaRule = str2double(app.SixSigmaRuleEditField.Value);
                    if strcmp(app.XvaluesEditField.Value,'auto')
                        app.x = [];
                    else
                        app.x = str2num(app.XvaluesEditField.Value);
                    end
                    if strcmp(app.VaRprobabilitiesEditField.Value,'auto')
                        app.prob = [];
                    else
                        app.prob = str2num(app.VaRprobabilitiesEditField.Value);
                    end                
            end
        end

        % Button pushed function: CALCULATEButton
        function Calculate(app, event)
            % (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
            % Ver.: 30-Jun-2018 11:40:35
            value = app.AlgorithmDropDown.Value;
            ValueChanged(app, event)
            switch value
                case app.AlgorithmDropDown.Items{1}
                    app.Result = cf2DistGP(app.cf,app.x,app.prob,app.options);
            %    case app.AlgorithmDropDown.Items{2}
            %        app.Result = cf2DistFFT(app.cf,app.x,app.prob,app.options);
                case app.AlgorithmDropDown.Items{2}
                    app.Result = [];
            end
            if ~isempty(app.Result)
                app.VaRquantilesLabel.Enable = 'on';
                app.VaRquantilesEditField.Enable = 'on';
                app.VaRquantilesEditField.Value = num2str(app.Result.qf);
                if ~isempty(app.cf)
                    T = ceil(32*app.Result.dt*10)/10;
                    tt = linspace(-T,T,501);
                    cft = app.cf(tt);
                    plot(app.UIAxes3,tt,real(cft),tt,imag(cft),'LineWidth',1.0);
                    plot(app.UIAxes,app.Result.x,app.Result.cdf,'LineWidth',1.0);
                    plot(app.UIAxes2,app.Result.x,app.Result.pdf,'LineWidth',1.0);
                end
                res = app.Result;
                if app.Table
                    T = table(app.Result.x(:),app.Result.pdf(:),app.Result.cdf(:),'VariableNames',{'x' 'PDF' 'CDF'});
                    writetable(T,'CRMToolTAB.xls')
                    if numel(app.Result.x) <= 101
                        disp(T)
                    end
                end
            else
                warning('VW:CRMTool:NothingToDo','NOTHING TO DO! Select the proper numerical inversion ALGORITHM!');
            end
            if app.Save
                res.cf  = app.cf;
                res.cfX = app.cfX;
                res.cfXParameter1 = app.cfXParameter1;
                res.cfXParameter2 = app.cfXParameter2;
                res.cfXParameter3 = app.cfXParameter3;
                res.cfNParameter1 = app.cfNParameter1;
                res.cfNParameter2 = app.cfNParameter2;
                res.cfNParameter3 = app.cfNParameter3;
                save CRMresult res;
            end
            clear res;
        end

        % Button pushed function: EXITButton
        function ExitFcn(app, event)
            % (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
            % Ver.: 30-Jun-2018 11:40:35
            delete(app)
            %warning('on','MATLAB:subscripting:noSubscriptsSpecified')
            %warning('on','MATLAB:Figure:FigureSavedToMATFile')
            %warning('on','MATLAB:ui:uifigure:UnsupportedAppDesignerFunctionality')
            %warning('on','MATLAB:lang:cannotClearExecutingFunction')
            warning('on','VW:CharFunTool:cf2DistGP')
        end

        % Value changed function: SaveSwitch
        function SaveFcn(app, event)
            % (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
            % Ver.: 30-Jun-2018 11:40:35
            value = app.SaveSwitch.Value;
            if strcmp(value,'On')
                app.Save = 1;
            else
                app.Save = 0;
            end
        end

        % Value changed function: TableSwitch
        function TableFcn(app, event)
            % (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
            % Ver.: 30-Jun-2018 11:40:35
            value = app.TableSwitch.Value;
            if strcmp(value,'On')
                app.Table = 1;
            else
                app.Table = 0;
            end
        end

        % Button pushed function: AboutButton
        function AboutFcn(app, event)
            % (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
            % Ver.: 30-Jun-2018 11:40:35
            About = [...
                'CRMTool is an Aggregate Loss Distribution Calculator, suggested to            ';
                'evaluate the Collective Risk Model (CRM) compound Aggregate Loss              ';
                'Distribution (ALD) and the associated Value at Risk (VaR / quantiles)         ';
                'by numerical inversion of its characteristic function (CF).                   ';
                '                                                                              ';
                'CRMTool is a part of the CharFunTool (The Characteristic Functions            ';
                'Toolbox) which consists of a set of algorithms for evaluating selected        ';
                'characteristic funcions and algorithms for numerical inversion of the         ';
                '(combined and/or compound) characteristic functions, used to evaluate the     ';
                'probability density function (PDF) and the cumulative distribution            ';
                'function (CDF).                                                               ';
                '                                                                              ';
                'The inversion algorithms are based on using simple trapezoidal rule           ';
                'for computing the integrals defined by the Gil-Pelaez formulae.               ';
                '                                                                              ';
                'CRMTool:   Aggregate Loss Distribution Calculator                             ';
                'Version:   2.0. (this version of CRMTool is part of CharFunTool)              ';
                'Date:      30-Jun-2018 11:40:35                                               ';
                'Copyright: (c) 2018 Viktor Witkovsky, Bratislava (witkovsky@gmail.com)        ';
                '                                                                              ';
                'All rights reserved.                                                          ';
                '                                                                              ';
                '    Redistribution and use in source and binary forms, with or without        ';
                '    modification, are permitted provided that the following conditions are    ';
                '    met:                                                                      ';
                '                                                                              ';
                '    * Redistributions of source code must retain the above copyright          ';
                '    notice, this list of conditions and the following disclaimer.             ';
                '    * Redistributions in binary form must reproduce the above copyright       ';
                '    notice, this list of conditions and the following disclaimer in the       ';
                '    documentation and/or other materials provided with the distribution       ';
                '                                                                              ';
                '    DISCLAIMER THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         ';
                '    CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,    ';
                '    BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS ';
                '    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT  ';
                '    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     ';
                '    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  ';
                '    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    ';
                '    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    ';
                '    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      ';
                '    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        ';
                '    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              '];
            disp(About)
            
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create CRMToolUIFigure
            app.CRMToolUIFigure = uifigure;
            app.CRMToolUIFigure.Color = [0.9412 0.9412 0.9412];
            app.CRMToolUIFigure.Position = [518 52 1000 760];
            app.CRMToolUIFigure.Name = 'CRMTool';

            % Create CollectiveRiskModelLabel
            app.CollectiveRiskModelLabel = uilabel(app.CRMToolUIFigure);
            app.CollectiveRiskModelLabel.VerticalAlignment = 'top';
            app.CollectiveRiskModelLabel.FontSize = 22;
            app.CollectiveRiskModelLabel.FontWeight = 'bold';
            app.CollectiveRiskModelLabel.FontColor = [0.7059 0.2118 0.2118];
            app.CollectiveRiskModelLabel.Position = [25 710 234 29];
            app.CollectiveRiskModelLabel.Text = 'Collective Risk Model';

            % Create SeverityDistributionPanel
            app.SeverityDistributionPanel = uipanel(app.CRMToolUIFigure);
            app.SeverityDistributionPanel.Title = 'Severity Distribution';
            app.SeverityDistributionPanel.Position = [25 601 320 100];

            % Create SeverityLabel
            app.SeverityLabel = uilabel(app.SeverityDistributionPanel);
            app.SeverityLabel.HorizontalAlignment = 'right';
            app.SeverityLabel.VerticalAlignment = 'top';
            app.SeverityLabel.Position = [19 47 52 15];
            app.SeverityLabel.Text = 'Severity:';

            % Create SeverityDropDown
            app.SeverityDropDown = uidropdown(app.SeverityDistributionPanel);
            app.SeverityDropDown.Items = {'Exponential', 'Fisher Snedecor F', 'Gamma', 'Generalized Pareto', 'Chi Squared', 'Inverse Gamma', 'Log Logistic', 'Log Normal', 'Pareto', 'Pearson Type V', 'Pearson Type VI', 'Weibull', 'Empirical', 'None'};
            app.SeverityDropDown.ValueChangedFcn = createCallbackFcn(app, @ValueChanged, true);
            app.SeverityDropDown.Position = [125 43 173 22];
            app.SeverityDropDown.Value = 'None';

            % Create FrequencyDistributionPanel
            app.FrequencyDistributionPanel = uipanel(app.CRMToolUIFigure);
            app.FrequencyDistributionPanel.Title = 'Frequency Distribution';
            app.FrequencyDistributionPanel.Position = [25 493 320 100];

            % Create FrequencyLabel
            app.FrequencyLabel = uilabel(app.FrequencyDistributionPanel);
            app.FrequencyLabel.HorizontalAlignment = 'right';
            app.FrequencyLabel.VerticalAlignment = 'top';
            app.FrequencyLabel.Position = [17 46 66 15];
            app.FrequencyLabel.Text = 'Frequency:';

            % Create FrequencyDropDown
            app.FrequencyDropDown = uidropdown(app.FrequencyDistributionPanel);
            app.FrequencyDropDown.Items = {'Binomial', 'Delaporte', 'Generalized Poisson', 'Geometric', 'Logarithmic', 'Negative Binomial', 'Poisson', 'Polya Eggenberger', 'Quinkert', 'Waring', 'Empirical', 'None'};
            app.FrequencyDropDown.ValueChangedFcn = createCallbackFcn(app, @ValueChanged, true);
            app.FrequencyDropDown.Position = [123 42 173 22];
            app.FrequencyDropDown.Value = 'None';

            % Create NumericalInversionAlgorithmPanel
            app.NumericalInversionAlgorithmPanel = uipanel(app.CRMToolUIFigure);
            app.NumericalInversionAlgorithmPanel.Title = 'Numerical Inversion Algorithm';
            app.NumericalInversionAlgorithmPanel.Position = [25 364 320 120];

            % Create AlgorithmLabel
            app.AlgorithmLabel = uilabel(app.NumericalInversionAlgorithmPanel);
            app.AlgorithmLabel.HorizontalAlignment = 'right';
            app.AlgorithmLabel.VerticalAlignment = 'top';
            app.AlgorithmLabel.Position = [30 60 60 15];
            app.AlgorithmLabel.Text = 'Algorithm:';

            % Create AlgorithmDropDown
            app.AlgorithmDropDown = uidropdown(app.NumericalInversionAlgorithmPanel);
            app.AlgorithmDropDown.Items = {'cf2DistGP', 'None'};
            app.AlgorithmDropDown.ValueChangedFcn = createCallbackFcn(app, @ValueChanged, true);
            app.AlgorithmDropDown.Position = [126 56 173 22];
            app.AlgorithmDropDown.Value = 'None';

            % Create OptionsPanel
            app.OptionsPanel = uipanel(app.CRMToolUIFigure);
            app.OptionsPanel.Title = 'Options';
            app.OptionsPanel.Position = [655 364 321 120];

            % Create NLabel
            app.NLabel = uilabel(app.OptionsPanel);
            app.NLabel.HorizontalAlignment = 'right';
            app.NLabel.VerticalAlignment = 'top';
            app.NLabel.Enable = 'off';
            app.NLabel.Position = [129 60 25 15];
            app.NLabel.Text = 'N:';

            % Create NEditField
            app.NEditField = uieditfield(app.OptionsPanel, 'text');
            app.NEditField.Enable = 'off';
            app.NEditField.Position = [169 56 146 22];
            app.NEditField.Value = '2^14';

            % Create SixSigmaRuleLabel
            app.SixSigmaRuleLabel = uilabel(app.OptionsPanel);
            app.SixSigmaRuleLabel.HorizontalAlignment = 'right';
            app.SixSigmaRuleLabel.VerticalAlignment = 'top';
            app.SixSigmaRuleLabel.Enable = 'off';
            app.SixSigmaRuleLabel.Position = [61 22 94 15];
            app.SixSigmaRuleLabel.Text = 'Six-Sigma-Rule:';

            % Create SixSigmaRuleEditField
            app.SixSigmaRuleEditField = uieditfield(app.OptionsPanel, 'text');
            app.SixSigmaRuleEditField.Enable = 'off';
            app.SixSigmaRuleEditField.Position = [170 18 145 22];
            app.SixSigmaRuleEditField.Value = '10';

            % Create SeverityDistributionParametersPanel
            app.SeverityDistributionParametersPanel = uipanel(app.CRMToolUIFigure);
            app.SeverityDistributionParametersPanel.Title = 'Severity Distribution Parameters';
            app.SeverityDistributionParametersPanel.Position = [345 601 630 100];

            % Create Parameter2EditFieldLabel
            app.Parameter2EditFieldLabel = uilabel(app.SeverityDistributionParametersPanel);
            app.Parameter2EditFieldLabel.HorizontalAlignment = 'right';
            app.Parameter2EditFieldLabel.VerticalAlignment = 'top';
            app.Parameter2EditFieldLabel.Enable = 'off';
            app.Parameter2EditFieldLabel.Position = [214 46 72 15];
            app.Parameter2EditFieldLabel.Text = 'Parameter 2';

            % Create Parameter2EditField
            app.Parameter2EditField = uieditfield(app.SeverityDistributionParametersPanel, 'text');
            app.Parameter2EditField.Enable = 'off';
            app.Parameter2EditField.Position = [301 42 115 22];

            % Create Parameter3EditFieldLabel
            app.Parameter3EditFieldLabel = uilabel(app.SeverityDistributionParametersPanel);
            app.Parameter3EditFieldLabel.HorizontalAlignment = 'right';
            app.Parameter3EditFieldLabel.VerticalAlignment = 'top';
            app.Parameter3EditFieldLabel.Enable = 'off';
            app.Parameter3EditFieldLabel.Position = [422 47 72 15];
            app.Parameter3EditFieldLabel.Text = 'Parameter 3';

            % Create Parameter3EditField
            app.Parameter3EditField = uieditfield(app.SeverityDistributionParametersPanel, 'text');
            app.Parameter3EditField.Enable = 'off';
            app.Parameter3EditField.Position = [509 43 115 22];

            % Create Label
            app.Label = uilabel(app.SeverityDistributionParametersPanel);
            app.Label.HorizontalAlignment = 'right';
            app.Label.VerticalAlignment = 'top';
            app.Label.Enable = 'off';
            app.Label.Position = [0 16 196 15];
            app.Label.Text = '';

            % Create Label_2
            app.Label_2 = uilabel(app.SeverityDistributionParametersPanel);
            app.Label_2.HorizontalAlignment = 'right';
            app.Label_2.VerticalAlignment = 'top';
            app.Label_2.Enable = 'off';
            app.Label_2.Position = [216 16 201 15];
            app.Label_2.Text = '';

            % Create Label_3
            app.Label_3 = uilabel(app.SeverityDistributionParametersPanel);
            app.Label_3.HorizontalAlignment = 'right';
            app.Label_3.VerticalAlignment = 'top';
            app.Label_3.Enable = 'off';
            app.Label_3.Position = [432 16 192 15];
            app.Label_3.Text = '';

            % Create Parameter1EditFieldLabel
            app.Parameter1EditFieldLabel = uilabel(app.SeverityDistributionParametersPanel);
            app.Parameter1EditFieldLabel.HorizontalAlignment = 'right';
            app.Parameter1EditFieldLabel.VerticalAlignment = 'top';
            app.Parameter1EditFieldLabel.Enable = 'off';
            app.Parameter1EditFieldLabel.Position = [0 46 72 15];
            app.Parameter1EditFieldLabel.Text = 'Parameter 1';

            % Create Parameter1EditField
            app.Parameter1EditField = uieditfield(app.SeverityDistributionParametersPanel, 'text');
            app.Parameter1EditField.Enable = 'off';
            app.Parameter1EditField.Position = [85 42 115 22];

            % Create FrequencyDistributionParametersPanel
            app.FrequencyDistributionParametersPanel = uipanel(app.CRMToolUIFigure);
            app.FrequencyDistributionParametersPanel.Title = 'Frequency Distribution Parameters';
            app.FrequencyDistributionParametersPanel.Position = [345 493 630 100];

            % Create Parameter1Label
            app.Parameter1Label = uilabel(app.FrequencyDistributionParametersPanel);
            app.Parameter1Label.HorizontalAlignment = 'right';
            app.Parameter1Label.VerticalAlignment = 'top';
            app.Parameter1Label.Enable = 'off';
            app.Parameter1Label.Position = [0 46 75 15];
            app.Parameter1Label.Text = ' Parameter 1';

            % Create Parameter1EditField_2
            app.Parameter1EditField_2 = uieditfield(app.FrequencyDistributionParametersPanel, 'text');
            app.Parameter1EditField_2.Enable = 'off';
            app.Parameter1EditField_2.Position = [81 42 119 22];

            % Create Parameter2Label
            app.Parameter2Label = uilabel(app.FrequencyDistributionParametersPanel);
            app.Parameter2Label.HorizontalAlignment = 'right';
            app.Parameter2Label.VerticalAlignment = 'top';
            app.Parameter2Label.Enable = 'off';
            app.Parameter2Label.Position = [214 50 72 15];
            app.Parameter2Label.Text = 'Parameter 2';

            % Create Parameter2EditField_2
            app.Parameter2EditField_2 = uieditfield(app.FrequencyDistributionParametersPanel, 'text');
            app.Parameter2EditField_2.Enable = 'off';
            app.Parameter2EditField_2.Position = [301 46 116 22];

            % Create Parameter3Label
            app.Parameter3Label = uilabel(app.FrequencyDistributionParametersPanel);
            app.Parameter3Label.HorizontalAlignment = 'right';
            app.Parameter3Label.VerticalAlignment = 'top';
            app.Parameter3Label.Enable = 'off';
            app.Parameter3Label.Position = [422 47 72 15];
            app.Parameter3Label.Text = 'Parameter 3';

            % Create Parameter3EditField_2
            app.Parameter3EditField_2 = uieditfield(app.FrequencyDistributionParametersPanel, 'text');
            app.Parameter3EditField_2.Enable = 'off';
            app.Parameter3EditField_2.Position = [509 43 115 22];

            % Create Label_4
            app.Label_4 = uilabel(app.FrequencyDistributionParametersPanel);
            app.Label_4.HorizontalAlignment = 'right';
            app.Label_4.VerticalAlignment = 'top';
            app.Label_4.Enable = 'off';
            app.Label_4.Position = [0 13 200 15];
            app.Label_4.Text = '';

            % Create Label_5
            app.Label_5 = uilabel(app.FrequencyDistributionParametersPanel);
            app.Label_5.HorizontalAlignment = 'right';
            app.Label_5.VerticalAlignment = 'top';
            app.Label_5.Enable = 'off';
            app.Label_5.Position = [214 13 202 15];
            app.Label_5.Text = '';

            % Create Label_6
            app.Label_6 = uilabel(app.FrequencyDistributionParametersPanel);
            app.Label_6.HorizontalAlignment = 'right';
            app.Label_6.VerticalAlignment = 'top';
            app.Label_6.Enable = 'off';
            app.Label_6.Position = [422 13 202 15];
            app.Label_6.Text = '';

            % Create InputArgumentsPanel
            app.InputArgumentsPanel = uipanel(app.CRMToolUIFigure);
            app.InputArgumentsPanel.Title = 'Input Arguments';
            app.InputArgumentsPanel.Position = [345 364 311 120];

            % Create VaRprobabilitiesEditFieldLabel
            app.VaRprobabilitiesEditFieldLabel = uilabel(app.InputArgumentsPanel);
            app.VaRprobabilitiesEditFieldLabel.HorizontalAlignment = 'right';
            app.VaRprobabilitiesEditFieldLabel.VerticalAlignment = 'top';
            app.VaRprobabilitiesEditFieldLabel.Enable = 'off';
            app.VaRprobabilitiesEditFieldLabel.Position = [15 22 102 15];
            app.VaRprobabilitiesEditFieldLabel.Text = 'VaR probabilities:';

            % Create VaRprobabilitiesEditField
            app.VaRprobabilitiesEditField = uieditfield(app.InputArgumentsPanel, 'text');
            app.VaRprobabilitiesEditField.Enable = 'off';
            app.VaRprobabilitiesEditField.Position = [127 18 166 22];
            app.VaRprobabilitiesEditField.Value = '[0.9, 0.95, 0.99]';

            % Create XvaluesLabel
            app.XvaluesLabel = uilabel(app.InputArgumentsPanel);
            app.XvaluesLabel.HorizontalAlignment = 'right';
            app.XvaluesLabel.VerticalAlignment = 'top';
            app.XvaluesLabel.Enable = 'off';
            app.XvaluesLabel.Position = [17 60 55 15];
            app.XvaluesLabel.Text = 'X values:';

            % Create XvaluesEditField
            app.XvaluesEditField = uieditfield(app.InputArgumentsPanel, 'text');
            app.XvaluesEditField.Enable = 'off';
            app.XvaluesEditField.Position = [129 56 164 22];
            app.XvaluesEditField.Value = 'auto';

            % Create ResultsPanel
            app.ResultsPanel = uipanel(app.CRMToolUIFigure);
            app.ResultsPanel.Title = 'Results';
            app.ResultsPanel.Position = [27 14 949 340];

            % Create UIAxes
            app.UIAxes = uiaxes(app.ResultsPanel);
            title(app.UIAxes, 'Aggregate Loss CDF')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'CDF')
            app.UIAxes.GridAlpha = 0.15;
            app.UIAxes.MinorGridAlpha = 0.25;
            app.UIAxes.Box = 'on';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.Position = [618 80 307 185];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.ResultsPanel);
            title(app.UIAxes2, 'Aggregate Loss PDF')
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'PDF')
            app.UIAxes2.GridAlpha = 0.15;
            app.UIAxes2.MinorGridAlpha = 0.25;
            app.UIAxes2.Box = 'on';
            app.UIAxes2.XGrid = 'on';
            app.UIAxes2.YGrid = 'on';
            app.UIAxes2.Position = [318 80 300 185];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.ResultsPanel);
            title(app.UIAxes3, 'Aggregate Loss CF')
            xlabel(app.UIAxes3, 't')
            ylabel(app.UIAxes3, 'CF')
            app.UIAxes3.GridAlpha = 0.15;
            app.UIAxes3.MinorGridAlpha = 0.25;
            app.UIAxes3.Box = 'on';
            app.UIAxes3.XGrid = 'on';
            app.UIAxes3.YGrid = 'on';
            app.UIAxes3.Position = [13 80 300 185];

            % Create CALCULATEButton
            app.CALCULATEButton = uibutton(app.ResultsPanel, 'push');
            app.CALCULATEButton.ButtonPushedFcn = createCallbackFcn(app, @Calculate, true);
            app.CALCULATEButton.Position = [61 41 126 22];
            app.CALCULATEButton.Text = 'CALCULATE';

            % Create EXITButton
            app.EXITButton = uibutton(app.ResultsPanel, 'push');
            app.EXITButton.ButtonPushedFcn = createCallbackFcn(app, @ExitFcn, true);
            app.EXITButton.Position = [787 41 126 22];
            app.EXITButton.Text = 'EXIT';

            % Create SaveSwitchLabel
            app.SaveSwitchLabel = uilabel(app.ResultsPanel);
            app.SaveSwitchLabel.HorizontalAlignment = 'center';
            app.SaveSwitchLabel.VerticalAlignment = 'top';
            app.SaveSwitchLabel.Position = [550 13 33 15];
            app.SaveSwitchLabel.Text = 'Save';

            % Create SaveSwitch
            app.SaveSwitch = uiswitch(app.ResultsPanel, 'slider');
            app.SaveSwitch.ValueChangedFcn = createCallbackFcn(app, @SaveFcn, true);
            app.SaveSwitch.Position = [544 43 45 20];

            % Create VaRquantilesLabel
            app.VaRquantilesLabel = uilabel(app.ResultsPanel);
            app.VaRquantilesLabel.HorizontalAlignment = 'right';
            app.VaRquantilesLabel.VerticalAlignment = 'top';
            app.VaRquantilesLabel.Enable = 'off';
            app.VaRquantilesLabel.Position = [28 287 85 15];
            app.VaRquantilesLabel.Text = 'VaR quantiles:';

            % Create VaRquantilesEditField
            app.VaRquantilesEditField = uieditfield(app.ResultsPanel, 'text');
            app.VaRquantilesEditField.Enable = 'off';
            app.VaRquantilesEditField.Position = [120 283 823 22];

            % Create TableSwitchLabel
            app.TableSwitchLabel = uilabel(app.ResultsPanel);
            app.TableSwitchLabel.HorizontalAlignment = 'center';
            app.TableSwitchLabel.VerticalAlignment = 'top';
            app.TableSwitchLabel.Position = [391.5 13 36 15];
            app.TableSwitchLabel.Text = 'Table';

            % Create TableSwitch
            app.TableSwitch = uiswitch(app.ResultsPanel, 'slider');
            app.TableSwitch.ValueChangedFcn = createCallbackFcn(app, @TableFcn, true);
            app.TableSwitch.Position = [387 43 45 20];

            % Create AggregateLossDistributionCalculatorLabel
            app.AggregateLossDistributionCalculatorLabel = uilabel(app.CRMToolUIFigure);
            app.AggregateLossDistributionCalculatorLabel.HorizontalAlignment = 'right';
            app.AggregateLossDistributionCalculatorLabel.VerticalAlignment = 'top';
            app.AggregateLossDistributionCalculatorLabel.FontSize = 24;
            app.AggregateLossDistributionCalculatorLabel.FontWeight = 'bold';
            app.AggregateLossDistributionCalculatorLabel.FontColor = [0.5255 0.5255 0.5255];
            app.AggregateLossDistributionCalculatorLabel.Position = [394 708 460 31];
            app.AggregateLossDistributionCalculatorLabel.Text = 'Aggregate Loss Distribution Calculator';

            % Create AboutButton
            app.AboutButton = uibutton(app.CRMToolUIFigure, 'push');
            app.AboutButton.ButtonPushedFcn = createCallbackFcn(app, @AboutFcn, true);
            app.AboutButton.Position = [876 712 100 22];
            app.AboutButton.Text = 'About';
        end
    end

    methods (Access = public)

        % Construct app
        function app = CRMTool

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.CRMToolUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.CRMToolUIFigure)
        end
    end
end