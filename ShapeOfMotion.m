classdef ShapeOfMotion
    %SHAPEOFMOTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N
        d
        d_aff
        global_camera_id
        num_tests_per_sigma
        transf_end_thresh
        max_icp_iterations
        num_edges_full 
        num_edges_triangular
        num_edges 
        procrustes_mode
        riem_grad_mode
        initguess_is_available
    end
    
    methods
%         function obj = ShapeOfMotion(inputArg1,inputArg2)
%             %SHAPEOFMOTION Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
%         end

        function obj = ShapeOfMotion(params_filename)
            %SHAPEOFMOTION Construct an instance of this class filling
            % parameters according to csv file named 'params_filename'
            % positional correspondence!
            [obj.N, obj.d, obj.global_camera_id, ... 
                obj.num_tests_per_sigma, obj.transf_end_thresh, obj.max_icp_iterations, ...
                obj.procrustes_mode, obj.riem_grad_mode, obj.initguess_is_available] = ...
            readvars(params_filename);
            obj.d_aff = obj.d + 1;
            obj.num_edges_full = obj.N*obj.N - obj.N;
            obj.num_edges_triangular = 0.5 * obj.N * (obj.N-1);
            obj.num_edges = obj.num_edges_full; %FIXME: when considering non-necessarily-full graph            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

