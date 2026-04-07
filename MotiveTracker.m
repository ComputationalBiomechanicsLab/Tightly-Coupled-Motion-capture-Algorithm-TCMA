classdef (Abstract) MotiveTracker < matlab.mixin.Heterogeneous
    
    properties (Description = 'Data')
        Position(:,3) double
        Time
    end
    
    properties
        Name {char,string}
        Type
    end
    
    %% Sealed methods
    methods (Sealed)
        function selection = select(objs,name)
            selection = objs(contains({objs.Name},name));
        end
    end
    
    %% Protected static methods
    methods (Static, Access = protected)
        function obj = getDefaultScalarElement()
            obj = MotiveMarker();
        end
    end
end