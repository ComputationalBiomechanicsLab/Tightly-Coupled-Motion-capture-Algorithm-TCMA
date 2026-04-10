classdef (Abstract) AbstractObject < handle & matlab.mixin.CustomDisplay
    % AbstractObject - Abstract class to define basic object functionality
    %
    
%     properties (Description = 'Object')
%         Name (1,:) {char,string}
%         Description {char,string}
%     end
    
    %% Constructor
    methods
        function obj = AbstractObject(varargin)
            % MyClass constructor
            if ~isempty(varargin)
                p = inputParser();
                props = properties(obj);
                for i = 1:numel(props)
                    prop = props{i};
                    p.addParameter(prop,[]);
                end
                p.parse(varargin{:});
                S = rmfield(p.Results,p.UsingDefaults);
                for f = fieldnames(S).'
                    obj.(f{:}) = S.(f{:});
                end
            end
        end
    end
    
    %% Protected methods
    methods (Access = protected)
        function pg = getPropertyGroups(obj)
            if ~isscalar(obj)
                pg = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else
                mp = metaclass(obj).PropertyList;
                mp = setdiff(mp,mp.findobj('Hidden',true));
                propDescription = {mp.Description};
                propNames = {mp.Name};
                [groupNames,~,groupIndex] = unique(propDescription,'stable');
                
                for i = 1:numel(groupNames)
                    if isempty(groupNames{i})
                        groupNames{i} = 'Other';
                    end
                    
                    groupTitle = sprintf('<strong>%s</strong>',groupNames{i});
                    pg(i) =  matlab.mixin.util.PropertyGroup('Title',groupTitle);
                    props = propNames(groupIndex == i);
                    vals = cellfun(@(p) obj.(p),props,'UniformOutput',false);
                    pg(i).PropertyList = cell2struct(vals,props,2);
                end
            end
        end
    end
    
end