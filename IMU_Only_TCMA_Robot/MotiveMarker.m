classdef MotiveMarker < MotiveTracker
    
    methods
        function newObj = transform(obj,H)
            
            % Nested call for arrays
            N = numel(obj);
            if N > 1
                newObj(N) = MotiveMarker();
                for i = 1:N
                    newObj(i) = obj(i).transform(H);
                end
                return
            end
            
            p = obj.Position;
            nt = size(p,1);
            
            pt = [p.'; ones(1,nt)];
            newPt = H * pt;
            newP = newPt(1:3,:).';
            
            newObj = MotiveMarker;
            newObj.Name = obj.Name;
            newObj.Time = obj.Time;
            newObj.Type = obj.Type;
            newObj.Position = newP;
        end
    end
end