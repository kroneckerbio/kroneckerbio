function handle = sdLinear(a, b)

handle = @sdHandle;

    function [sigma dsigmady d2sigmady2] = sdHandle(t, yInd, yVal)
        sigma = a + b .* yVal;
        
        if nargout >= 2
            % First derivative
            dsigmady = b + zeros(size(yVal));
            
            if nargout >=3
                % Second derivative
                d2sigmady2 = zeros(size(yVal));
            end
        end
    end
end