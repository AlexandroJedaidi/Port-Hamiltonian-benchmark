function [h] = cubic_spline(s)
% cubic spline function for DWE in order to calculate q_i

h = 0;

if (0<=s) && (s<=1)
    
    h = 1 - 1.5*(s^2) + 0.75*(s^3);
    
elseif (1<s) && (s<=2)
        
    h = 0.25*(2 - s)^3;
   
end

end

