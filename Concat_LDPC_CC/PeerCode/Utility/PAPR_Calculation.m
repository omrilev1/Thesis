function [PAPR] = PAPR_Calculation(Constellation_Mapping_Array)

    PAPR = max(abs(Constellation_Mapping_Array).^2)/mean(abs(Constellation_Mapping_Array).^2);

end

