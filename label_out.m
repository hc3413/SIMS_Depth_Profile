function [] = label_out(xn)

if xn == 0
    ylabel('Counts')
    xlabel('Sputter time')
elseif xn == 1
    ylabel('Relative oxygen concentration (%)')
    xlabel('Sputter time')
else
    ylabel('Relative oxygen concentration (%)')
    xlabel('Depth (nm')
end

end

