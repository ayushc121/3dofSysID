function AC_VER2_plotting(oe_y, oe_PE, GT_y, GT_PE, u_cmd, u_act, gust, suffix)

time = GT_y(:, 5);

folder = "plots";
if ~exist(folder, 'dir')
    mkdir(folder);
end

figure
subplot(2, 2, 1);
    hold on
    plot(time, oe_y(:, 1), 'r');
    plot(time, GT_y(:, 1), 'b');
    hold off
    title('V');
    xlabel("Time, seconds")    
    ylabel("V, m/s")
    legend('Identified Model', 'Ground Truth')
subplot(2, 2, 2);
    hold on
    plot(time, oe_y(:, 2), 'r');
    plot(time, GT_y(:, 2), 'b');
    hold off
    title('Alpha');
    xlabel("Time, seconds")
    ylabel("Alpha, rad")
    legend('Identified Model', 'Ground Truth')
subplot(2, 2, 3);
    hold on
    plot(time, oe_y(:, 3), 'r');
    plot(time, GT_y(:, 3), 'b');
    hold off
    title('Gamma');
    xlabel("Time, seconds")
    ylabel("Gamma, rad")
    legend('Identified Model', 'Ground Truth')
subplot(2, 2, 4);
    hold on
    plot(time, oe_y(:, 4), 'r');
    plot(time, GT_y(:, 4), 'b');
    hold off
    title('Q');
    xlabel("Time, seconds")
    ylabel("Q, rad/s")
    legend('Identified Model', 'Ground Truth')
sgtitle('Truth and Modeled State over Time, ' + suffix)

filename = "state_comparison_" + suffix + ".png";
exportgraphics(gcf, fullfile(folder, filename), 'Resolution', 300);



figure
subplot(3, 1, 1)
    hold on
    plot(time, u_cmd(:, 1), 'r');
    plot(time, u_act(:, 1), 'b');
    hold off
    title('Commanded and Applied Thrust Inputs');
    xlabel("Time, seconds")
    ylabel("Thrust Input")
    legend('Commanded Thrust Input', 'Applied Thrust Input')
subplot(3, 1, 2)
    hold on
    plot(time, u_cmd(:, 2), 'r');
    plot(time, u_act(:, 2), 'b');
    hold off
    title('Commanded and Applied Elevator Inputs');
    xlabel("Time, seconds")
    ylabel("Elevator Input")
    legend('Commanded Elevator Input', 'Applied Elevator Input')
subplot(3, 1, 3)
    plot(time, gust, 'b');
    title('Gust');
    xlabel("Time, seconds")
    ylabel("Gust")
    legend('Truth Gust State')
sgtitle('Commanded/Applied Inputs and Gust over Time, ' + suffix)

filename = "inputs_" + suffix + ".png";
exportgraphics(gcf, fullfile(folder, filename), 'Resolution', 300);

end