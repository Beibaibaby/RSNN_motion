I_i=ones(10000,1);

delta_t=0.1;
v0=-70;
v_rest=-65;
v_thresh=-50;
v_reset=-65;
v_spike=40;
tau=20;

firingrate = [];
for i = linspace(0,2,100)
    I=i*I_i;
    [spike_array, spike_timestamps, potential] = lif_neuron(I, delta_t, v0, v_rest, v_thresh,v_reset, v_spike, tau);
    %disp(size(spike_timestamps));
    %disp(spike_timestamps);
    %disp(potential)
    firingrate(end+1) = size(spike_timestamps,2);
end
%[spike_array, spike_timestamps, potential] = lif_neuron(I, delta_t, v0, v_rest, v_thresh,v_reset, v_spike, tau);
%disp(spike_array);
%disp( potential);
%disp(spike_timestamps);
%plot(potential)

plot(linspace(0,2,100),firingrate)