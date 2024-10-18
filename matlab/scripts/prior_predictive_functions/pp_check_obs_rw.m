function omega2 = pp_check_obs_rw(u,n, percep,optimal)

for i=1:n
    worked =0;
    
    while ~worked
        try
            result{i}= tapas_sampleModel(u,percep);
            worked=1;
            omega2(i) = result{1, i}.p_prc.al;
        end
    end    
           
    t = ones(1,size(result{i}.u,1));
    
    ts = cumsum(t);
    ts = [0, ts];
    plot(ts, [0.5; result{1, i}.traj.v], 'c', 'LineWidth', 0.5);
    hold on
end

optimal1= tapas_sampleModel(u,optimal);
plot(ts, [0.5; optimal1.traj.v], 'r', 'LineWidth', 1);


plot(u, '.', 'Color', [0 0.6 0], 'MarkerSize', 11)
xlabel('Trial number')
ylabel('u')
axis([1, 320, -0.1, 1.1])

% tapas_hgf_binary_pu_tgi_plotTraj(ans)
