function kap = pp_check_kap(u,n, percep,optimal)

for i=1:n
    worked =0;
    
    while ~worked
        try
            result{i}= tapas_sampleModel(u,percep);
            worked=1;
            kap(i) = result{1,i}.p_prc.ka(2);
        end
    end    
           
    t = ones(1,size(result{i}.u,1));
    
    ts = cumsum(t);
    ts = [0, ts];
    plot(ts, [tapas_sgm(result{i}.p_prc.mu_0(2), 1); tapas_sgm(result{i}.traj.mu(:,2), 1)], 'c', 'LineWidth', 0.5);
    hold on
end

optimal1= tapas_sampleModel(u,optimal);
plot(ts, [tapas_sgm(optimal1.p_prc.mu_0(2), 1); tapas_sgm(optimal1.traj.mu(:,2), 1)], 'r', 'LineWidth', 1);


plot(u, '.', 'Color', [0 0.6 0], 'MarkerSize', 11)
xlabel('Trial number')
ylabel('u')
axis([1, 320, -0.1, 1.1])

% tapas_hgf_binary_pu_tgi_plotTraj(ans)
