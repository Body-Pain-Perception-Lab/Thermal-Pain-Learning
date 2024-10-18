function dn = pp_check(u,n, percep, respons,optimal)

for i=1:n
    worked =0;
    
    while ~worked
        try
            result{i}= tapas_sampleModel(u,percep,respons);
            worked=1;
            dn(i) = result{1,i}.p_obs.ze;
            resp(:,i) = result{1, i}.y;
        end
    end    
end
for i = 1:306
if mean(resp(i,:)) > 0.5
    resps(i,1) = 1
else
    resps(i,1) = 0
end
end
    t = ones(1,size(result{1}.u,1));
    
    ts = cumsum(t);
    ts = [0, ts];
optimal1= tapas_sampleModel(u,optimal);
plot(ts, [tapas_sgm(optimal1.p_prc.mu_0(2), 1); tapas_sgm(optimal1.traj.mu(:,2), 1)], 'r', 'LineWidth', 1);
hold on
plot(u, '.', 'Color', [0 0.6 0], 'MarkerSize', 11)
xlabel('Trial number')
ylabel('u')
axis([1, 320, -0.1, 1.1])
hold on
resps = resps -0.5; resps = 1.16 *resps; resps = resps +0.5
plot(resps, '.', 'Color', [1 0.7 0])
