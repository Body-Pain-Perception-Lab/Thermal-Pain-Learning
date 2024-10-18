function [traj, infStates] = tapas_hhmm(r, p, varargin)
% Estimates a hierarchical hidden Markov model (HHMM)
%
% This function can be called in two ways:
% 
% (1) hhmm(r, p)
%   
%     where r is the structure generated by fitModel and p is the parameter vector in native space;
%
% (2) hhmm(r, ptrans, 'trans')
% 
%     where r is the structure generated by fitModel, ptrans is the parameter vector in
%     transformed space, and 'trans' is a flag indicating this.
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2013 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Transform paramaters back to their native space if needed
if ~isempty(varargin) && strcmp(varargin{1},'trans');
    [p pstruct] = hhmm_transp(r, p);
end

% Fixed configuration elements
%
% Number of possible outcomes
m = r.c_prc.n_outcomes;

% Get model tree (N is for node)
N = pstruct.N;

% Check constraints
if ~isempty(N{1}.V)
    error('tapas:hgf:hhmm:IllegEntryProbRoot', 'Illegal entry probability for root node.');
end

for id = 1:length(N)
    if isempty(N{id}.A) == isempty(N{id}.B)
        error('tapas:hgf:hhmm:IllegCombOfAB', 'Illegal combination of A and B for node no. %d.', id);
    end
    
    if length(N{id}.children(:)) ~= size(N{id}.A,2)
        error('tapas:hgf:hhmm:NumOfChildIncons', 'Number of children inconsistent with A for node no. %d.', id);
    end
    
    if ~isempty(N{id}.A) && any(sum(N{id}.A,2)>1)
        error('tapas:hgf:hhmm:IllegA', 'Illegal transition matrix A for node no. %d: row sums have to be less than or equal to 1.', id);
    end
    
    if ~isempty(N{id}.A)
        for cid = N{id}.children
            cidx = find(N{id}.children==cid);
            if ~isempty(N{cid}.children) && N{id}.A(cidx,cidx) ~= 0
                error('tapas:hgf:hhmm:IllegASelf', 'Illegal transition matrix A for node no. %d: only production nodes may have self-transitions.', id);
            end
        end
    end
    
    if ~isempty(N{id}.B) && sum(N{id}.B(:))~=1
        error('tapas:hgf:hhmm:IllegB', 'Illegal outcome contingency vector B for node no. %d.', id);
    end
    
    if ~isempty(N{id}.children)
        Vsum = 0;
        for cid = N{id}.children
            Vsum = Vsum + N{cid}.V;
        end
        
        if Vsum ~= 1
            error('tapas:hgf:hhmm:IllegV', 'Illegal vertical transition probabilities V from node no. %d.', id);
        end
    end
end

% Flatten the tree into one large transition matrix
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Find production nodes
pn = [];
for id = 1:length(N)
    if isempty(N{id}.children)
        pn = [pn, id];
    end
end
flatDim = length(pn);

% Initialize outcome contingency matrix
Bflat = NaN(m,flatDim);

% Fill Bflat
for i = 1:flatDim
    Bflat(:,i) = N{pn(i)}.B';
end

% Check Bflat
if any(~isfinite(Bflat(:)))
    error('tapas:hgf:hhmm:NoOutConMat', 'Could not construct outcome contingency matrix.');
end

% Initialize flattened transitioned matrix
Aflat = NaN(flatDim);

% Fill Aflat
for i = 1:flatDim
    for j = 1:flatDim
        if N{pn(i)}.parent == N{pn(j)}.parent
            % If the production nodes are siblings, read out
            % their parent's A matrix
            pid = N{pn(i)}.parent;
            idx = find(N{pid}.children==pn(i));
            jdx = find(N{pid}.children==pn(j));
            Aflat(i,j) = N{pid}.A(idx,jdx);
        else
            % Otherwise, determine their lowest common ancestor
            % and use that to calculate the transition probability
            caid = ca(N,pn(i),pn(j));
            tp = 1;

            nid  = pn(i);
            pid  = N{nid}.parent;
            nidx = find(N{pid}.children==nid);

            % Move up to one node below lowest common ancestor from start node pn(i)
            while pid ~= caid
                aend = 1-sum(N{pid}.A(nidx,:));
                tp = tp*aend;

                nid = pid;
                pid = N{nid}.parent;
                nidx = find(N{pid}.children==nid);
            end

            % Do the horizontal transition to the ancestral line of target node pn(j)
            ancj = anc(N,pn(j));
            while ancj(1) ~= caid
                ancj(1) = [];
            end
            ancj(1) = [];

            caidxi = nidx;
            caidxj = find(N{caid}.children==ancj(1));
            tp = tp*N{caid}.A(caidxi,caidxj);
            
            % Go down to the target node pn(j)
            ancj(1) = [];
            while ~isempty(ancj)
                tp = tp*N{ancj(1)}.V;
                ancj(1) = [];
            end

            Aflat(i,j) = tp;
        end
    end
end

% Check Aflat
if any(~isfinite(Aflat(:)))
    error('tapas:hgf:hhmm:NoFlatTransMat', 'Could not flatten transition matrix.');
end

% Calculate prior probabilities of production nodes
pnp = NaN(1,flatDim);
for i = 1:flatDim
    anci = anc(N,pn(i));
    anci(1) = [];
    p = 1;
    while ~isempty(anci)
        p = p*N{anci(1)}.V;
        anci(1) = [];
    end
    pnp(i) = p;
end

% Check pnp
if sum(pnp) ~= 1
    error('tapas:hgf:hhmm:NoPriorProdNodes', 'Cannot calculate prior probabilities of production nodes.');
end

% Input and number of trials
u = r.u(:,1);
n = length(u);

% Initialize alpha-prime
alpr = NaN(n,flatDim);

% alpr(1,:)
altmp = pnp.*Bflat(u(1),:);
llh = sum(altmp);
alpr(1,:) = altmp./llh;

% Pass through alpha-prime update loop
for k = 2:1:n
    if not(ismember(k, r.ign))
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Effect of input u(k)
        %%%%%%%%%%%%%%%%%%%%%%
        
        altmp = Bflat(u(k),:).*(alpr(k-1,:)*Aflat);
        llh = sum(altmp);
        alpr(k,:) = altmp./llh;
    else
        alpr(k,:) = alpr(k-1,:);
    end
end

% Predicted states
alprhat = [pnp; alpr];
alprhat(end,:) = [];

% Create result data structure
traj = struct;

traj.alpr    = alpr;
traj.alprhat = alprhat;

% Create matrix needed by observation model
infStates = traj.alpr;

end % function hhmm

% ----------------------------------------------------------------------------------------
% Find lowest common ancestor of nodes
function ca = ca(N,ida,idb)
    % Find ancestors of ida and idb
    anca = anc(N,ida);
    ancb = anc(N,idb);
    
    % Determine lowest common ancestor
    ca = NaN;
    aa = anca(1);
    ab = ancb(1);
    while aa == ab && ~isempty(anca) && ~isempty(ancb)
        ca = aa;
        anca(1) = [];
        ancb(1) = [];
        aa = anca(1);
        ab = ancb(1);
    end
    
    if aa == ab
        ca = aa;
    end
end

% ----------------------------------------------------------------------------------------
% Find ancestors of a node
function anc = anc(N,id)
    anc = id;
    idt = N{id}.parent;
    while ~isempty(idt)
        anc = [idt, anc];
        idt = N{idt}.parent;
    end
end
