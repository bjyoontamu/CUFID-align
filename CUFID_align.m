function net_align = CUFID_align(Net_id_list, input_folder, id_flag, out_file)
%-------------------------------------------------------------------------
% 1. CUFID_align finds the pairwise alignment (one-to-one mapping) of biological networks
%
%
% 2. Input and output arguments of the function
% Input:
%       Net_id_list  : List of the name for the PPI networks
%       input_folder : File path for the input data
%       id_flag      : If the nodes are in numeric format ("species id+number"),
%                      set id_flag = 1
%                      Otherwise, set id_flag = 0
%       out_file     : Name of the output file
% Output:
%       net_align    : Alignment results
%
% 3. How to run using the test data
% Please type the following commands in the matlab command window
%
%       Net_id_list = {'a', 'b'};
%       input_folder = './test';
%       id_flag = 1;
%       out_file = './test/test_out.txt';
%       net_align = CUFID_align(Net_id_list, input_folder, id_flag, out_file);
%
%
% 4. File format for the input and output
% Input file format
%       - Network file:
%               It supports tab-delimited file format.
%               If two proteins have a interaction,
%               a1  a2
%               a1  a3
%               a4  a2
%               Note that it also supports the edge weight:
%               a1  a2  0.7
%               a1  a3  0.5
%               a4  a2  0.3
%       - Similarity score file:
%               Elements in the first column indicate the nodes in the one network
%               Elements in the second column indicate the nodes in another network
%               Value in the third column indicates their node simialrity score.
%
%               As an example, a-b.sim may be as follows:
%               a1  b1  100
%               a1  b3  78
%               a2  b2  80
%               a2  b4  45
%               a3  b4  60
%
% Output file format
%       Each row in the output file is the aligned node pair.
%               a1  b1
%               a2  b2
%               a3  b4
%
%
%
% 5. Reference
% For more information on the algorithms, please see:
%       H. Jeong, X. Qian, and B.-J. Yoon (2016), Effective comparative
%       analysis of protein-protein interaction networks by measuring
%       the steady-state network flow using a Markov model,
%       BMC Bioinformatics, In Press.
%
% 6. Contact: bjyoon@ece.tamu.edu
%-------------------------------------------------------------------------



%--------------------------------------------------------------------------
% read networks
%--------------------------------------------------------------------------
[G S L edges nodes Sims] = ReadNetworks(input_folder, Net_id_list,id_flag);
%--------------------------------------------------------------------------




tic
%--------------------------------------------------------------------------
% construct integrated network
%--------------------------------------------------------------------------
A = G;  
for ind =1:length(L)    % make a stochastic matrix
    A{ind} = spdiags(1./sum(A{ind},2), 0, L(ind), L(ind)) * A{ind};
end

S12 = (spdiags(1./sum(S{1,2},2), 0, L(1), L(1)) * S{1,2});
S21 = ( S{1,2} * spdiags(1./sum(S{1,2})', 0, L(2), L(2)) );

T = [A{1} S12; S21' A{2}];
TT = spdiags(1./sum(T,2), 0, size(T,1), size(T,2)) * T;




%--------------------------------------------------------------------------
% Compute steady-state probability
%--------------------------------------------------------------------------
Pi = Steady(TT,repmat(1/size(TT,1),1,size(TT,1)));

PiX = Pi(1:L(1));
PiY = Pi(L(1) + 1 : L(1)+L(2));

F1 = spdiags(PiX, 0, L(1), L(1)) * S12;

F2 = bsxfun(@times, PiY', S21);
FlowM = F1 + F2;



init_ind = find(FlowM);
K=1;
alpha = 0.90;
for ii = 1:K
    FlowM = alpha*FlowM + (1-alpha)*A{1}*FlowM*A{2}';
end
FlowM_b = FlowM;
th = prctile(FlowM(find(FlowM)), 90);   % compute threshold
FlowM(find(FlowM<th)) = 0;
FlowM(init_ind) = FlowM_b(init_ind);
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
% Construct maximum-weighted bipartite matching
%--------------------------------------------------------------------------
[val m1 m2] = bipartite_matching(FlowM);


t2=toc;
sprintf('Elapsed time for Alignment: %f seconds',t2)



%--------------------------------------------------------------------------
% Output
%--------------------------------------------------------------------------
fod = fopen(out_file,'w');
for ind = 1:length(m1)

    if id_flag==0
        fprintf(fod,'%s %s \n', cell2mat(nodes{1}(m1(ind))), cell2mat(nodes{2}(m2(ind))));
        alignment{ind} = { cell2mat(nodes{1}(m1(ind))), cell2mat(nodes{2}(m2(ind))) };
    else
        fprintf(fod,'%s %s\n',strcat(Net_id_list{1},num2str(m1(ind))),strcat(Net_id_list{2},num2str(m2(ind))));
        alignment{ind} = { strcat(Net_id_list{1},num2str(m1(ind))), strcat(Net_id_list{2},num2str(m2(ind))) };
    end;
end
fclose(fod);
net_align = alignment;




end

function [G S L edges nodes Sims] = ReadNetworks(input_folder, Net_id_list,id_flag)


M=length(Net_id_list);
nets=cell(1,M);
nodes=cell(1,M);
edges=cell(1,M);
Sims=cell(M,M);
Net_id_list=lower(Net_id_list);

for i=1:M
    fid = fopen(strcat(input_folder,'/',upper(Net_id_list{i}),'.net'));
    if id_flag==0
        nets{i} = textscan(fid,'%s\t%s\t%s', 1000000);
    else
        nets{i} = textscan(fid,strcat(Net_id_list{i},'%d\t',Net_id_list{i},'%d\t%s'), 1000000);
    end;
    fclose(fid);
    if isempty(nets{i}{3}{1})
        nets{i}(3)=[];
    end;
end;

for i=1:M
    for j=i+1:M
        fid = fopen(strcat(input_folder,'/',upper(Net_id_list{i}),'-',upper(Net_id_list{j}),'.sim'));
        if id_flag==0
            Sims{i,j} = textscan(fid,'%s\t%s\t%f', 10000000);
        else
            Sims{i,j}= textscan(fid,strcat(Net_id_list{i},'%d\t',Net_id_list{j},'%d\t%f'), 1000000);
        end;
        fclose(fid);
    end;
end;

if (id_flag==0)
for i = 1: M
    nodes{i} = unique([nets{i}{1}; nets{i}{2}]);
    for j = i + 1: M
        nodes{i} = unique([nodes{i}; Sims{i, j}{1}]);
    end
    
    for j = 1: i - 1
        nodes{i} = unique([nodes{i}; Sims{j, i}{2}]);
    end
    ids = cell(1, M);
    for j = 1: M
        ids{j} = [];
    end
    
    for j = i + 1: M
        [~, tempList] = ismember(Sims{i, j}{1}, nodes{i});
        Sims{i, j}{1} = {tempList};
    end
    
    for j = 1: i - 1
        [~, tempList] = ismember(Sims{j, i}{2}, nodes{i});
        Sims{j, i}{2} = {tempList};
    end
end

for i = 1: M
    edges{i} = zeros(length(nets{i}{1}), 3);
    [~, edges{i}(:, 1)] = ismember(nets{i}{1}, nodes{i});
    [~, edges{i}(:, 2)] = ismember(nets{i}{2}, nodes{i});
    
    if length(nets{i}) == 2
        edges{i}(:, 3) = 1;
    else
        edges{i}(:, 3) = str2double(nets{i}{3});
    end
end

for i = 1: M
    for j = i + 1: M
        Sims{i, j}=[cell2mat(Sims{i, j}{1}), cell2mat(Sims{i, j}{2}), (Sims{i, j}{3})];
    end
end
else
    nnodes=cell(1,M);
    for i=1:M
        nnodes{i}=max([nets{i}{1};nets{i}{2}]);
        for j=i+1:M
            nnodes{i}=max([nnodes{i};Sims{i,j}{1}]);
        end;
        for j=1:i-1
            nnodes{i}=max([nnodes{i};Sims{j,i}{2}]);
        end;
        
    end;
    
    for i=1:M
        nodes{i}=1:nnodes{i};
    end;
    for i=1:M
        if length(nets{i})==2
            edges{i}=[nets{i}{1},nets{i}{2},ones(length(nets{i}{1}),1)];
        else
            edges{i}=[nets{i}{1},nets{i}{2},str2double(nets{i}{3})];
        end;
    end;
    
    
    for i=1:M
        for j=i+1:M
            Sims{i,j}=[double(Sims{i,j}{1}),double(Sims{i,j}{2}), double(Sims{i,j}{3})];
        end;
    end;
    
    
end;

%---initializtion
M=length(edges);
L = zeros(1,M);
for i=1:M
    L(i)=length(nodes{i});
end;
G=cell(1,M);
for i=1:M
    ee = double(edges{i});
    Q = sparse(ee(:, 1), ee(:, 2), ee(:, 3), L(i), L(i));
    indices  = sub2ind(size(Q), ee(:, 2), ee(:, 1));
    Q(indices) = ee(:, 3);
    G{i} = Q;
end;
S=cell(M,M);
for i=1:M
    for j=i+1:M
        if isempty(Sims{i,j})
            Sims{i,j} = [1 1 0];
        end;
        S{i,j}=sparse(double(Sims{i,j}(:,1)),double(Sims{i,j}(:,2)),double(Sims{i,j}(:,3)),L(i),L(j));
    end;
end;


end % ends of function


%----------------------------------------------------------
%Computing steady state distribution
%----------------------------------------------------------
function Pi=Steady(A,p)
err=1;
J=p';
cnt=0;
while ((err>.00001)&(cnt<100))
    J0=J;
    JJ=(J'*A)';
    J=JJ/norm(JJ,1);
    err=norm(J0-J,2);
    cnt=cnt+1;
    if (cnt>1000)
        aaa=1;
    end;
end;
Pi=J;

end



function [val m1 m2 mi]=bipartite_matching(varargin)
% BIPARTITE_MATCHING Solve a maximum weight bipartite matching problem
%
% [val m1 m2]=bipartite_matching(A) for a rectangular matrix A 
% [val m1 m2 mi]=bipartite_matching(x,ei,ej,n,m) for a matrix stored
% in triplet format.  This call also returns a matching indicator mi so
% that val = x'*mi.
%
% The maximum weight bipartite matching problem tries to pick out elements
% from A such that each row and column get only a single non-zero but the
% sum of all the chosen elements is as large as possible.
%
% This function is slightly atypical for a graph library, because it will
% be primarily used on rectangular inputs.  However, these rectangular
% inputs model bipartite graphs and we take advantage of that stucture in
% this code.  The underlying graph adjency matrix is 
%   G = spaugment(A,0); 
% where A is the rectangular input to the bipartite_matching function.
%
% Matlab already has the dmperm function that computes a maximum
% cardinality matching between the rows and the columns.  This function
% gives us the maximum weight matching instead.  For unweighted graphs, the
% two functions are equivalent.
%
% Note: If ei and ej contain duplicate edges, the results of this function
% are incorrect.
%
% See also DMPERM
%
% Example:
%   A = rand(10,8); % bipartite matching between random data
%   [val mi mj] = bipartite_matching(A);
%   val

% David F. Gleich and Ying Wang
% Copyright, Stanford University, 2008-2009
% Computational Approaches to Digital Stewardship

% 2008-04-24: Initial coding (copy from Ying Wang matching_sparse_mex.cpp)
% 2008-11-15: Added triplet input/output
% 2009-04-30: Modified for gaimc library
% 2009-05-15: Fixed error with empty inputs and triple added example.

[rp ci ai tripi n m] = bipartite_matching_setup(varargin{:});

if isempty(tripi)
    error(nargoutchk(0,3,nargout,'struct'));
else    
    error(nargoutchk(0,4,nargout,'struct'));
end


if ~isempty(tripi) && nargout>3
    [val m1 m2 mi] = bipartite_matching_primal_dual(rp, ci, ai, tripi, n, m);
else
    [val m1 m2] = bipartite_matching_primal_dual(rp, ci, ai, tripi, n, m);
end
end

function [rp ci ai tripi n m]= bipartite_matching_setup(A,ei,ej,n,m)
% convert the input

if nargin == 1
    if isstruct(A)
        [nzi nzj nzv]=csr_to_sparse(A.rp,A.ci,A.ai);
    else
        [nzi nzj nzv]=find(A); 
    end
    [n m]=size(A);
    triplet = 0;
elseif nargin >= 3 && nargin <= 5    
    nzi = ei;
    nzj = ej;
    nzv = A;
    if ~exist('n','var') || isempty(n), n = max(nzi); end
    if ~exist('m','var') || isempty(m), m = max(nzj); end
    triplet = 1;
else    
    error(nargchk(3,5,nargin,'struct'));
end
nedges = length(nzi);

rp = ones(n+1,1); % csr matrix with extra edges
ci = zeros(nedges+n,1);
ai = zeros(nedges+n,1);
if triplet, tripi = zeros(nedges+n,1); % triplet index
else tripi = [];
end

%
% 1. build csr representation with a set of extra edges from vertex i to
% vertex m+i
%
rp(1)=0;
for i=1:nedges
    rp(nzi(i)+1)=rp(nzi(i)+1)+1;
end
rp=cumsum(rp); 
for i=1:nedges
    if triplet, tripi(rp(nzi(i))+1)=i; end % triplet index
    ai(rp(nzi(i))+1)=nzv(i);
    ci(rp(nzi(i))+1)=nzj(i);
    rp(nzi(i))=rp(nzi(i))+1;
end
for i=1:n % add the extra edges
    if triplet, tripi(rp(i)+1)=-1; end % triplet index
    ai(rp(i)+1)=0;
    ci(rp(i)+1)=m+i;
    rp(i)=rp(i)+1;
end
% restore the row pointer array
for i=n:-1:1
    rp(i+1)=rp(i);
end
rp(1)=0;
rp=rp+1;

%
% 1a. check for duplicates in the data
%
colind = false(m+n,1);
for i=1:n
    for rpi=rp(i):rp(i+1)-1
        if colind(ci(rpi)), error('bipartite_matching:duplicateEdge',...
            'duplicate edge detected (%i,%i)',i,ci(rpi)); 
        end
        colind(ci(rpi))=1;
    end
    for rpi=rp(i):rp(i+1)-1, colind(ci(rpi))=0; end % reset indicator
end

end

function [val m1 m2 mi]=bipartite_matching_primal_dual(...
                            rp, ci, ai, tripi, n, m)
% BIPARTITE_MATCHING_PRIMAL_DUAL                         

alpha=zeros(n,1); % variables used for the primal-dual algorithm
beta=zeros(n+m,1);
queue=zeros(n,1);
t=zeros(n+m,1);
match1=zeros(n,1);
match2=zeros(n+m,1);
tmod = zeros(n+m,1);
ntmod=0;


% 
% initialize the primal and dual variables
%
for i=1:n
    for rpi=rp(i):rp(i+1)-1
        if ai(rpi) > alpha(i), alpha(i)=ai(rpi); end
    end
end
% dual variables (beta) are initialized to 0 already
% match1 and match2 are both 0, which indicates no matches
i=1;
while i<=n
    % repeat the problem for n stages
    
    % clear t(j)
    for j=1:ntmod, t(tmod(j))=0; end
    ntmod=0;
    

    % add i to the stack
    head=1; tail=1;
    queue(head)=i; % add i to the head of the queue
    while head <= tail && match1(i)==0
        k=queue(head);
        for rpi=rp(k):rp(k+1)-1
            j = ci(rpi);
            if ai(rpi) < alpha(k)+beta(j) - 1e-8, continue; end % skip if tight
            if t(j)==0,
                tail=tail+1; queue(tail)=match2(j);
                t(j)=k;
                ntmod=ntmod+1; tmod(ntmod)=j;
                if match2(j)<1,
                    while j>0, 
                        match2(j)=t(j);
                        k=t(j);
                        temp=match1(k);
                        match1(k)=j;
                        j=temp;
                    end
                    break; % we found an alternating path
                end
            end
        end
        head=head+1;
    end
    
    if match1(i) < 1, % still not matched, so update primal, dual and repeat
        theta=inf;
        for j=1:head-1
            t1=queue(j);
            for rpi=rp(t1):rp(t1+1)-1
                t2=ci(rpi);
                if t(t2) == 0 && alpha(t1) + beta(t2) - ai(rpi) < theta,
                    theta = alpha(t1) + beta(t2) - ai(rpi);
                end
            end
        end
        
        for j=1:head-1, alpha(queue(j)) = alpha(queue(j)) - theta; end
        
        for j=1:ntmod, beta(tmod(j)) = beta(tmod(j)) + theta; end
            
        continue;
    end
        
    i=i+1; % increment i
end

val=0;
for i=1:n
    for rpi=rp(i):rp(i+1)-1
        if ci(rpi)==match1(i), val=val+ai(rpi); end
    end
end
noute = 0; % count number of output edges
for i=1:n
    if match1(i)<=m, noute=noute+1; end
end
m1=zeros(noute,1); m2=m1; % copy over the 0 array
noute=1;
for i=1:n
    if match1(i)<=m, m1(noute)=i; m2(noute)=match1(i);noute=noute+1; end
end

if nargout>3
    mi= false(length(tripi)-n,1);
    for i=1:n
        for rpi=rp(i):rp(i+1)-1
            if match1(i)<=m && ci(rpi)==match1(i), mi(tripi(rpi))=1; end
        end
    end
end



end
