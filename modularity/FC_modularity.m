function FC_modularity(nsubs,in_dir,outfile)
addpath(genpath("BCT/"))

seed=123;
rng(seed,'twister')
network='example';
if strcmp(network,'example')
    comb_start = [1 7 11];
    comb_end = [6 10 14];
    names={'RP','CC','EP'};
    ROI_count=size(comb_start,2);
end


fcount=nsubs;
int_ROI_count=ROI_count;
cols=2+ROI_count+(ROI_count*(ROI_count-1)/2);
results=table('Size',[fcount cols],'VariableTypes',['string', repmat({'double'},1,cols-1)]);
results.Properties.VariableNames{1}='subjid';
results.Properties.VariableNames{2}='modularity';
count = 3;
for i=1:size(names,2)
    results.Properties.VariableNames{count}=sprintf('Rec.%s',names{i});
    count=count+1;
end
for i=1:size(names,2)
    for j=(i+1):size(names,2)
        results.Properties.VariableNames{count}=sprintf('Int.%sx%s',names{i},names{j});
        count=count+1;
    end
end


modularity = zeros(fcount,1);
M_data_all = zeros(fcount,int_ROI_count,int_ROI_count);
M_dist_all = zeros(fcount,int_ROI_count,int_ROI_count);
A_mat=zeros(fcount,int_ROI_count,int_ROI_count);

for i=1:fcount
    filename = sprintf('%s/resultsROI_Subject%03d_Condition001.mat',in_dir,i);
    infile=load(filename);

    W=infile.Z;
    W(1:1+size(W,1):end) = 0;
    M_data=zeros(100,int_ROI_count,int_ROI_count);
    Q_final=0;
    for c=1:100                     % Loop 100 times to generate a stable estimate
        n  = size(W,1);             % number of nodes
        M  = 1:n;                   % initial community affiliations
        Q0 = -1; Q1 = 0;            % initialize modularity values
        while Q1-Q0>1e-5           % while modularity increases
            Q0 = Q1;                % perform community detection
            [M, Q1] = community_louvain(W, [], M, 'negative_asym');
        end
        if Q1>Q_final               %saves the optimal Q value across all 100 runs
            Q_final=Q1;
        end
        
        % From group assignment, generate a full matrix with temp(a,b)=1 if
        % a and b are assigned to the same group and 0 otherwise
        temp=zeros(size(M,1)); 
        for a=1:size(M,1)
            for b=1:size(M,1)
                if M(a)==M(b)
                    temp(a,b)=1;
                end
            end
        end
        
        % Collapse full matrix into an affinity matrix based on groups
        for m=1:ROI_count
            for n=1:ROI_count
                M_data(c,m,n)=sum(temp(comb_start(m):comb_end(m),comb_start(n):comb_end(n)),'all');
            end
        end

    end
    M_data_all(i,:,:)=squeeze(mean(M_data,1));
    
    
    % iterate 1000 times with randomized connectivity matrices to generate
    % a null distribution
    M_dist=zeros(1000,int_ROI_count,int_ROI_count);
    Q_dist=zeros(1000,1);
    for c=1:1000
        temp=tril(W,-1);
        W_lin=temp(temp~=0);
        W_lin=W_lin(randperm(length(W_lin)));                    

        W_shuffle = tril(ones(size(W,2)),-1);
        W_shuffle(W_shuffle>0) = W_lin;
        W_shuffle = W_shuffle + W_shuffle';


        n  = size(W_shuffle,1);             % number of nodes
        M  = 1:n;                   % initial community affiliations
        Q0 = -1; Q1 = 0;            % initialize modularity values
        while Q1-Q0>1e-5           % while modularity increases
            Q0 = Q1;                % perform community detection
            [M, Q1] = community_louvain(W_shuffle, [], M, 'negative_sym');
        end
        Q_dist(c)=Q1;
        temp=zeros(size(M,1));
        for a=1:size(M,1)
            for b=1:size(M,1)
                if M(a)==M(b)
                    temp(a,b)=1;
                end
            end
        end
        for m=1:ROI_count
            for n=1:ROI_count
                M_dist(c,m,n)=sum(temp(comb_start(m):comb_end(m),comb_start(n):comb_end(n)),'all');
            end
        end
    end

    M_dist_all(i,:,:)=squeeze(mean(M_dist,1));
    
    %Normalize individual estimates with the individual's null distribtuion
    modularity(i)=(Q_final-mean(Q_dist))/std(Q_dist);
    for m=1:ROI_count
        for n=1:ROI_count
            A_mat(i,m,n)=(M_data_all(i,m,n)-mean(M_dist(:,m,n),1))/std(M_dist(:,m,n));
        end
    end
end

% save data to file
for i=1:fcount
    Q_factor=modularity(i);
    Allegiance_matrix=squeeze(A_mat(i,:,:));
    %save(file_out, 'Q_factor','Allegiance_matrix');

    temp=tril(ones(size(Allegiance_matrix)),-1);

    results.subjid(i)=sprintf('subj_%03d',i);
    
    results{i,3:2+ROI_count}=diag(Allegiance_matrix)';
    results{i,3+ROI_count:end}=Allegiance_matrix(temp~=0)';
end
results.modularity=modularity;

writetable(results,outfile)

end
