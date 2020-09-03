% clickClust - matlab wrapper for click.exe
% Written by Adi Maron-Katz Jan 2016
% clusters numeric data based on R. Sharan, A. Maron-Katz and R. Shamir, "CLICK and EXPANDER: A system for clustering and visualizing gene expression data," Bioinformatics, vol. 19, no. 14, 2003. 
% input: data - a numeric matrix with n rows and m columns (rows will be clustered)
% Output: clust - a numeric vector of length m that indicates the cluster assignment for each row in data (0 indicates no assignment)
%         centers - a matrix of size cXm where c is the number of clusters. Each row is the centroid of one cluster
%		  stds - a matrix of size cXm where c is the number of clusters. Each row holds the standard deviation from the mean (centroid) in one cluster
function [ clust,centers,stds ] = clickClust( data )
rowIndices=(1:size(data,1))';
inputFileName='data4click.orig';
firstLine=size(data);
dlmwrite(inputFileName,firstLine, ' ');
dlmwrite(inputFileName,[rowIndices,data],'-append','delimiter','\t');
[status,cmdout] = dos('click.exe clickParams.txt','-echo');

sol=dlmread('./clickRes.res.sol');
clust=sol(:,2);
clusts=unique(clust);
clusts(find(clusts==0))=[];
centers=zeros(length(clusts),size(data,2));
stds=zeros(length(clusts),size(data,2));
for c=1:length(clusts)
    pos=find(clust==clusts(c));
    data_c=data(pos,:);
    centers(c,:)=mean(data_c);
    stds(c,:)=std(data_c);
end
delete('data4click.*');
delete('clickRes.anl');
delete('clickRes.res.sol');
end

