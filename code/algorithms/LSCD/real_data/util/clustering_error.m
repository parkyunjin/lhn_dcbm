function [error]=clustering_error(id,label,k)
n=length(label);
labelset=unique(label);
idset=unique(id);
error=zeros(k,k);
for i=1:k  %%id
    for j=1:k  %%label
        for node=1:n
            if (id(node)==i & label(node)~=labelset(j))
                error(i,j)=error(i,j)+1;
            end
        end
    end
end
