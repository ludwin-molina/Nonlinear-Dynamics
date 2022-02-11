function [s,ix,iy,D]=dtw_mvx(A,B,d)
%INPUTS
%A is the time series 1
%B is the time series 2
%d is the selected distance. It could be 'squared','absolute' or
%'euclidean'
%Outputs
%s is the calculated distance between the time series
%ix is the vector order of the time series 1
%iy is the vector order of the time series 2

A=normalize(A);
B=normalize(B);
N=size(A,1);
M=size(B,1);
DTW=zeros(N,M);

for i=1:N
    for j=1:M
        DTW(i,j)= Inf;
    end
end

switch d
    case 'squared'
        for i=1:N
            for j=1:M
                cost=0;
                for k=1:size(A,2)
                    cost= cost+(A(i,k)-B(j,k))^2;
                end
                if and(i==1,j==1)
                    DTW(i,j)=cost;
                elseif i==1
                    DTW(i,j)=cost+DTW(i,j-1);
                elseif j==1
                    DTW(i,j)=cost+DTW(i-1,j);
                else
                    DTW(i,j)=cost+min([DTW(i-1,j),DTW(i,j-1),DTW(i-1,j-1)]);
                end
            end
        end
    case 'absolute'
        for i=1:N
            for j=1:M
                cost=0;
                for k=1:size(A,2)
                    cost= cost+abs(A(i,k)-B(j,k));
                end
                if and(i==1,j==1)
                    DTW(i,j)=cost;
                elseif i==1
                    DTW(i,j)=cost+DTW(i,j-1);
                elseif j==1
                    DTW(i,j)=cost+DTW(i-1,j);
                else
                    DTW(i,j)=cost+min([DTW(i-1,j),DTW(i,j-1),DTW(i-1,j-1)]);
                end
            end
        end
        
    case 'euclidean'
        for i=1:N
            for j=1:M
                cost=0;
                for k=1:size(A,2)
                    cost= cost+(A(i,k)-B(j,k))^2;
                end
                cost=sqrt(cost);
                if and(i==1,j==1)
                    DTW(i,j)=cost;
                elseif i==1
                    DTW(i,j)=cost+DTW(i,j-1);
                elseif j==1
                    DTW(i,j)=cost+DTW(i-1,j);
                else
                    DTW(i,j)=cost+min([DTW(i-1,j),DTW(i,j-1),DTW(i-1,j-1)]);
                end
            end
        end
        
    otherwise
        disp('please define a measurement')
        
end


ix(1)=N;
iy(1)=M;
l1=N;
l2=M;
n=1;
while or(l1~=1,l2~=1)
   n=n+1;
   if l1==1
       d1=DTW(1,l2);
       d2=DTW(l1,l2-1);
       d3=DTW(1,l2-1);
       ix(n)=l1;
       iy(n)=l2-1;
   elseif l2==1
       d1=DTW(l1-1,l2);
       d2=DTW(l1,l2);
       d3=DTW(l1-1,l2);
       ix(n)=l1-1;
       iy(n)=l2;
   else
       d1=DTW(l1-1,l2);
       d2=DTW(l1,l2-1);
       d3=DTW(l1-1,l2-1);
       if and(d3<=d1,d3<=d2)
        ix(n)=l1-1;
        iy(n)=l2-1;
        
       elseif and(d2<=d1,d2<=d3)
        ix(n)=l1;
        iy(n)=l2-1;
        
       elseif and(d1<=d2,d1<=d3)
        ix(n)=l1-1;
        iy(n)=l2;
        
       end
   end
    l1=ix(n);
    l2=iy(n);
end

ix=fliplr(ix)';
iy=fliplr(iy)';

s=DTW(N,M);

D=zeros(1,length(ix));

for i=1:length(ix)
    dx=0;
    for k=1:size(A,2)
    dx=dx+(A(ix(i),k)-B(iy(i),k))^2;
    end
    D(i)=sqrt(dx);
end

ix=ix';
iy=iy';
end

