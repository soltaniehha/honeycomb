function []=display_contour(b)
c=zeros(cieling(size(b,2)/10));
for i=0:size(c)
    for j=0:size(c)
        for k=1:10
            for l=1:10
                c(i+1,j+1)=c(i+1,j+1)+b(i*10+k,j*10+l);
            end
        end
    end
end
plot(c);