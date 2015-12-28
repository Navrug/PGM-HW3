function [ labels, centroids, distortion ] = k_means( data, k )
%K-MEANS Summary of this function goes here
%   Detailed explanation goes here
    n = size(data,1);
    centroids = data(randi(n,1,k),:);
    labels = zeros(1,n);
    old_labels = ones(1,n);
%     count = 0
    while not (isequal(labels, old_labels))
%       count = count +1
       distortion = 0;
       old_labels = labels;
       labels = ones(1,n);
       % Computing labels
       for i = 1:n
          best = norm(data(i,:) - centroids(1,:));
          for j = 2:k
              diff = norm(data(i,:) - centroids(j,:));
              if diff < best
                 best = diff;
                 labels(i) = j;
              end
          end
          distortion = distortion + best;
       end
       % Computing centers
       for j = 1:k
          centroids(k,:) = mean(data(labels==k, :));
       end
    end
end

