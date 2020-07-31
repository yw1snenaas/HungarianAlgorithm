function [ cost, x ] = HungarianAlgorithm( C )
% The implementation of the Hungarian algorithm.
% Reference:
% https://csclab.murraystate.edu/~bob.pilgrim/445/munkres_old.html
% C: cost matrix
[n,m] = size(C);
k = min(n,m);
starMatrix = zeros(n,m);
primedMatrix = zeros(n,m);
coveredRow = zeros(n,1);
coveredColumn = zeros(1,m);
% Step1: For each row of the matrix, find the smallest element and subtract
% it from every element in its row. Go to Step 2.
C1 = zeros(n,m);
parfor i = 1:n
    tempV = min(C(i,:));
    C1(i,:) = C(i,:) - tempV;
end
% Step2: Find a zero (Z) in the resulting matrix. If there is no starred
% zero in its row or column, star (Z). Repeat for each element in the
% matrix. Go to Step 3.
for i = 1:n
    for j = 1:m
        if C1(i,j) == 0 && (sum(starMatrix(i,:))==0) && (sum(starMatrix(:,j))==0)
            starMatrix(i,j) = 1;
        end
    end
end
flag = zeros(1,5);
flag(1) = 1;
while flag(5) == 0
    % Step3: Cover each column containing a starred zero. If K columns are
    % covered the starred zeros describe a complete set of unique
    % assignments. In this case, Go to DOWN, otherwise, Go to Step 4.
    if flag(1) == 1
        parfor j = 1:m
            temp = sum(starMatrix(:,j));
            if temp > 0
                coveredColumn(j) = 1;
            end
        end
        if sum(coveredColumn) == k
            flag = zeros(1,5);
            flag(5) = 1;
        else
            flag = zeros(1,5);
            flag(2) = 1;
        end
    end
    % Step4: Find a noncovered zero and prime it. If there is no starred
    % zero in the row containing this primed zero. Go to Step 5. Otherwise,
    % cover this row and uncover the column containing the starred zero.
    % Continue in this manner until there are no uncovered zeros left. Save
    % the smallest uncovered value and Go to Step 6.
    if flag(2) == 1
        tempC = C1 + ones(n,1) * coveredColumn + coveredRow * ones(1,m);
        [idx1, idx2] = find(tempC == 0);
        for j = 1:length(idx1)
            primeIdx1 = idx1(j);
            primeIdx2 = idx2(j);
            primedMatrix(primeIdx1, primeIdx2) = 1;
            if sum(starMatrix(primeIdx1,:)) == 0
                % Go to Step 5:
                flag = zeros(1,5);
                flag(3) = 1;
                break;
            else
                coveredRow(primeIdx1) = 1;
                idx = find(starMatrix(primeIdx1,:) == 1);
                coveredColumn(idx) = 0;
            end
        end
        if flag(3) == 0
            m1 = coveredRow * ones(1,m) * inf;
            m1(find(isnan(m1)==1)) = 0;
            m2 = ones(n,1) * coveredColumn * inf;
            m2(find(isnan(m2)==1)) = 0;
            tempC = C1 + m1 + m2;
            minValue = min(min(tempC));
            [smallestIdx1, smallestIdx2] = find(tempC == minValue);
            smallestIdx1 = smallestIdx1(1);
            smallestIdx2 = smallestIdx2(1);
            % Go to Step 6:
            flag = zeros(1,5);
            flag(4) = 1;
        end
    end
    if flag(3) == 1
        % Step5: Construct a series of alternating primed and starred zeros
        % as follows. Let Z0 represent the uncovered primed zero found in
        % Step 4. Let Z1 denote the starred zero in the column of Z0 (if
        % any). Let Z2 denote the primed zero in the row of Z1 (there will
        % always be one). Continue until the series terminates at a primed
        % zero that has no starred zero in its column. Unstar each starred
        % zero of the series, star each primed zero of the series, erase
        % all primes and uncover every line in the matrix. Return to Step
        % 3.
        tempFlag = true;
        Z = [primeIdx1,primeIdx2];
        while tempFlag
            starIdx1 = find(starMatrix(:,primeIdx2) == 1);
            starMatrix(primeIdx1, primeIdx2) = 1;
            if isempty(starIdx1)
                tempFlag = false;
            else
                tempZ = [starIdx1(1),primeIdx2];
                Z = [Z;tempZ];
                starMatrix(tempZ(1),tempZ(2)) = 0;
                primeIdx1 = tempZ(1);
                primeIdx2 = find(primedMatrix(primeIdx1,:) == 1);
                primeIdx2 = primeIdx2(1);
                tempZ = [primeIdx1, primeIdx2];
                Z = [Z;tempZ];
            end
        end
        primedMatrix = zeros(n,m);
        coveredRow = zeros(n,1);
        coveredColumn = zeros(1,m);
        % Return to Step 3:
        flag = zeros(1,5);
        flag(1) = 1;
    end
    if flag(4) == 1
        % Step6: Add the value found in Step 4 to every element of each
        % covered row, and subtract it from every element of each uncovered
        % column. Return to Step 4 without altering any stars, primes, or
        % covered lines.
        C1 = C1 + coveredRow * ones(1,m) * minValue;
        uncoveredColumn = ones(1,m) - coveredColumn;
        C1 = C1 - ones(n,1) * uncoveredColumn * minValue;
        % Return to Step 4:
        flag = zeros(1,5);
        flag(2) = 1;
    end
end
% DONE: Assignment pairs are indicated by the positions of the starred
% zeros in the cost matrix. If C(i,j) is a starred zeros, then the element
% associated with row i is assigned to the element associated with column
% j.
x = starMatrix;
cost = 0;
[idx1,idx2] = find(starMatrix == 1);
for j = 1:length(idx1)
    cost = cost + C(idx1(j),idx2(j));
end


end

