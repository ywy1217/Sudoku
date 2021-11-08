function candidate_set = Sudoku_solver_v0(obj)
    obj.solution = obj.puzzle;
    % flag for solved numbers
    fix = obj.solution ~= 0;
%     figure;image(fix*255);grid off;axis image;
    % prep initial feasible set
%     fprintf('prep initial candidate set..\n');
%     tic;
    candidate_set = struct();
    count = 1;
    for i_col = 1:9
        for i_row = 1:9
            candidate_set(count).i_row = i_row;
            candidate_set(count).i_col = i_col;
            if fix(i_row, i_col)
                candidate_set(count).numbers = obj.solution(i_row, i_col);
                candidate_set(count).cardinality = 1;
            else
                % update feasible_set based on fixed grid in the puzzle
                candidate_set(count).numbers = 1:9;
                row_elements = obj.solution(i_row, :);
                col_elements = obj.solution(:, i_col);
                blk_elements = obj.solution(sub2block_sub(i_row), sub2block_sub(i_col));
                fix_elements = unique(nonzeros(cat(1, row_elements(:),col_elements(:), blk_elements(:) )));
                candidate_set(count).numbers( fix_elements ) = [];
                candidate_set(count).cardinality = length(candidate_set(count).numbers);
            end
            count = count + 1;
        end
    end
    candidate_set = reshape(candidate_set, 9, 9); 

%     fprintf('Done! (%.3fs)\n', toc);
    
    % prep indices for all 27 regions
    region_indices_all = extractIdx_27regions;
    
    % stack for backtracking
    candidate_stack = struct();
    stack_depth = 0;
    max_depth = 0;
    
    [sorted_candidate_set, ~] = sort_candiate_set(candidate_set, fix);
    n_unresolved_grid = length(sorted_candidate_set);
    shortest_len = 1;
    empty_flag = false;
    while( n_unresolved_grid > 0)
        % update candidate set until guessing is needed
%         fprintf('Update candidate_set until all singleton sets are addressed...\n');
%         tic;
        singleton_count = sum([sorted_candidate_set.cardinality]==1, 'all');
        while(singleton_count>0 && ~empty_flag)
            for idx = 1 : singleton_count
                i_row = sorted_candidate_set(idx).i_row;
                i_col = sorted_candidate_set(idx).i_col;
                num = sorted_candidate_set(idx).numbers;

                fix(i_row, i_col) = true;
                obj.solution(i_row, i_col) = num;
                candidate_set = update_candidate_set(candidate_set, num, i_row, i_col);
                if ~all([candidate_set.cardinality],'all')
                    empty_flag = true;
                    break;
                end
            end
            % update sorted set
            [sorted_candidate_set, ~] = sort_candiate_set(candidate_set, fix);
            singleton_count = sum([sorted_candidate_set.cardinality]==1, 'all');
        end
%         fprintf('Done! (%.3fs)\n', toc);
        n_unresolved_grid = length(sorted_candidate_set);
        
        %check intermitent results
%         obj.print_sudoku(true);

%         fprintf('Excluding candidates breaking the rules...\n');
%         tic;
        if n_unresolved_grid > 0
            shortest_len = sorted_candidate_set(1).cardinality;
            assert(shortest_len ~= 1);
        end
        
        for idx = 1:27
            region_indices_all{idx} = region_indices_all{idx} ( ~fix( region_indices_all{idx} ) );
            preserve_order = false;
            while(~preserve_order && shortest_len ~= 1)
                
                sorted_indices_per_region = sort_region(candidate_set, region_indices_all{idx});
                n_unfixed_grid_per_region = length(sorted_indices_per_region);
                
                region_indices_all{idx} = sorted_indices_per_region;
                preserve_order = true;
                region_set = candidate_set( sorted_indices_per_region );
                
                % remove unfeasible numbers from candidate set
                for k = 1:n_unfixed_grid_per_region-1
                    unique_set = unique([region_set(1:k).numbers]);
                    if length(unique_set) <= k
                        
                        for m = k+1:n_unfixed_grid_per_region
                            p = sorted_indices_per_region(m);
                            numbers_grid = candidate_set( p ).numbers;
                            [candidate_set( p ).numbers, i_rem] = setdiff(numbers_grid, unique_set);
                            new_len = length(i_rem);
                            if length(i_rem) < length(numbers_grid)
                                candidate_set(p).cardinality = new_len;
                                preserve_order = false;
                            end
                            if new_len < shortest_len
                                shortest_len = new_len;
                            end
                        end

                        if ~preserve_order
                            break;
                        end
                    end
                end
            end
            
            if shortest_len == 1
                break;
            end
        end
        
        %check intermitent results
%         obj.print_sudoku(true);

        % This section is much slower given that it goes through all the numbers
%         for idx = 1:n_unresolved_grid
%             if shortest_len == 1
%                 break
%             end
%             
%             i_row = sorted_candidate_set(idx).i_row;
%             i_col = sorted_candidate_set(idx).i_col;
%             for guess_num = [sorted_candidate_set(idx).numbers]
%                 if ~check_feasibility(candidate_set, guess_num, i_row, i_col)
%                     fprintf("Found unfeasible candidate (%d, %d, %d)!\n",i_row, i_col, guess_num);
%                     candidate_set(i_row, i_col).numbers = setdiff(candidate_set(i_row, i_col).numbers, guess_num);
%                     new_len = length(candidate_set(i_row, i_col).numbers);
%                     candidate_set(i_row, i_col).cardinality = new_len;
%                     if new_len == 1
%                         shortest_len = 1;
%                         break;
%                     end    
%                 end
%             end
%             
%         end
%         fprintf('Done! (%.3fs)\n', toc);
        
        if n_unresolved_grid > 0 && ~empty_flag
            [sorted_candidate_set, ~] = sort_candiate_set(candidate_set, fix);
            if shortest_len == 1
%                 fprintf('Found new singleton sets... \n');
            else
                if stack_depth == 0
                    fprintf('Guessing is needed... \n');
                end
                % recursive backtracking
                i_row = sorted_candidate_set(1).i_row;
                i_col = sorted_candidate_set(1).i_col;
                num = sorted_candidate_set(1).numbers(1);
                assert(length(num) == 1);
                
                % push
                stack_depth = stack_depth + 1;
                max_depth = max(max_depth, stack_depth);
                candidate_stack(stack_depth).i_row = i_row;
                candidate_stack(stack_depth).i_col = i_col;
                candidate_stack(stack_depth).guess_num = num;
                candidate_stack(stack_depth).candidate_set = candidate_set;
                candidate_stack(stack_depth).fix = fix;
                candidate_stack(stack_depth).solution = obj.solution;
                
%                 fix(i_row, i_col) = true;
%                 obj.solution(i_row, i_col) = num;
                candidate_set(i_row, i_col).numbers = num;
                candidate_set(i_row, i_col).cardinality = 1;
                shortest_len = 1;
%                 candidate_set = update_candidate_set(candidate_set, num, i_row, i_col);
                [sorted_candidate_set, ~] = sort_candiate_set(candidate_set, fix);
            end            
        else
            if all(obj.solution, 'all') 
                fprintf('Found a solution!\n');
            else
                if stack_depth == 0
                    fprintf('This puzzle doesn''t have a solution!\n');
                else
                    % pop
                    obj.solution = candidate_stack(stack_depth).solution;
                    fix = candidate_stack(stack_depth).fix;
                    candidate_set = candidate_stack(stack_depth).candidate_set;
                    i_row = candidate_stack(stack_depth).i_row;
                    i_col = candidate_stack(stack_depth).i_col;
                    guess_num = candidate_stack(stack_depth).guess_num;
                    stack_depth = stack_depth - 1;
                    
                    candidate_set(i_row, i_col).numbers = setdiff( candidate_set(i_row, i_col).numbers, guess_num );
                    candidate_set(i_row, i_col).cardinality = length(candidate_set(i_row, i_col).numbers);
                    assert(candidate_set(i_row, i_col).cardinality > 0);
                    [sorted_candidate_set, ~] = sort_candiate_set(candidate_set, fix);
                    
                    empty_flag = false;
                end
            end
        end

    end
    
    % Check solution complies with Sudoku rules
    check_rules(obj.solution);
end

function sorted_indices = sort_region(set, region_indices)
% heuristic sorting for removing unfeasible candidates easily
if isempty(region_indices)
    sorted_indices = region_indices;
    return;
end
region_set = set(region_indices);
n_unresolved_grids_region = length(region_indices);
% sort by likelihood (average occurence) - descending
occur = zeros(9,1);
for k = [region_set.numbers]
    occur(k) = occur(k) + 1;
end
likelihood = zeros(n_unresolved_grids_region,1);
for k = 1:n_unresolved_grids_region
    likelihood(k) = mean( occur(region_set(k).numbers), 'all');
end
[~, sort_likelihood] = sort(likelihood, 'descend');
% region_indices = region_indices(sort_likelihood);
% region_set = region_set(sort_likelihood);

sorted_indices = region_indices(sort_likelihood);
% sort by cardinality - ascending (preserves likelihood order)
% [~, sort_card] = sort([region_set.cardinality]);
% sorted_indices = region_indices(sort_card);

end

function flag = check_feasibility(candidate_set, guess_num, i_row, i_col)
% check whether a number in the candidate set is feasible
block_r = ceil(i_row/3);
block_c = ceil(i_col/3);
flag = true;
for idx = 1:9
    % check all rows
    if idx ~= i_row
        flag = ismember(guess_num, [ candidate_set(idx, (1:9)~=i_col).numbers ]);
    else
        flag = length(setdiff(unique([ candidate_set(idx, (1:9)~=i_col).numbers ]), guess_num)) == 8;
    end
    if flag == false
        return;
    end
    % check all cols
    if idx ~= i_col
        flag = ismember(guess_num, [ candidate_set((1:9)~=i_row, idx).numbers ]);
    else
        flag = length(setdiff(unique([ candidate_set((1:9)~=i_row, idx).numbers ]), guess_num)) == 8;
    end
    if flag == false
        return;
    end
    
end
% check all affected blocks 
if flag
    block_r_range = (block_r-1)*3+(1:3);
    block_r_range_diff = block_r_range(block_r_range ~= i_row);
    block_c_range = (block_c-1)*3+(1:3);
    block_c_range_diff = block_c_range(block_c_range ~= i_col);
    
    flag = length(setdiff(unique([ candidate_set( block_r_range_diff, block_c_range ).numbers, ...
        candidate_set( block_r_range, block_c_range_diff ).numbers]), guess_num)) == 8;
    if flag
        for idx = 1:3
            % same block row
            if idx ~= block_c
                flag = ismember(guess_num, [ candidate_set( block_r_range_diff, (idx-1)*3+(1:3) ).numbers] );
            end
            if flag == false
                return;
            end
            % same block col
            if idx ~= block_r
                flag = ismember(guess_num, [ candidate_set( (idx-1)*3+(1:3), block_c_range_diff ).numbers] );
            end
            if flag == false
                return;
            end
        end
    end
end


end

function check_rules(mat)
    all_indices = extractIdx_27regions;
    % check uniqueness
    for idx = 1:length(all_indices)
        v = mat(all_indices{idx});
        v = v(:);
        v = v(v~=0);
        assert( length(v) == length(unique(v)) );
    end
    % all cells are resolved
    assert(all(mat,'all'));
end

function new_set = update_candidate_set(set, num, i_row, i_col)
    
    % remove num from i_row
    for idx = 1:9
        if idx ~= i_col
            remove_idx = set( i_row+(idx-1)*9 ).numbers == num;
            if any(remove_idx)
                set( i_row+(idx-1)*9 ).numbers = set( i_row+(idx-1)*9 ).numbers(~remove_idx);
                set( i_row+(idx-1)*9 ).cardinality = set( i_row+(idx-1)*9 ).cardinality - 1;
            end
        end
    end
    % remove num from i_col
    for idx = 1:9
        if idx ~= i_row
            remove_idx = set( idx, i_col ).numbers == num;
            if any(remove_idx)
                set( idx, i_col).numbers = set( idx, i_col ).numbers(~remove_idx);
                set( idx+(i_col-1)*9 ).cardinality = set( idx+(i_col-1)*9 ).cardinality - 1;
            end
        end
    end
    % remove num from the block containing (i_row, i_col)
    for idx_1 = sub2block_sub(i_row)
        for idx_2 = sub2block_sub(i_col)
            if idx_1 ~= i_row && idx_2 ~= i_col
                remove_idx = set(idx_1, idx_2).numbers == num;
                if any(remove_idx)
                    set(idx_1, idx_2).numbers = set(idx_1, idx_2).numbers(~remove_idx);
                    set(idx_1, idx_2).cardinality = set( idx_1, idx_2 ).cardinality - 1;
                end
            end
        end
    end
    
    new_set = set;
end

function [sorted_set, original_idx] = sort_candiate_set(set, fix)
% take true candidate set (excluding solved grid)
z = (1:numel(set))';
original_idx = z(~fix);
set = set(~fix);

% sort the candidate set based on cardinality
[~, idx_permute] = sort([set.cardinality]);
original_idx = original_idx( idx_permute );
sorted_set = set(idx_permute);
end

% function solution = apply_singleton(set)
%     singleton_idx = [set.cardinality] == 1;
%     solution = zeros(9,9);
%     solution(singleton_idx) = [set(singleton_idx).numbers];
% end

function block_sub = sub2block_sub(sub)
    block_sub = 3*(ceil(sub/3)-1) + (1:3); 
end

function all_indices = extractIdx_27regions

all_indices = cell(27, 1);
for idx = 1:9
    % all rows
    all_indices{idx} = idx + 9*(0:8)';
    % all cols
    all_indices{idx+9} = 9*(idx-1) + (1:9)';
    % all blocks
    [block_r, block_c] = ind2sub([3 3], idx);
    block_r_range = (block_r-1)*3+(1:3);
    block_c_range = (block_c-1)*3+(1:3);
    [block_r_range, block_c_range]=ndgrid(block_r_range, block_c_range);
    all_indices{idx+18} = sub2ind([9 9], block_r_range(:), block_c_range(:));
end

end
