function basins = p_myfunc_estimateBasin(average_vecs, cell_centers)
    [m, n] = size(average_vecs);
    basins = zeros(m, n); % 各セルが属するbasinを示す行列
    visited = false(m, n); % 走査済みセルを示すマトリックス
    
    % ベクトルの長さを計算
    vector_lengths = cellfun(@(vec) ifelse(isempty(vec), inf, norm(vec)), average_vecs);
    
    % ベクトルの長さが小さい順にソート
    [sorted_lengths, indices] = sort(vector_lengths(:));
    [sorted_rows, sorted_cols] = ind2sub([m, n], indices); % インデックスを行列の座標に変換
    
    % 離散化された各セルを順に走査
    basin_id = 1;
    for k = 1:length(sorted_lengths)
        % 変化ベクトルの長さが短い順に、そのセルの位置をrow行col列とする
        row = sorted_rows(k);
        col = sorted_cols(k);
        
        % すでに他のbasinに入っている場合、スキップ
        if visited(row, col) || isempty(average_vecs{row, col})
            continue;
        end
        
        % 走査対象の始め（現在のセル）
        targets = [row, col];
        
        l = 0;
        invalid = 0;
        while ~isempty(targets)
            current = targets(1, :);
            targets(1, :) = []; % 走査対象から削除
            
            % 走査済みの場合はスキップ
            if visited(current(1), current(2))
                continue;
            end
            
            % 現在のセルに basin_id を記録
            basins(current(1), current(2)) = basin_id;
            visited(current(1), current(2)) = true;
            
            % 現在のセルの周囲をチェック
            neighbors = [];
            if current(1) > 1
                neighbors = [neighbors; current(1)-1, current(2)];
            end
            if current(1) < m
                neighbors = [neighbors; current(1)+1, current(2)];
            end
            if current(2) > 1
                neighbors = [neighbors; current(1), current(2)-1];
            end
            if current(2) < n
                neighbors = [neighbors; current(1), current(2)+1];
            end
            
            % 周囲のセルを走査対象に追加
            targets_additional = [];
            for i = 1:size(neighbors, 1)
                neighbor_row = neighbors(i, 1);
                neighbor_col = neighbors(i, 2);
                
                % 隣接するセルがすでに他のbasinに属していたら、そのセルは走査対象外（ここ要考察）
                if visited(neighbor_row, neighbor_col)
                    continue;
                % 隣接するセルに変化ベクトルがない場合、
                elseif isempty(average_vecs{neighbor_row, neighbor_col})
                    % それが初回（l=0）であればcurrent周りのbasinは存在しないとする。(invalid=1)
                    if l == 0
                        invalid = 1;
                        break;
                    end
                    % それが初回でなければ（l>0）であればそのセルは走査対象外
                    continue;
                end
                
                % 隣接するセルにおける変化ベクトル（neighbor_vec）と、隣接するセルの中心とcurrentの中心を結ぶベクトル（neighbor_ideal_vec）
                current_center = cell_centers{current(1), current(2)};
                neighbor_center = cell_centers{neighbor_row, neighbor_col};
                neighbor_ideal_vec = current_center - neighbor_center;
                neighbor_vec = average_vecs{neighbor_row, neighbor_col};
                
                % ベクトルの内積を計算し、cosが0を超える場合、隣接するセルと同じbasinにcurrentを追加
                if dot(neighbor_ideal_vec, neighbor_vec) > 0
                    if l == 0
                        targets_additional = [targets_additional; neighbor_row, neighbor_col];
                    else
                        targets = [targets; neighbor_row, neighbor_col];
                    end
                % もしそうでなければ、初回であればcurrent周りのbasinは存在しないとする。(invalid=1)
                elseif l == 0
                    invalid = 1;
                    break;
                end
            end
            % current周りのbasinは存在しない場合、currentに記録していたbasinとvisitedを初期化する
            if l == 0
                if invalid == 1
                    basins(current(1), current(2)) = 0;
                    visited(current(1), current(2)) = false;
                else
                    targets = [targets; targets_additional];
                end
            end
            l = l+1;
        end
        basin_id = basin_id + 1; % 新しいbasin id
    end
end

% 条件分岐用の関数
function result = ifelse(condition, true_value, false_value)
  if condition
      result = true_value;
  else
      result = false_value;
  end
end