#
#  Copyright (c) 2023 Institute of Communication Networks (ComNets),
#                     Hamburg University of Technology (TUHH),
#                     https://www.tuhh.de/comnets
#  Copyright (c) 2023 Leonard Fisser <leonard.fisser@tuhh.de>
# 
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 
# Note: The code looks very over-complicated, but this is due to the fact that it was originally designed to work on GPU CUDA cores (which I think it sitll does work)
# Not all programming features are available on CUDA cores and therefore most operations are broken down into basic operations such as for loops and if statements.
# 给定的代码执行二进制矩阵操作，重点在于将矩阵转换为其行简化阶梯形（Reduced Row Echelon Form, RREF），并处理特定行的这种转换

# 关键部分的详细解析
# 1 前向消元 
#   主元寻找和消除：
#           对每列，代码识别具有主元的第一行（前导1）。
#           然后通过减去（或二进制加法）主元行来消除下方的条目。
# 2 主元排序
#   冒泡排序：
#       代码根据主元位置对行进行排序，确保具有较高前导项的行位于具有较低前导项的行之上。
#       这种排序确保矩阵遵循 RREF 结构。
# 3 后向消元
#   主元上方消元：
#       从最后一行开始向上移动，代码消除主元上方的条目，确保每个主元在其列中是唯一的非零条目。
# 4 零行移动到底部
#    零行识别和移动：
#       代码识别全零条目的行。
#       将这些行移动到底部，保持 RREF 结构。

# 这个函数使用高斯消元法将三维矩阵 A 转换为其行简化阶梯形（RREF）。函数操作于由 rx_id 标识的特定矩阵切片。
# A: 输入的三维矩阵。
# rx_id: 要转换的矩阵切片的索引。
# index: 可选参数，默认值为 1，指定第三维索引。


# 清除内容的时候,传入 状态A 和 接收节点号
function bin_mat_rref!(A, rx_id, index=1)
    """Performs Gaussian elimination to bring A(:,:,index) into Reduced-Row-Echolon-Form. 行阶梯形矩阵（Row-Echelon Form）"""
    num_colons = size(A, 2)
    # Forward pass 这一步逐行识别和消除主元（leading entries）以形成上三角矩阵
    curr_find = 0
    for col_id in 1:num_colons # for every coloumn, find the first row which has a pivot
        for row_id_search = curr_find+1:num_colons # as A is in upper diagonal form, we may only select a new pivot row after the current one
            if @inbounds A[col_id, row_id_search, rx_id, index] == 1
                curr_find = row_id_search
                break
            end
        end
        if curr_find > 0
            for row_id in curr_find+1:num_colons # iterate over all rows which are below the current pivot row
                if @inbounds A[col_id, row_id, rx_id, index] == 1 # if we have a 1 in the coloum, we need to substract the pivot row
                    for k in 1:num_colons
                        @inbounds A[k, row_id, rx_id, index] = (A[k, row_id, rx_id, index] + A[k, curr_find, rx_id, index]) % 2
                    end
                end
            end
        end
    end
    # pivot sorting 确保主元按特定顺序排列，将具有高主元的行放在具有低主元的行之上。
    # check all coloumns
    for col_id in 1:num_colons
        # check all rows
        for row_id in 1:num_colons
            # if the current row has a a one at the current coloumn
            @inbounds if A[col_id, row_id, rx_id, index] == 1
                replace_row = 0
                # check all rows above the current row, and check if we need to change positions
                for row_id_search in 1:row_id
                    # check every coloum entry
                    for col_id_search in 1:num_colons
                        # if we have a higher leading entry (or entries), this is the position we need to change with
                        @inbounds if A[col_id_search, row_id, rx_id, index] > A[col_id_search, row_id_search, rx_id, index]
                            replace_row = row_id_search
                            break
                        end
                        # if the current coloum had a entry which we did not have, we can skip the search to the next row (e.g. exit this col_id search)
                        @inbounds if A[col_id_search, row_id_search, rx_id, index] == 1
                            break
                        end
                    end
                    # If we found a row, we copy and move to sort the next row
                    if replace_row != 0
                        for col_id_copy in 1:num_colons
                            @inbounds element_copy = A[col_id_copy, replace_row, rx_id, index]
                            @inbounds A[col_id_copy, replace_row, rx_id, index] = A[col_id_copy, row_id, rx_id, index]
                            @inbounds A[col_id_copy, row_id, rx_id, index] = element_copy
                        end
                        break
                    end
                end
            end
        end
    end

    # backward pass ：消除主元上方的条目以形成 RREF。
    for row_id in num_colons:-1:1
        # check if the current row contains a leading 1 at or after identiy position
        leading_one_position = -1
        for col_id in row_id:num_colons
            @inbounds if A[col_id, row_id, rx_id, index] == 1
                leading_one_position = col_id
                break
            end
        end
        if leading_one_position > 0
            # check remaining rows for a pivot at the same position
            for to_reduce_row in row_id-1:-1:1
                # eliminate pivot by adding the current row
                @inbounds if A[leading_one_position, to_reduce_row, rx_id, index] == 1
                    for k in to_reduce_row:num_colons
                        @inbounds A[k, to_reduce_row, rx_id, index] = (A[k, to_reduce_row, rx_id, index] + A[k, row_id, rx_id, index]) % 2
                    end
                end
            end
        end
    end
    # put all zero lines to bottom 将全零行移动到底部。
    for row_id in 1:num_colons
        # check if current row is only zeros
        zero_row = true
        for col_id in 1:num_colons
            @inbounds if A[col_id, row_id, rx_id, index] > 0
                zero_row = false
                break
            end
        end
        if zero_row
            # search the remaining rows for a row with entries
            for row_id_search in row_id+1:num_colons
                zero_row_search = true
                for col_id_search in 1:num_colons
                    if @inbounds A[col_id_search, row_id_search, rx_id, index] > 0
                        zero_row_search = false
                        break
                    end
                end
                if !zero_row_search # found non zero row below a zero row, we need to copy
                    for col_id_copy in 1:num_colons
                        @inbounds A[col_id_copy, row_id, rx_id, index] = A[col_id_copy, row_id_search, rx_id, index]
                        @inbounds A[col_id_copy, row_id_search, rx_id, index] = 0
                    end
                    break
                end
            end
        end
    end
    return A
end

# 这个函数针对添加到三维矩阵 A 中的特定行（new_row）执行高斯消元，重点在于将这行新行纳入现有的 RREF。
# 参数：
# A: 输入的三维矩阵。
# rx_id: 要转换的矩阵切片的索引。
# new_row: 新行的行索引。
# index: 可选参数，默认值为 1，指定第三维索引。

# 按列求和为0的1个列号为  new_row
function bin_mat_rref_line!(A, rx_id, new_row, index=1)
    """Performs Gaussian elimination to bring A(:,:,index) into Reduced-Row-Echolon-Form. (枢轴是所在列的唯一非0元素)
       Works only on eliminating information introduced by the new_row index row.
    """
    # 高斯消元法简化，新列引入信息的（state第一个求和为0的列，加上payload，之后，这个列索引用于这个函数的一个参数）
    num_colons = size(A, 2)
    # Forward pass 识别新行的主元并消除其下方的条目
    new_row_pivot = 0
    found_first_unhandled_pivot = false
    for col_id in 1:num_colons # for every coloumn, find the first row which has a pivot 对每个行，找到有一个枢纽的第一列，
        curr_find = 0
        if @inbounds A[col_id, new_row, rx_id, index] == 1 # if we have a 1 in the coloum, we need to substract the pivot row 如果行中有1，需要去掉这一枢纽列
            found_necessary_pivot = false

            # 如1 发1 给了3， 第一个求和不为0的列是，第二列new_row=2，加上负载 1 0 0
            # 此时进来的节点3是
            # 0 1 0
            # 0 0 0
            # 1 0 0

            #节点3 上面满足==1时候，col_id=1 即发出内容1；已知 new_row=2 
            for row_id_search = curr_find+1:new_row-1 
                # as A is in upper diagonal form, we may only select a new pivot row after the current one
                # 对于A的上对角形式，对当前的1，在其之后选择一个新的枢纽
                # curr_find+1=1  new_row-1=1

                # 就等于1

                correct_pivot = true
                @inbounds if A[col_id, row_id_search, rx_id, index] == 1 # candidate row
                    # 如果A[1,1,3]是1 （这里不满足）

                    for c_r = 1:col_id-1 # check for leading zeros
                        @inbounds if A[c_r, row_id_search, rx_id, index] == 1 # bad candidate
                            correct_pivot = false
                            break
                        end
                    end
                else
                    continue
                end


                if correct_pivot # good candidate  满足
                    curr_find = row_id_search
                    # 这里 curr_find = 1 
                    found_necessary_pivot = true
                    for k in col_id:num_colons
                        @inbounds A[k, new_row, rx_id, index] = (A[k, new_row, rx_id, index] + A[k, curr_find, rx_id, index]) % 2
                    end
                    break # we handled this pivot, nothing to be done anymore
                end
            end


            if !found_first_unhandled_pivot
                if !found_necessary_pivot # We found a pivot which we weren't able to resolve, so this will be the leading element in the new line, we save it
                    found_first_unhandled_pivot = true
                    new_row_pivot = col_id
                end
            end


        end


    end

    # pivot sorting 确保新行按主元正确排序。
    # 1) find the pivot element of the new row
    # We logged this already during the forward pass step: new_row_pivot
    # Row was elimiated by forward pass, nothing to be done anymore
    if new_row_pivot == 0
        return nothing
    end

    # 2) find row in which the new row needs to be placed
    # & 3) bubble sort the row to the new location
    new_row_new_position = new_row
    row_pivot = 0
    element_copy = 0
    for row_id in new_row-1:-1:1
        row_pivot = 0
        for col_id in 1:num_colons
            @inbounds if A[col_id, row_id, rx_id, index] == 1
                row_pivot = col_id
                break
            end
        end
        if row_pivot > new_row_pivot # switch rows
            for col_id_copy in 1:num_colons
                element_copy = A[col_id_copy, row_id, rx_id, index]
                @inbounds A[col_id_copy, row_id, rx_id, index] = A[col_id_copy, new_row_new_position, rx_id, index]
                @inbounds A[col_id_copy, new_row_new_position, rx_id, index] = element_copy
            end
            new_row_new_position = row_id
        else # we are done
            break
        end
    end

    # backward pass 消除主元上方的条目。
    for row_id in new_row:-1:1
        # check if the current row contains a leading 1 at or after identity position
        row_pivot = 0
        for col_id in row_id:num_colons
            if A[col_id, row_id, rx_id, index] == 1
                row_pivot = col_id
                break
            end
        end
        if row_pivot > 0
            # check remaining rows for a pivot at the same position
            for to_reduce_row in row_id-1:-1:1
                # eliminate pivot by adding the current row
                if A[row_pivot, to_reduce_row, rx_id, index] == 1
                    for k in row_pivot:num_colons
                        @inbounds A[k, to_reduce_row, rx_id, index] = (A[k, to_reduce_row, rx_id, index] + A[k, row_id, rx_id, index]) % 2
                    end
                end
            end
        end
    end
    # put all zero lines to bottom 包括新行在内的全零行移动到底部。
    for row_id in 1:new_row-1
        # check if current row is only zeros
        zero_row = true
        for col_id in 1:num_colons
            @inbounds if A[col_id, row_id, rx_id, index] > 0
                zero_row = false
                break
            end
        end
        if zero_row
            # search the remaining rows for a row with entries
            for row_id_search in row_id+1:new_row
                zero_row_search = true
                for col_id_search in 1:num_colons
                    if @inbounds A[col_id_search, row_id_search, rx_id, index] > 0
                        zero_row_search = false
                        break
                    end
                end
                if !zero_row_search # found non zero row below a zero row, we need to copy
                    for col_id_copy in 1:num_colons
                        @inbounds A[col_id_copy, row_id, rx_id, index] = A[col_id_copy, row_id_search, rx_id, index]
                        @inbounds A[col_id_copy, row_id_search, rx_id, index] = 0
                    end
                    break
                end
            end
        end
    end

    return A
end
