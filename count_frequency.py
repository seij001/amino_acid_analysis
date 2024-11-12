from collections import defaultdict
from pathlib import Path

def read_and_count_continuous_sequences(file_path, output_file):
    """Read the file, sort by Residue Number, count continuous sequences of length 1 to 6, and save the results to a text file."""
    sequence_counts = defaultdict(int)
    
    # 读取文件并解析行
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # 提取氨基酸序列及Residue Number
    amino_acid_data = []
    for line in lines[2:]:  # 跳过前两行标题
        parts = line.split()
        
        # 确保格式正确，过滤掉异常行
        if len(parts) < 4 or not parts[2].isdigit():
            print(f"Skipping line due to incorrect format: {line.strip()}")
            continue
        
        residue = parts[0]  # 氨基酸类型
        residue_number = int(parts[2])  # Residue Number
        amino_acid_data.append((residue, residue_number))
    
    # 按 Residue Number 排序
    amino_acid_data.sort(key=lambda x: x[1])

    # 提取排序后的序列并判断连续性
    amino_acid_sequence = []
    previous_residue_number = None
    for residue, residue_number in amino_acid_data:
        if previous_residue_number is None or residue_number == previous_residue_number + 1:
            amino_acid_sequence.append((residue, residue_number))
        else:
            amino_acid_sequence.append(None)  # 插入None来中断不连续的序列
        previous_residue_number = residue_number

    # 统计连续序列的频次
    for i in range(len(amino_acid_sequence)):
        if amino_acid_sequence[i] is None:
            continue
        for length in range(1, 7):  # 生成长度从1到6的序列
            if i + length <= len(amino_acid_sequence) and all(amino_acid_sequence[j] is not None for j in range(i, i + length)):
                seq = '-'.join(residue for residue, _ in amino_acid_sequence[i:i+length])
                sequence_counts[seq] += 1

    # 将结果保存到文本文件
    with open(output_file, 'w') as out_file:
        for seq, count in sequence_counts.items():
            out_file.write(f"{seq} {count}\n")

    print(f"Results saved to: {output_file}")

# 指定文件路径
input_file_path = Path("D:/abroad/Duke/CBB520/assginment3/amino_acid_analysis/extracted/test.txt")
output_file_path = Path("D:/abroad/Duke/CBB520/assginment3/amino_acid_analysis/continuous_sequence_counts_test.txt")

# 运行函数
read_and_count_continuous_sequences(input_file_path, output_file_path)
