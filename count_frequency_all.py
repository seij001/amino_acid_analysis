from collections import defaultdict
from pathlib import Path

def read_and_count_continuous_sequences_from_directory(input_directory, output_file):
    """Read all files from the directory, count continuous sequences of length 1 to 6, and save the merged results to a text file."""
    sequence_counts = defaultdict(int)
    
    # 遍历目录中的所有文件
    input_directory = Path(input_directory)
    for file_path in input_directory.glob("*.txt"):  # 假设文件后缀是 .txt
        print(f"Processing file: {file_path.name}")
        
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        # 提取氨基酸序列及Residue Number
        amino_acid_data = []
        for line in lines[2:]:  # 跳过前两行标题
            parts = line.split()
            
            # 确保格式正确，过滤掉异常行
            if len(parts) < 4 or not parts[2].isdigit():
                print(f"Skipping line due to incorrect format in {file_path.name}: {line.strip()}")
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

    # 将结果保存到合并后的输出文件
    with open(output_file, 'w') as out_file:
        for seq, count in sequence_counts.items():
            out_file.write(f"{seq} {count}\n")
    
    print(f"Merged results saved to: {output_file}")

# 指定文件夹路径和输出文件路径
input_directory_path = "D:/abroad/Duke/CBB520/assginment3/amino_acid_analysis/extracted"
output_file_path = "D:/abroad/Duke/CBB520/assginment3/amino_acid_analysis/merged_sequence_counts.txt"

# 运行函数
read_and_count_continuous_sequences_from_directory(input_directory_path, output_file_path)
