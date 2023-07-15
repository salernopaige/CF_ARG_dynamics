import os

# Specify the root directory
root_dir = '/dartfs/rc/lab/R/RossB/OtherLabsPublished/Li2021/ReadData'

# Get a list of all subdirectories in the root directory
subdirs = [d for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))]

# Iterate over each subdirectory
for subdir in subdirs:
    # Get the path to the subdirectory
    subdir_path = os.path.join(root_dir, subdir)
    
    # Get a list of all files in the subdirectory
    files = os.listdir(subdir_path)
    # Iterate over each file in the subdirectory
    for filename in files:
        # Get the file extension
        file_ext = os.path.splitext(filename)[1]
        # Check if the file ends with "1.fastq.gz" or "2.fastq.gz"
        if file_ext == '.gz' and (filename.endswith('_1.fastq.gz') or filename.endswith('_2.fastq.gz')):
            # Generate the new filename
            new_filename = filename.replace('_1.fastq.gz', '_R1.fastq.gz').replace('_2.fastq.gz', '_R2.fastq.gz')
            # Rename the file
            old_file_path = os.path.join(subdir_path, filename)
            new_file_path = os.path.join(subdir_path, new_filename)
            os.rename(old_file_path, new_file_path)
            print(f'Renamed file: {filename} -> {new_filename}')
