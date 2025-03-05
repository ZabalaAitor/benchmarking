import os
import sys

def filter_circles(input_dir, output_dir, split_read_col=5, use_split_read=True, split_read_threshold=2, remove_duplicates=True):
    """
    Filters and optionally deduplicates circles from BED files.

    Processes each BED file in the input directory, removing duplicate circles based on either split read count 
    (if `use_split_read` is True) or circle length (if `use_split_read` is False). Saves the results to the output directory.
    Only processes circles with split read count greater than or equal to `split_read_threshold` if `use_split_read` is True.
    If `remove_duplicates` is False, only filters based on split read count or length without removing duplicates.

    Parameters:
        input_dir (str): Directory with input BED files.
        output_dir (str): Directory for saving filtered BED files.
        split_read_col (int): Column index for split read count (1-based).
        use_split_read (bool): True to filter by split read count, False by circle length.
        split_read_threshold (int): Minimum split read count to include a circle if `use_split_read` is True.
        remove_duplicates (bool): True to remove duplicates, False to only filter based on split reads or length.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize variable to count duplicates
    total_duplicates = 0
    
    # Iterate over each bed file in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.bed'):
            input_path = os.path.join(input_dir, file_name)
            output_path = os.path.join(output_dir, file_name.replace('.bed', '.bed'))
            
            # Load all circles from the file
            circles = []
            with open(input_path, 'r') as input_file:
                for line in input_file:
                    fields = line.strip().split('\t')
                    if use_split_read:
                        split_read_count = int(fields[split_read_col - 1])
                        if split_read_count >= split_read_threshold:
                            circles.append((fields[0], int(fields[1]), int(fields[2]), split_read_count))
                    else:
                        circles.append((fields[0], int(fields[1]), int(fields[2])))
            
            # Remove duplicates based on overlapping and split-read counts if remove_duplicates is True
            if remove_duplicates:
                total_duplicates_in_file = 0
                unique_circles = []
                circles.sort(key=lambda x: (x[0], x[1], x[2]))  # Sort by chromosome, start, end
                
                i = 0
                while i < len(circles):
                    current_circle = circles[i]
                    overlapping_group = [current_circle]
                    
                    # Collect all overlapping circles
                    j = i + 1
                    while j < len(circles) and circles[j][0] == current_circle[0] and is_overlap(current_circle, circles[j]):
                        overlapping_group.append(circles[j])
                        j += 1
                    
                    if use_split_read:
                        # Keep the circle with the highest split-read count in the overlapping group
                        best_circle = max(overlapping_group, key=lambda x: x[3])
                    else:
                        # Keep the circle with the maximum length (end - start)
                        best_circle = max(overlapping_group, key=lambda x: x[2] - x[1])

                    unique_circles.append(best_circle)
                    
                    # Count duplicates (all but the best circle in the group)
                    total_duplicates_in_file += len(overlapping_group) - 1
                    total_duplicates += len(overlapping_group) - 1
                    
                    # Move to the next non-overlapping circle
                    i = j
            else:
                unique_circles = circles  # If not removing duplicates, all circles are considered unique
            
            # Write unique circles to output file
            with open(output_path, 'w') as output_file:
                for circle in unique_circles:
                    if use_split_read:
                        output_file.write(f"{circle[0]}\t{circle[1]}\t{circle[2]}\t{circle[3]}\n")
                    else:
                        output_file.write(f"{circle[0]}\t{circle[1]}\t{circle[2]}\n")
            
            if remove_duplicates:
                print(f"{total_duplicates_in_file} duplicates removed from {file_name}")
    
    if remove_duplicates:
        print(f"Total duplicates removed from all files: {total_duplicates}")

def is_overlap(circle1, circle2):
    """
    Check if circle1 overlaps with circle2 based on their start and end positions.
    """
    chr1, start1, end1 = circle1[:3]
    chr2, start2, end2 = circle2[:3]
    
    # Overlap occurs if the start of one circle is within the other circle's range
    return chr1 == chr2 and (start1 <= start2 <= end1 or start1 <= end2 <= end1 or (start1 <= start2 and end1 >= end2))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python filter_circles.py <input_directory> <output_directory> [split_read_col=5] [use_split_read=True] [split_read_threshold=2] [remove_duplicates=True]")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    
    # Set default values
    split_read_col = 5
    use_split_read = True
    split_read_threshold = 2
    remove_duplicates = True
    
    # Parse additional arguments
    for arg in sys.argv[3:]:
        if arg.startswith("split_read_col="):
            split_read_col = int(arg.split("=")[1])
        elif arg.startswith("use_split_read="):
            use_split_read = arg.split("=")[1].lower() == 'true'
        elif arg.startswith("split_read_threshold="):
            split_read_threshold = int(arg.split("=")[1])
        elif arg.startswith("remove_duplicates="):
            remove_duplicates = arg.split("=")[1].lower() == 'true'

    filter_circles(input_directory, output_directory, split_read_col, use_split_read, split_read_threshold, remove_duplicates)
