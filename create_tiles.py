import os
import sys

""" creates .top and .conf files for a tile binding simulation """
def get_top_file(top_file, num_tiles):
    """ lengths of all tile properties must be uniform """
    
    with open(top_file, 'r') as f:
        lines = f.readlines()[1:]

        # bases
        strand_1 = [x.split()[1] for x in lines if x[0] == '1']
        strand_2 = [x.split()[1] for x in lines if x[0] == '2']

    assert len(strand_1) == len(strand_2)
    strands = [strand_1, strand_2]
    strand_length = len(strand_2)

    output = f'{num_tiles*2*strand_length} {num_tiles*2}\n'

    for tile_num in range(num_tiles):
        for i, strand in enumerate(strands):
            for j, base in enumerate(strand):

                strand_index = i + 1 + tile_num * 2
                base_index = j + i * strand_length + tile_num * 2 * strand_length

                next_base = base_index + 1 if j < strand_length - 1 else -1
                prev_base = base_index - 1 if j != 0 else -1

                output += f'{strand_index} {base} {prev_base} {next_base}\n'


    return output

def get_conf_file(conf_file, num_tiles, strand_length):
    """ returns the contents of the conf file with num_strands """
    
    with open(conf_file, 'r') as f:
        lines = f.readlines()

    # box_size = max([float(x) for x in lines[1].split('=')[-1].split()])
    box_size = 1000
        
    output = f't = 0\nb = {box_size}.00 {box_size}.00 {box_size}.00\nE = 0.00 0.00 0.00\n'

    positions = lines[3:]

    assert len(positions) == 2*strand_length

    strand_1 = get_normalized_pos([[float(x) for x in line.split()] for line in positions[:strand_length]])
    strand_2 = get_normalized_pos([[float(x) for x in line.split()] for line in positions[strand_length:]])

    strand_confs = [strand_1, strand_2]

    
    # 8 works for strand_length of 40, 3.2 works for strand_length of 16
    centers = generate_xyz(box_size, x_spacing=strand_length/5, y_spacing=100, z_spacing=6)

    for tile_num in range(num_tiles):
        for i, strand in enumerate(strand_confs):
            center = next(centers)
            print(center)
            for j, base in enumerate(strand):
                position = base[:]
                position[0] += center[0]
                position[1] += center[1]
                position[2] += center[2]
                output += get_base_conf_str(position) +'\n'

    return output

def get_normalized_pos(lines):
    """ lines is a list of list of floats """
    centroid = [0,0,0]
    for line in lines:
        centroid[0] += line[0]
        centroid[1] += line[1]
        centroid[2] += line[2]

    centroid[0] /= len(lines)
    centroid[1] /= len(lines)
    centroid[2] /= len(lines)

    for line in lines: # lists are mutable
        line[0] -= centroid[0]
        line[1] -= centroid[1]
        line[2] -= centroid[2]

    print(centroid)
    print(centroid)
    return lines

def get_base_conf_str(position):
    """ position is list of floats """
    return ' '.join(['{:.2f}'.format(x) for x in position])

def generate_xyz(box_size, x_spacing, y_spacing, z_spacing):
    is_offset = True
    box_size = box_size/4
    for i in range(int(box_size/x_spacing)-1):
        for j in range(int(box_size/y_spacing)-1):
            for k in range(int(box_size/z_spacing)-1):

                center = [i*x_spacing + box_size/2, j*y_spacing + box_size/2, k*z_spacing + box_size/2]

                if is_offset:
                    center[0] += x_spacing/2
                
                is_offset = not is_offset
                    
                yield center



if __name__ == '__main__':
    """ usage:
            python create_files.py name num_tiles output_dir
        Example:
            python create_files.py 4 100 .
    """
    output_dir = sys.argv[3]
    name = os.path.join(output_dir, sys.argv[1])
    top_file = f'{name}.top'
    conf_file = f'{name}.conf'

    num_tiles = int(sys.argv[2])
    
    box_size = 200

    out_name = f'{sys.argv[1]}_{num_tiles}'
    out_top = f'{out_name}.top'
    out_conf = f'{out_name}.conf'

    with open(out_top, 'w') as f:
        f.write(get_top_file(top_file, num_tiles))

    with open(out_conf, 'w') as f:
        f.write(get_conf_file(conf_file, num_tiles, int(sys.argv[1])*4))

    print(top_file, conf_file)