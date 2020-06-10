import sys

with open(sys.argv[1]) as f_input, open(sys.argv[2], 'w') as f_output:
    block = []

    for line in f_input:
        if line.startswith('>'):
            if block:
                f_output.write(''.join(block) + '\n')
                block = []
            f_output.write(line)
        else:
            block.append(line.strip())

    if block:
        f_output.write(''.join(block) + '\n')
